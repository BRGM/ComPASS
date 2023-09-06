import time, sys, yaml
import numpy as np
import petsc4py

petsc4py.init(sys.argv)
from petsc4py import PETSc

try:
    casename = sys.argv[1]
except IndexError:
    print("Syntax :\npython3 load_and_solve.py <path/to/a/dumped/linear/system>")
    exit()


class PressureAMG(object):
    """The PC object used by CPR-AMG to apply an AMG procedure on the pressure unknowns"""

    def __init__(self, A, pressure_IS, amg_type, comm=PETSc.COMM_WORLD):
        self.pressure_IS = pressure_IS
        pressure_mat = A.createSubMatrix(pressure_IS, pressure_IS)
        self.amg_pc = PETSc.PC().create(comm=comm)
        self.amg_pc.setOperators(pressure_mat, pressure_mat)

        # PETSc offers several implementations for the AMG procedure
        if amg_type not in ("gamg", "hypre"):
            mpi.master_print(
                f"Unknown AMG type in CPR-AMG preconditioner : {amg_type}; supported types : hypre, gamg\nUsing default hypre"
            )
            amg_type = "hypre"

        if amg_type == "hypre":
            try:
                self.amg_pc.setType(PETSc.PC.Type.HYPRE)
                PETSc.Options().setValue(
                    "-pc_hypre_boomeramg_strong_threshold",
                    0.5,
                )
            except PETSc.Error:  # If HYPRE is not available use default AMG from PETSc
                mpi.master_print(
                    "Hypre BoomerAMG is not available, using PETSc's GAMG procedure instead"
                )
        elif amg_type == "gamg":
            self.amg_pc.setType(PETSc.PC.Type.GAMG)
            self.amg_pc.setGAMGType("agg")

        self.amg_pc.setFromOptions()
        self.amg_pc.setUp()

    def setUp(self, pc):
        pass

    def apply(self, pc, x, y):
        y.set(0.0)
        apply_subvec = x.getSubVector(self.pressure_IS)
        return_subvec = y.getSubVector(self.pressure_IS)
        self.amg_pc.apply(apply_subvec, return_subvec)
        # x.view()
        y.restoreSubVector(self.pressure_IS, return_subvec)


class CPRAMG(PETSc.PC):
    """PETSc implementation of the CPR-AMG preconditioning procedure"""

    def __init__(self, A, part_data, amg_type="hypre", comm=PETSc.COMM_WORLD):
        # CPR-AMG is a multiplicative composite PC :
        # CPR-AMG = M2*(I - A*M1) + M1
        # Where M1 is an AMG procedure on the pressure unknowns,
        # set using the Petsc fieldsplit type.
        # and M2 = ILU(A)

        self.create(comm)
        self.setOperators(A, A)
        self.setType(PETSc.PC.Type.COMPOSITE)
        self.setCompositeType(PETSc.PC.CompositeType.MULTIPLICATIVE)
        if PETSc.Sys.getVersion() >= (3, 15):
            ## This function's name changed in version 3.15
            addCompositePCType = self.addCompositePCType
        else:
            addCompositePCType = self.addCompositePC
        addCompositePCType(PETSc.PC.Type.PYTHON)
        addCompositePCType(PETSc.PC.Type.BJACOBI)
        self.setUp()

        # We now define M1 as a Python Type PC applying an AMG v-cycle procedure
        # on the pressure and well unknowns, and returning 0 for every other unknowns
        # We start by building the PETSc.IS which stores the global row indices
        # of the pressure unknowns of the current proc
        (sizes, my_rowstart, block_size, n_wells) = part_data
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        non_well_nrowl = n_rowl - n_wells
        p_indices = my_rowstart + np.arange(
            non_well_nrowl, step=block_size, dtype="int32"
        )
        # Adding the trailing well equations if there are any (we want wells to be
        # included in the AMG procedure)
        p_indices = np.concatenate(
            (
                p_indices,
                my_rowstart
                + np.arange(start=non_well_nrowl, stop=n_rowl, dtype="int32"),
            )
        )

        p_IS = PETSc.IS().createGeneral(p_indices, comm=comm)

        self.amg_type = amg_type
        amg_pc = self.getCompositePC(0)
        amg_pc.setType(PETSc.PC.Type.PYTHON)
        amg_pc.setPythonContext(PressureAMG(A, p_IS, self.amg_type, comm=comm))
        self.setUp()

    def __repr__(self):
        return f"CPR-AMG preconditioner using {self.amg_type} AMG procedure on the pressure unknowns"


comm = PETSc.COMM_WORLD
size = comm.getSize()
rank = comm.getRank()

# Loading the partitioning data for the input linear system. This file is
# dumped with the linear solver when performing a binary_dump with a linear solver
# of the new implementation
with open(casename + "part_data.yaml", "r") as f:
    all_part = yaml.safe_load(f)

my_part = all_part[f"proc {rank}"]
sizes = my_part["local_number_of_rows"], all_part["total_number_of_rows"]
my_rowstart = my_part["global_index_of_first_row"]
block_size = all_part["block_size"]
n_wells = my_part["local_number_of_wells"]
part_data = (sizes, my_rowstart, block_size, n_wells)


# Loading the linear systems
# A: Linear operator
# b: Right hand side
# x: First guess and solution vector
viewer = PETSc.Viewer().createBinary(casename + "A.dat", "r", comm=comm)
A = PETSc.Mat().createAIJ((sizes, sizes), comm=comm)
A.load(viewer)
print(f"[{rank}] A loaded, sizes =", sizes)
x = A.createVecRight()
x.set(0.0)
viewer = PETSc.Viewer().createBinary(casename + "RHS.dat", "r", comm=comm)
b = PETSc.Vec().createMPI(sizes, comm=comm)
b.load(viewer)
viewer.destroy()

# KSP object definition
ksp = PETSc.KSP().create(comm=comm)
ksp.setOperators(A, A)
ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
ksp.setTolerances(rtol=1e-6, atol=1.0e-20, max_it=1000)
ksp.setPC(CPRAMG(A, part_data, amg_type="hypre"))
PETSc.Options().setValue("-ksp_monitor", None)
ksp.setFromOptions()
ksp.setUp()

# Solve
start = time.time()
ksp.solve(b, x)
end = time.time()

PETSc.Sys.Print("*" * 30)
PETSc.Sys.Print(
    f"Solve time : {end-start}\nTime per iteration : {(end-start)/ksp.getIterationNumber()}"
)
