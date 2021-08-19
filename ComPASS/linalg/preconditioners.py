from .__init__ import PETSc, mpi, np


class NonePC(PETSc.PC):
    def __init__(self, linear_system, comm=PETSc.COMM_WORLD):
        self.create(comm)
        self.setOperators(linear_system.A, linear_system.A)
        self.setType(PETSc.PC.Type.NONE)


class PressureAMG:
    def __init__(self, A, pressure_IS, amg_type, comm=PETSc.COMM_WORLD):
        self.pressure_IS = pressure_IS
        self.amg_pc = PETSc.PC().create(comm=comm)
        self.A = A

        # PETSc offers several implementations for the AMG procedure
        # Following code block sets it up
        if amg_type not in ("gamg", "hypre"):
            mpi.master_print(
                f"Unknown AMG type in CPR-AMG preconditioner : {amg_type}; supported types : hypre, gamg\nUsing default hypre"
            )
            amg_type = "hypre"

        if amg_type == "hypre":
            try:
                self.amg_pc.setType(PETSc.PC.Type.HYPRE)
                PETSc.Options().setValue(  # petsc4py API doesn't allow us to set this parameter another way
                    "-pc_hypre_boomeramg_strong_threshold", 0.5,
                )
            except PETSc.Error:  # If PETSc is not compiled with HYPRE use default AMG from PETSc
                mpi.master_print(
                    "Hypre BoomerAMG is not available, using PETSc's GAMG procedure instead"
                )
                amg_type = "gamg"
        if amg_type == "gamg":
            self.amg_pc.setType(PETSc.PC.Type.GAMG)
            self.amg_pc.setGAMGType("agg")
        self.amg_pc.setFromOptions()

    def setUp(self, pc):
        # If Rp is the restriction matrix to the pressure unknowns
        # then pressure_mat = Ap = Rp*A*RpT
        # and amg_pc = AMG(Ap)
        pressure_mat = self.A.createSubMatrix(self.pressure_IS, self.pressure_IS)
        self.amg_pc.setOperators(pressure_mat, pressure_mat)
        self.amg_pc.setUp()

    def apply(self, pc, x, y):
        y.set(0.0)
        # Retrieve a view to the input and output vectors,
        # restricted to the pressure unknowns
        apply_subvec = x.getSubVector(self.pressure_IS)
        return_subvec = y.getSubVector(self.pressure_IS)
        self.amg_pc.apply(apply_subvec, return_subvec)
        y.restoreSubVector(self.pressure_IS, return_subvec)


class BlockJacobi(PETSc.PC):
    def __init__(self, linear_system, comm=PETSc.COMM_WORLD):
        self.create(comm)
        self.setOperators(linear_system.A, linear_system.A)
        self.setType(PETSc.PC.Type.BJACOBI)
        self.setFactorLevels(1)

    def __repr__(self):
        return "Block Jacobi performing ILU(1) on each block"


class CPRAMG(PETSc.PC):
    """ PETSc implementation of the CPR-AMG preconditioning procedure"""

    def __init__(self, linear_system, amg_type="hypre", comm=PETSc.COMM_WORLD):
        # CPR-AMG is a multiplicative composite PC :
        # CPR-AMG = M2*(I - A*M1) + M1
        # Where M1 is an AMG procedure on the pressure unknowns
        # and M2 = BlockJacobi(A)
        self.create(comm)
        self.setOperators(linear_system.A, linear_system.A)
        self.setType(PETSc.PC.Type.COMPOSITE)
        self.setCompositeType(PETSc.PC.CompositeType.MULTIPLICATIVE)
        PETSc_version = PETSc.Sys.getVersion()
        if PETSc_version[0] >= 3 and PETSc_version[1] > 14:
            ## This function's name changed in version 3.15
            addCompositePCType = self.addCompositePCType
        else:
            addCompositePCType = self.addCompositePC
        addCompositePCType(PETSc.PC.Type.PYTHON)
        addCompositePCType(PETSc.PC.Type.BJACOBI)
        self.setUp()

        # We now define M1 as a Python Type PC applying an AMG v-cycle procedure
        # on the pressure and well unknowns, and returning 0 for all other unknowns
        # We start by building the PETSc.IS which stores the global row indices
        # of the pressure unknowns owned by the current proc
        # The pressure unknown is always the first of each block
        lsbuilder = linear_system.lsbuilder
        block_size = lsbuilder.get_block_size()
        (sizes, d_nnz, o_nnz) = lsbuilder.get_non_zeros()
        n_rowl, n_rowg = sizes
        n_wells = lsbuilder.get_n_wells()
        non_well_nrowl = n_rowl - n_wells
        my_rowstart = lsbuilder.get_rowstart(mpi.proc_rank)
        p_ids = my_rowstart + np.arange(non_well_nrowl, step=block_size, dtype="int32")

        # We now append the trailing well indices if there are any
        # because we want wells to be included in the AMG procedure
        w_ids = my_rowstart + np.arange(non_well_nrowl, n_rowl, dtype="int32")
        p_ids = np.concatenate((p_ids, w_ids))

        # The actual IndexSet object
        p_IS = PETSc.IS().createGeneral(p_ids, comm=comm)

        # Following code block is just the syntax for setting
        # p_sub_pc = M1, e.g. the first sub PC of CPR-AMG
        p_sub_pc = self.getCompositePC(0)
        p_pc_context = PressureAMG(linear_system.A, p_IS, amg_type, comm=comm)
        p_sub_pc.setPythonContext(p_pc_context)
        self.setUp()

    def __repr__(self):
        amg_type = self.getCompositePC(0).getPythonContext().amg_pc.getType()
        return f"CPR-AMG preconditioner using {amg_type} AMG procedure on the pressure unknowns"
