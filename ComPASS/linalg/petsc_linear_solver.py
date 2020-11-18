#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .__init__ import *
from .solver import *
from .exceptions import IterativeSolverFailure, DirectSolverFailure


class PetscLinearSystem:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self, simulation):
        """
        :param simulation: an initialised simulation object to get the data from
        """
        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()

        self.lsbuilder = simulation.LinearSystemBuilder()
        (sizes, d_nnz, o_nnz) = self.lsbuilder.get_non_zeros()
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))
        self.x.createMPI((n_rowl, n_rowg))
        self.RHS.createMPI((n_rowl, n_rowg))

    def dump_ascii(self, basename, comm=PETSc.COMM_WORLD):
        """
        Writes the linear system (Matrix, solution and RHS) in three different files in ASCII format

        :param basename: common part of the file names
        :comm: MPI communicator
        """
        makeviewer = PETSc.Viewer().createASCII

        def dump_item(item, name):
            viewer = makeviewer(basename + name + ".dat", "w", comm)
            item.view(viewer)
            viewer.destroy

        dump_item(self.A, "A")
        dump_item(self.RHS, "RHS")
        dump_item(self.x, "x")

    def dump_binary(self, basename="", comm=PETSc.COMM_WORLD):
        """
                Writes the linear system (Matrix, solution and RHS) in three different files in binary format

                :param basename: common part of the file names
                :comm: MPI communicator
                """
        makeviewer = PETSc.Viewer().createBinary

        def dump_item(item, name):
            viewer = makeviewer(basename + name + ".dat", "w", comm)
            item.view(viewer)
            viewer.destroy

        dump_item(self.A, "A")
        dump_item(self.RHS, "RHS")
        dump_item(self.x, "x")

    def check_residual_norm(self):

        """
        Displays the residual norm (1-Norm, 2-Norm and infinity norm) for convergence check
        """

        y = self.RHS.duplicate()
        self.RHS.copy(y)  # y = b
        y.scale(-1.0)  # y = -b
        self.A.multAdd(self.x, y, y)  # y = Ax-b
        mpi.master_print("Linear solution check ||Ax-b||")
        norm = y.norm(PETSc.NormType.NORM_1)
        mpi.master_print("  1-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_2)
        mpi.master_print("  2-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_INFINITY)
        mpi.master_print("  Infinity norm", norm)

    def set_from_jacobian(self):

        self.lsbuilder.set_AMPI(self.A)
        self.lsbuilder.set_RHS(self.RHS)


class PetscIterativeSolver(IterativeSolver):

    """
    A structure that holds an iterative PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, settings, activate_cpramg=True, comm=PETSc.COMM_WORLD,
    ):

        self.activate_cpramg = activate_cpramg
        super().__init__(linear_system, settings)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        self.pc = self.ksp.getPC()
        if self.activate_cpramg:
            self.set_cpramg_pc(comm)
        else:
            self.activate_cpramg = False
            self.pc.setFactorLevels(1)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        self.tolerance, self.max_iterations, self.restart_size = self.settings[:]
        self.ksp.setFromOptions()

    tolerance = property(
        fget=lambda self: self.ksp.rtol,
        fset=lambda self, value: self.ksp.setTolerances(rtol=value),
        doc="Relative decrease in the residual norm required for convergence",
    )
    max_iterations = property(
        fget=lambda self: self.ksp.max_it,
        fset=lambda self, value: self.ksp.setTolerances(max_it=value),
        doc="Maximum number of iterations accepted before convergence failure",
    )
    restart_size = property(
        # restart_size is not available in petsc4py,
        # but can be set using the command line petsc option -ksp_gmres_restart restart
        fget=None,
        fset=lambda self, value: self.ksp.setGMRESRestart(value),
        doc="Number of iterations at which GMRES restarts",
    )

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        ksp_reason = self.ksp.getConvergedReason()
        nit = self.ksp.getIterationNumber()

        if ksp_reason < 0:
            self.number_of_unsuccessful_iterations += nit
            raise IterativeSolverFailure(
                "Petsc KSP object returned error code: %s" % ksp_reason, nit
            )
        else:
            self.number_of_successful_iterations += nit

        return self.linear_system.x, nit, ksp_reason

    def set_cpramg_pc(self, comm):

        # CPR-AMG is a multiplicative composite PC :
        # CPR-AMG = M2*(I - A*M1) + M1
        # Where M1 is an AMG procedure on the pressure unknowns,
        # set using the Petsc fieldsplit type.
        # and M2 = ILU(A)
        cpramg_pc = self.pc
        cpramg_pc.setType(PETSc.PC.Type.COMPOSITE)
        cpramg_pc.setCompositeType(PETSc.PC.CompositeType.MULTIPLICATIVE)
        cpramg_pc.addCompositePC(PETSc.PC.Type.FIELDSPLIT)
        cpramg_pc.addCompositePC(PETSc.PC.Type.BJACOBI)
        cpramg_pc.setUp()

        # We now define M1 as an additive fieldsplit pc
        # M1 = M11 + M12
        # with M11 an AMG v-cycle procedure on the pressure field
        # and M12 = 0 a Null PC which returns zero on the T/s field
        block_size = self.linear_system.lsbuilder.get_block_size()
        (sizes, d_nnz, o_nnz) = self.linear_system.lsbuilder.get_non_zeros()
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes
        n_wells = self.linear_system.lsbuilder.get_n_wells()
        non_well_nrowl = n_rowl - n_wells
        fs_pc = cpramg_pc.getCompositePC(0)
        p_indices = np.arange(non_well_nrowl, step=block_size, dtype="int32")
        p_indices = np.concatenate(
            (p_indices, np.arange(start=non_well_nrowl, stop=n_rowl, dtype="int32"))
        )


        rest_size = non_well_nrowl - (non_well_nrowl // block_size)
        rest_indices = np.zeros(rest_size, dtype="int32")
        for i in range(non_well_nrowl // block_size):
            for j in range(block_size - 1):
                k = i * (block_size - 1) + j
                rest_indices[k] = k + i + 1


        # The Index Set for the pressure field
        p_IS = PETSc.IS().createGeneral(p_indices, comm=comm)
        # The Index Set for the temperature and saturation field
        rest_IS = PETSc.IS().createGeneral(rest_indices, comm=comm)
        fs_pc.setFieldSplitIS(("pressure", p_IS), ("rest", rest_IS))
        fs_pc.setFieldSplitType(PETSc.PC.CompositeType.ADDITIVE)
        sub_ksp_list = fs_pc.getFieldSplitSubKSP()

        # Algebraic multigrid procedure on the pressure field
        pressure_ksp = sub_ksp_list[0]
        pressure_ksp.setType(PETSc.KSP.Type.PREONLY)
        pressure_pc = pressure_ksp.getPC()
        pressure_pc.setType(PETSc.PC.Type.HYPRE)
        PETSc.Options().setValue(
            "-sub_0_fieldsplit_pressure_pc_hypre_boomeramg_strong_threshold", 0.5
        )

        class NullPC(object):
            """ A PC-Python Context which returns zero. Used on the T/s field """

            def setUp(self, pc):
                pass

            def apply(self, pc, x, y):
                y.set(0.0)

        rest_ksp = sub_ksp_list[1]
        rest_ksp.setType(PETSc.KSP.Type.PREONLY)
        null_pc = rest_ksp.getPC()
        null_pc.setType(PETSc.PC.Type.PYTHON)
        null_pc.setPythonContext(NullPC())

    def __str__(self):

        cpramg_description = (
            "activated" if self.activate_cpramg == True else "not activated"
        )
        return f"{super().__str__()}\n   petsc4py new implementation\n   Settings : {self.settings}\n   CPR-AMG : {cpramg_description}"


class PetscDirectSolver(DirectSolver):
    """
    A structure that holds a direct PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, comm=PETSc.COMM_WORLD,
    ):

        super().__init__(linear_system)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        self.ksp.setType("preonly")
        self.ksp.getPC().setType("lu")
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        ksp_reason = self.ksp.getConvergedReason()

        if ksp_reason < 0:
            raise DirectSolverFailure(
                "Petsc KSP object returned error code: %s" % ksp_reason
            )

        return self.linear_system.x, None, ksp_reason

    def __str__(self):
        return f"{super().__str__()}\n   petsc4py new implementation"
