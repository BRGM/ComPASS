import petsc4py
import sys

petsc4py.init(sys.argv)
from . import mpi
from petsc4py import PETSc
from collections import namedtuple

IterativeSolverSettings = namedtuple(
    "IterativeSolverSettings",
    ["tolerance", "max_iterations", "gmres_restart"],
    defaults=[None],
)


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


class LinearSolver:
    """
    A base structure to hold the common parameters of ComPASS linear solvers
    """

    def __init__(self, linear_system):
        """
        :param linear_system: the linear system structure
        """
        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.linear_system = linear_system


class PetscLinearSolver(LinearSolver):
    """
    A base structure to manage the common objects and methods of PETSc linear solvers
    """

    def __init__(self, linear_system, comm):

        """
        :param comm: MPI communicator
        """
        super().__init__(linear_system)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        reason = self.ksp.getConvergedReason()

        return reason

    def get_iteration_number(self):

        return self.ksp.getIterationNumber()


class PetscIterativeSolver(PetscLinearSolver):

    """
    A structure that holds an iterative PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, settings, activate_cpramg=True, comm=PETSc.COMM_WORLD,
    ):
        """
        :param settings: An IterativeSolverSettings object containing the wanted parameters for iterative solving
        """

        super().__init__(linear_system, comm)
        self.activate_cpramg = activate_cpramg
        self.activate_direct_solver = False
        self.ksp.getPC().setFactorLevels(1)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        self.tolerance, self.max_iterations, self.gmres_restart = settings[:]

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
    gmres_restart = property(
        fget=None,  # gmres_restart is not available in petsc4py  :(
        fset=lambda self, value: self.ksp.setGMRESRestart(value),
        doc="Number of iterations at which GMRES restarts",
    )


class PetscDirectSolver(PetscLinearSolver):
    """
    A structure that holds a direct PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, comm=PETSc.COMM_WORLD,
    ):

        super().__init__(linear_system, comm)
        self.activate_direct_solver = True
        self.ksp.setType("preonly")
        self.ksp.getPC().setType("lu")
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
