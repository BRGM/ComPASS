import petsc4py
import sys

petsc4py.init(sys.argv)
from . import mpi
from petsc4py import PETSc


class PetscLinearSystem:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self):

        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()

    def create(self, nonzeros):

        """
        Allocate memory for the PETSc structures using the 'create' functions
        """

        (sizes, d_nnz, o_nnz) = nonzeros
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))
        self.x.createMPI((n_rowl, n_rowg))
        self.RHS.createMPI((n_rowl, n_rowg))

    def _dump(self, basename, makeviewer, comm):
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


class LinearSolver:

    """
    A structure that holds the PETSc KSP Object to solve the linear system
    """

    def __init__(
        self,
        simulation,
        tol=1e-6,
        maxit=150,
        restart=None,
        activate_cpramg=True,
        activate_direct_solver=False,
        comm=PETSc.COMM_WORLD,
    ):
        """
        :param simulation: an initialised simulation object to get the data from
        :param tol: relative tolerance (for iterative solvers).
        :param maxit: maximum number of iterations (for iterative solvers).
        :param restart: number of iterations before a restart (for gmres like iterative solvers).
        :param activate_cpramg: CPR-AMG preconditioner activator
        :param activate_direct_solver: direct solver activator
        :param comm: MPI communicator
        """
        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.last_residual_history = []
        self.activate_cpramg = activate_cpramg
        self.activate_direct_solver = activate_direct_solver
        self.rtol = tol
        self.maxit = maxit
        self.restart = restart

        self.linear_system = PetscLinearSystem()
        self.lsbuilder = simulation.LinearSystemBuilder()
        self.linear_system.create(self.lsbuilder.get_non_zeros())

        self.ksp = PETSc.KSP().create(comm=comm)
        self.set_parameters(self.rtol, self.maxit, self.restart)
        self.last_residual_history = []
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        if self.activate_direct_solver:
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.getPC().setFactorLevels(1)
        self.ksp.setFromOptions()

        def monitor(ksp, it, rnorm):
            self.last_residual_history.append((it, rnorm))

        self.ksp.setMonitor(monitor)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)

    def setup_system_from_jacobian(self):

        self.lsbuilder.set_AMPI(self.linear_system.A)
        self.lsbuilder.set_RHS(self.linear_system.RHS)

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        reason = self.ksp.getConvergedReason()

        return reason

    def get_iteration_number(self):

        return self.ksp.getIterationNumber()

    def check_residual_norm(self):

        self.linear_system.check_residual_norm()

    def dump_system(self, basename="", binary=False, comm=PETSc.COMM_WORLD):
        """
        Writes the linear system (Matrix, solution and RHS) in three different files

        :param basename: common part of the file names
        :param binary: if False, writes the files in ASCII format (good for text display)
                       if True, uses the binary format of PETSc (faster, good for loading the PETSc Object in an external script)
        :comm: MPI communicator
        """
        makeviewer = (
            PETSc.Viewer().createBinary if binary else PETSc.Viewer().createASCII
        )
        self.linear_system._dump(basename, makeviewer, comm)

    def get_solution(self):

        return self.linear_system.x

    def set_parameters(self, tol=None, maxit=None, restart=None):

        self.rtol = tol or self.rtol
        self.maxit = maxit or self.maxit
        self.restart = restart or self.maxit
        if tol or maxit:
            self.ksp.setTolerances(rtol=self.rtol, max_it=self.maxit)
        if restart:
            self.ksp.setGMRESRestart(self.restart)
