import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc
from . import mpi
from ._kernel import get_kernel
from . import mpi


class LinearSystem:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self):

        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()

    def createAMPI(self, nonzeros):

        (sizes, d_nnz, o_nnz) = nonzeros
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))

    def createVecs(self, sizes):

        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        self.x.createMPI((n_rowl, n_rowg))
        self.RHS.createMPI((n_rowl, n_rowg))

    def setAMPI(self, cpp_lsystem):

        cpp_lsystem.set_AMPI_cpp(self.A)

    def setRHS(self, cpp_lsystem):

        cpp_lsystem.set_RHS_cpp(self.RHS)

    def dump_LinearSystem(self, basename="", binary=False):

        if binary:
            viewer = PETSc.Viewer().createBinary(
                basename + "A_binary.dat", "w", PETSc.COMM_WORLD
            )
            self.A.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createBinary(
                basename + "b_binary.dat", "w", PETSc.COMM_WORLD
            )
            self.RHS.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createBinary(
                basename + "x.dat", "w", PETSc.COMM_WORLD
            )
            self.x.view(viewer)
            viewer.destroy()

        else:
            viewer = PETSc.Viewer().createASCII(
                basename + "A.dat", "w", PETSc.COMM_WORLD
            )
            self.A.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createASCII(
                basename + "b.dat", "w", PETSc.COMM_WORLD
            )
            self.RHS.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createASCII(
                basename + "x.dat", "w", PETSc.COMM_WORLD
            )
            self.x.view(viewer)
            viewer.destroy()


class LinearSolver:

    """
    A stucture that holds the PETSc KSP Object to solve the linear system
    """

    def __init__(self, tol, maxit, restart=None):
        """
        :param tol: relative tolerance (for iterative solvers).
        :param maxit: maximum number of iterations (for iterative solvers).
        :param restart: number of iterations before a restart (for gmres like iterative solvers).
        """
        self.lsystem = LinearSystem()
        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.last_residual_history = []
        self.activate_direct_solver = False
        self.rtol = tol
        self.maxit = maxit
        self.restart = restart

        comm = PETSc.COMM_WORLD
        self.ksp = PETSc.KSP().create(comm=comm)

    def initialize(self, cpp_lsystem):

        self.reset(self.rtol, self.maxit, self.restart)

        def monitor(ksp, it, rnorm):
            self.last_residual_history.append((it, rnorm))

        self.ksp.setMonitor(monitor)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        nonzeros = cpp_lsystem.get_non_zeros()
        self.lsystem.createAMPI(nonzeros)
        self.lsystem.createVecs(nonzeros[0])

    def setUp(self, cpp_lsystem):

        self.lsystem.setAMPI(cpp_lsystem)
        self.lsystem.setRHS(cpp_lsystem)
        self.ksp.setOperators(self.lsystem.A, self.lsystem.A)
        if self.activate_direct_solver:
            self.ksp.setType("preonly")
            self.ksp.getPC().setType("lu")
        else:
            self.ksp.getPC().setFactorLevels(1)
        self.ksp.setFromOptions()

    def solve(self):

        self.last_residual_history = []
        # print('Direct Solver', self.activate_direct_solver)
        # print('KSP Type :', self.ksp.type)
        # print('PS Type : ', self.ksp.getPC().type)
        self.ksp.solve(self.lsystem.RHS, self.lsystem.x)
        reason = self.ksp.getConvergedReason()

        return reason

    def checkSolution(self):

        y = self.lsystem.RHS.duplicate()
        self.lsystem.RHS.copy(y)  # y = b
        y.scale(-1.0)  # y = -b
        self.lsystem.A.multAdd(self.lsystem.x, y, y)  # y = Ax-b
        mpi.master_print("Linear solution check ||Ax-b||")
        norm = y.norm(PETSc.NormType.NORM_1)
        mpi.master_print("  1-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_2)
        mpi.master_print("  2-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_INFINITY)
        mpi.master_print("  Infinity norm", norm)

    def reset(self, tol, maxit, restart=None):

        self.rtol = tol
        self.maxit = maxit
        self.restart = restart or maxit
        self.ksp.setTolerances(rtol=self.rtol, max_it=self.maxit)
        self.ksp.setGMRESRestart(self.restart)
        self.ksp.setFromOptions()

    def getSolution(self):

        return self.lsystem.x
