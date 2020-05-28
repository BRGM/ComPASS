import petsc4py
import sys

petsc4py.init(sys.argv)
from petsc4py import PETSc
from . import mpi


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

        (sizes, d_nnz, o_nnz) = nonzeros
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))
        self.x.createMPI((n_rowl, n_rowg))
        self.RHS.createMPI((n_rowl, n_rowg))

    def dump(self, basename="", binary=False):

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

    def checkSolution(self):

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
    A stucture that holds the PETSc KSP Object to solve the linear system
    """

    def __init__(self, simulation, tol=1e-6, maxit=150, restart=None):
        """
        :param tol: relative tolerance (for iterative solvers).
        :param maxit: maximum number of iterations (for iterative solvers).
        :param restart: number of iterations before a restart (for gmres like iterative solvers).
        """
        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.last_residual_history = []
        self.activate_direct_solver = False
        self.rtol = tol
        self.maxit = maxit
        self.restart = restart

        self.plsystem = PetscLinearSystem()
        self.lsbuilder = simulation.LinearSystemBuilder()
        self.plsystem.create(self.lsbuilder.get_non_zeros())

        comm = PETSc.COMM_WORLD
        self.ksp = PETSc.KSP().create(comm=comm)
        self.set(self.rtol, self.maxit, self.restart)
        self.last_residual_history = []
        self.ksp.setOperators(self.plsystem.A, self.plsystem.A)
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

    def set_values(self):

        self.lsbuilder.set_AMPI(self.plsystem.A)
        self.lsbuilder.set_RHS(self.plsystem.RHS)

    def solve(self):

        self.ksp.solve(self.plsystem.RHS, self.plsystem.x)
        reason = self.ksp.getConvergedReason()

        return reason

    def get_iteration_number(self):

        return self.ksp.getIterationNumber()

    def check_solution(self):

        self.plsystem.checkSolution()

    def dump_system(self, basename=""):

        self.plsystem.dump(basename=basename)

    def get_solution(self):

        return self.plsystem.x

    def set(self, tol=None, maxit=None, restart=None):

        self.rtol = tol or self.rtol
        self.maxit = maxit or self.maxit
        self.restart = restart or self.maxit
        if tol or maxit:
            self.ksp.setTolerances(rtol=self.rtol, max_it=self.maxit)
        if restart:
            self.ksp.setGMRESRestart(self.restart)
