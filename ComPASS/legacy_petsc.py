from ._kernel import get_kernel
import petsc4py
import sys

petsc4py.init()
from . import mpi
from petsc4py import PETSc


class LegacyLinearSystem:
    """
    A ghost structure used to mimic the PetscLinearSystem class
    """

    def __init__(self, simulation):

        x = PETSc.Vec()
        x.createMPI(
            (
                simulation.info.system.local_nb_cols,
                simulation.info.system.global_nb_cols,
            )
        )
        x.set(0)
        x.assemblyBegin()
        x.assemblyEnd()
        self.x = x
        self.kernel = get_kernel()

    def check_residual_norm(self):

        self.kernel.SolvePetsc_check_solution(self.x)

    def set_from_jacobian(self):

        self.kernel.SolvePetsc_SetUp()

    def dump_ascii(self, basename="", comm=PETSc.COMM_WORLD):

        """
        Writes the linear system (Matrix, solution and RHS) in three different files in ASCII format

        :param basename: common part of the file names
        :comm: MPI communicator
        """

        self.kernel.SolvePetsc_dump_system(basename)
        viewer = PETSc.Viewer().createASCII(basename + "x" + ".dat", "w", comm)
        self.x.view(viewer)
        viewer.destroy

    def dump_binary(self, basename="", comm=PETSc.COMM_WORLD):

        mpi.master_print(
            "Binary_dump is not available in the legacy linear solver\nPerforming an ASCII dump instead"
        )
        self.dump_ascii(basename, comm=PETSc.COMM_WORLD)


class LegacyLinearSolver:
    """
    A structure used to call the fortran
    core functions for linear system solving
    """

    def __init__(
        self,
        simulation,
        tol=1e-6,
        maxit=150,
        restart=None,
        activate_cpramg=True,
        activate_direct_solver=False,
        comm=None,
    ):

        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.activate_cpramg = activate_cpramg
        self.activate_direct_solver = activate_direct_solver
        self.tol = tol
        self.maxit = maxit
        self.restart = restart or maxit
        self.linear_system = LegacyLinearSystem(simulation)
        self.kernel = get_kernel()
        self.kernel.SolvePetsc_Init(
            self.maxit, self.tol, self.activate_cpramg, self.activate_direct_solver
        )

    def solve(self):

        return self.kernel.SolvePetsc_ksp_solve(self.linear_system.x)

    def get_iteration_number(self):

        return self.kernel.SolvePetsc_KspSolveIterationNumber()

    def set_parameters(self, tol=None, maxit=None, restart=None):

        self.tol = tol or self.tol
        self.maxit = maxit or self.maxit
        self.restart = restart or self.restart
        if (tol or maxit or restart) and (not self.activate_direct_solver):
            self.kernel.SolvePetsc_Ksp_configuration(self.tol, self.maxit, self.restart)
