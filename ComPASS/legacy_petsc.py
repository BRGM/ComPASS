from ._kernel import get_kernel
from .simulation import state
import petsc4py
import sys

petsc4py.init()
from petsc4py import PETSc


class LegacyLinearSolver:
    """
    A structure used to call the fortran
    core functions for linear system solving
    """

    def __init__(
        self,
        tol=1e-6,
        maxit=150,
        restart=None,
        activate_cpramg=True,
        activate_direct_solver=False,
    ):

        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.activate_cpramg = activate_cpramg
        self.activate_direct_solver = activate_direct_solver
        self.tol = tol
        self.maxit = maxit
        self.restart = restart or maxit
        self.kernel = get_kernel()

        x = PETSc.Vec()
        x.createMPI((state.info.system.local_nb_cols, state.info.system.global_nb_cols))
        x.set(0)
        x.assemblyBegin()
        x.assemblyEnd()
        self.x = x

        self.kernel.SolvePetsc_Init(
            self.maxit, self.tol, self.activate_cpramg, self.activate_direct_solver
        )

    def set_values(self):

        self.kernel.SolvePetsc_SetUp()

    def solve(self):

        return self.kernel.SolvePetsc_ksp_solve(self.x)

    def get_iteration_number(self):

        return self.kernel.SolvePetsc_KspSolveIterationNumber()

    def check_solution(self):

        self.kernel.SolvePetsc_check_solution(self.x)

    def dump_system(self, basename=""):

        self.kernel.SolvePetsc_dump_system(basename)

    def get_solution(self):

        return self.x

    def set(self, tol=None, maxit=None, restart=None):

        self.tol = tol or self.tol
        self.maxit = maxit or self.maxit
        self.restart = restart or self.restart
        if (tol or maxit or restart) and (not self.activate_direct_solver):
            self.kernel.SolvePetsc_Ksp_configuration(self.tol, self.maxit, self.restart)
