from petsc4py import PETSc


class PetscElements:
    def __init__(self, system):
        # CHECKME: is there a way to check PETSc has been initialized here
        assert system.kernel is not None
        self.system = system
        x = PETSc.Vec()
        x.createMPI((system.local_nb_cols, system.global_nb_cols))
        x.set(0)
        x.assemblyBegin()
        x.assemblyEnd()
        self.x = x

    def solve(self):
        return self.system.kernel.SolvePetsc_ksp_solve(self.x)

    def check_solution(self):
        self.system.kernel.SolvePetsc_check_solution(self.x)
