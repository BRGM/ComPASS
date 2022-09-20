import numpy as np
from petsc4py import PETSc


class Synchronizer:
    def __init__(self, system, only_mswells=False):
        # CHECKME: is there a way to check PETSc has been initialized here
        assert system.kernel is not None
        self.system = system
        if not only_mswells:
            self.retrieve_fortran = system.kernel.SyncPetsc_GetSolNodeFracWellMSWell
        else:
            self.retrieve_fortran = system.kernel.SyncPetscMSWells_GetSolMSWell
        nbrows, nbcols = system.local_system_size
        gnbrows, gnbcols = system.global_system_size
        M = PETSc.Mat()
        sizes = (nbrows, gnbrows), (nbcols, gnbcols)
        csr = np.arange(nbrows + 1, dtype=np.int32), system.colnum(), np.ones(nbrows)
        M.createAIJ(sizes, csr=csr)
        M.assemblyBegin()
        M.assemblyEnd()
        x = PETSc.Vec()
        x.createMPI((nbrows, gnbrows))
        x.set(0)
        x.assemblyBegin()
        x.assemblyEnd()
        self.M_s = M
        self.x_s = x

    def synchronize(self, x):
        self.M_s.mult(x, self.x_s)

    # CHECKME: might be elsewhere...
    def retrieve_solution(self, newton_increments):
        self.retrieve_fortran(self.x_s, newton_increments)
