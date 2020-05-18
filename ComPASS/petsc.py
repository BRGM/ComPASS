from petsc4py import PETSc
import numpy as np
import cProfile
from . import mpi


class LinearSystem:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self, system):
        # CHECKME: is there a way to check PETSc has been initialized here
        assert system.kernel is not None
        self.system = system
        comm = PETSc.COMM_WORLD
        self.x = PETSc.Vec().createMPI((system.local_nb_cols, system.global_nb_cols))
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        self.ksp.setTolerances(rtol=1e-6, max_it=150)

        def monitor(ksp, it, rnorm):
            self.residual_history.append((it, rnorm))

        self.ksp.setMonitor(monitor)
        # self.viewer = PETSc.Viewer().createASCII("ksp.dat", "a", PETSc.COMM_WORLD)
        self.residual_history = []

    def setUp(self, simulation):

        """fortran"""
        # self.system.kernel.SolvePetsc_SetUp() # For now linear solving is still operated by fortran SolvePetsc

        """self"""
        self.createAMPI(simulation)
        self.setAMPI(simulation)
        self.setVecs(simulation)
        self.ksp.setOperators(self.A, self.A)
        self.ksp.setFromOptions()

    def solve(self):

        """ Fortran """
        # reason = self.system.kernel.SolvePetsc_ksp_solve(self.x)
        # self.system.kernel.SolvePetsc_dump_system("t1/f90")

        """ self """
        self.residual_history = []
        self.ksp.solve(self.RHS, self.x)
        reason = self.ksp.getConvergedReason()
        # self.ksp.view(self.viewer)
        # mpi.master_print("lin_residuals :", self.residual_history[:])
        # self.dump_LinearSystem(basename="", binary = True)
        # 1/0

        return reason

    def check_solution(self):

        self.system.kernel.SolvePetsc_check_solution(self.x)

    def createAMPI(self, simulation):

        sizes, d_nnz, o_nnz = simulation.get_AMPI_nnz_cpp()
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes
        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))

    def setAMPI(self, simulation):

        simulation.set_AMPI_cpp(self.A)

    def setVecs(self, simulation):

        assert self.A.isAssembled(), "A must be assembled before calling setVecs"

        self.RHS = self.A.createVecs(side="right")
        self.x = self.A.createVecs(side="left")

        simulation.set_RHS_cpp(self.RHS)

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
