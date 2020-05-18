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
        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setTolerances(rtol=1e-6, max_it=150)

    def setUp(self, simulation):

        self.createAMPI(simulation)
        self.setAMPI(simulation)
        self.setVecs(simulation)
        self.ksp.setOperators(self.A, self.A)
        self.ksp.setFromOptions()

    def SolvePetsc_solve(self):

        return self.system.kernel.SolvePetsc_ksp_solve(self.x)

    def solve(self):

        self.ksp.solve(self.RHS, self.x)
        return self.ksp.getConvergedReason()

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

    def dump_LinearSystem(self, binary=False):

        if binary:
            viewer = PETSc.Viewer().createBinary("A_binary.dat", "w", PETSc.COMM_WORLD)
            self.A.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createBinary("b_binary.dat", "w", PETSc.COMM_WORLD)
            self.RHS.view(viewer)
            viewer.destroy()
        else:
            viewer = PETSc.Viewer().createASCII("A.dat", "w", PETSc.COMM_WORLD)
            self.A.view(viewer)
            viewer.destroy()

            viewer = PETSc.Viewer().createASCII("b.dat", "w", PETSc.COMM_WORLD)
            self.RHS.view(viewer)
            viewer.destroy()
