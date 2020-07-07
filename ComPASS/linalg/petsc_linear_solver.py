#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import petsc4py
import sys

petsc4py.init(sys.argv)
from .. import mpi
from .solver import *
from petsc4py import PETSc
from .exceptions import IterativeSolverFailure, DirectSolverFailure


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


class PetscIterativeSolver(IterativeSolver):

    """
    A structure that holds an iterative PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, settings, comm=PETSc.COMM_WORLD,
    ):

        super().__init__(linear_system, settings)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        self.ksp.getPC().setFactorLevels(1)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        self.tolerance, self.max_iterations, self.restart_size = self.settings[:]
        self.ksp.setFromOptions()

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
    restart_size = property(
        # restart_size is not available in petsc4py,
        # but can be set using the command line petsc option -ksp_gmres_restart restart
        fget=None,
        fset=lambda self, value: self.ksp.setGMRESRestart(value),
        doc="Number of iterations at which GMRES restarts",
    )

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        ksp_reason = self.ksp.getConvergedReason()
        nit = self.ksp.getIterationNumber()

        if ksp_reason < 0:
            self.number_of_unsuccessful_iterations += nit
            raise IterativeSolverFailure(
                "Petsc KSP object returned error code: %s" % ksp_reason, nit
            )
        else:
            self.number_of_successful_iterations += nit

        return self.linear_system.x, nit, ksp_reason


class PetscDirectSolver(DirectSolver):
    """
    A structure that holds a direct PETSc KSP Object to solve the linear system
    """

    def __init__(
        self, linear_system, comm=PETSc.COMM_WORLD,
    ):

        super().__init__(linear_system)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        self.ksp.setType("preonly")
        self.ksp.getPC().setType("lu")
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        ksp_reason = self.ksp.getConvergedReason()

        if ksp_reason < 0:
            raise DirectSolverFailure(
                "Petsc KSP object returned error code: %s" % ksp_reason
            )

        return self.linear_system.x, None, ksp_reason
