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
from petsc4py import PETSc
from .solver import IterativeSolver, DirectSolver, IterativeSolverSettings
from .._kernel import get_kernel
from .exceptions import IterativeSolverFailure, DirectSolverFailure


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


class LegacyIterativeSolver(IterativeSolver):
    """
    A structure used to call the fortran
    core functions for linear system solving
    """

    def __init__(
        self, linear_system, settings, activate_cpramg=True,
    ):
        """
        :param settings: An IterativeSolverSettings object containing the wanted parameters for iterative solving
        """
        super().__init__(linear_system, settings)
        self.kernel = get_kernel()
        self.activate_cpramg = activate_cpramg
        self.settings = (
            self.settings._asdict()
        )  # Is there a better way than storing it as a dictionary ?
        self.kernel.SolvePetsc_Init(
            settings.max_iterations, settings.tolerance, self.activate_cpramg, False,
        )

    def make_setter(name):
        def setter(self, value):
            self.settings[name] = value
            self.kernel.SolvePetsc_Ksp_configuration(
                self.settings["tolerance"],
                self.settings["max_iterations"],
                self.settings["restart_size"],
            )

        return setter

    tolerance = property(
        fget=lambda self: self.settings["tolerance"],
        fset=make_setter("tolerance"),
        doc="Relative decrease in the residual norm required for convergence",
    )
    max_iterations = property(
        fget=lambda self: self.settings["max_iterations"],
        fset=make_setter("max_iterations"),
        doc="Maximum number of iterations accepted before convergence failure",
    )
    restart_size = property(
        fget=lambda self: self.settings["restart_size"],
        fset=make_setter("restart_size"),
        doc="Number of iterations at which GMRES restarts",
    )

    def solve(self):

        ksp_reason = self.kernel.SolvePetsc_ksp_solve(self.linear_system.x)
        nit = self.kernel.SolvePetsc_KspSolveIterationNumber()

        if ksp_reason < 0:
            self.number_of_unsuccessful_iterations += nit
            raise IterativeSolverFailure(
                f"Legacy KSP object returned error code: {ksp_reason}", nit
            )
        else:
            self.number_of_successful_iterations += nit

        return self.linear_system.x, nit, ksp_reason


class LegacyDirectSolver(DirectSolver):
    def __init__(
        self, linear_system, comm=None,
    ):
        super().__init__(linear_system)
        self.linear_system = linear_system
        self.kernel = get_kernel()
        self.activate_cpramg = False
        self.kernel.SolvePetsc_Init(0, 0.0, False, True)

    def solve(self):

        ksp_reason = self.kernel.SolvePetsc_ksp_solve(self.linear_system.x)
        nit = None

        if ksp_reason < 0:
            raise DirectSolverFailure(
                f"Legacy KSP object returned error code: {ksp_reason}"
            )

        return self.linear_system.x, nit, ksp_reason
