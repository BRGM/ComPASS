#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .__init__ import *
from .solver import IterativeSolver, DirectSolver, IterativeSolverSettings
from .._kernel import get_kernel
from .. import mpi
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

        mpi.master_print(">> Linear system dump")
        self.kernel.SolvePetsc_dump_system(basename + "/")
        viewer = PETSc.Viewer().createASCII(f"{basename}/x.dat", "w", comm)
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
        self,
        linear_system,
        settings,
        activate_cpramg=True,
    ):
        """
        :param settings: An IterativeSolverSettings object containing the wanted parameters for iterative solving
        """
        self.kernel = get_kernel()
        self.activate_cpramg = activate_cpramg
        super().__init__(linear_system, settings)
        self.kernel.SolvePetsc_Init(
            settings.max_iterations,
            settings.tolerance,
            self.activate_cpramg,
            False,
        )
        # Fortran doesn't set restart at initialization
        self.restart_size = self.my_settings.restart_size

    def set_settings(self, settings):
        self.my_settings = settings
        self.kernel.SolvePetsc_Ksp_configuration(
            self.my_settings.tolerance,
            self.my_settings.max_iterations,
            self.my_settings.restart_size,
        )

    settings = property(
        fget=lambda self: self.my_settings,
        fset=lambda self, value: self.set_settings(value),
        doc="Iterative solver parameters",
    )
    method = property(
        fget=lambda self: self.my_settings.method,
        fset=lambda self, value: mpi.master_print(
            "Method parameter of the legacy solver cannot be reset after initialization"
        ),
        doc="Iterative method used for linear solving",
    )
    tolerance = property(
        fget=lambda self: self.my_settings.tolerance,
        fset=lambda self, value: self.set_settings(
            IterativeSolverSettings(
                method=self.method,
                tolerance=value,
                max_iterations=self.max_iterations,
                restart_size=self.restart_size,
            )
        ),
        doc="Relative decrease in the residual norm required for convergence",
    )
    max_iterations = property(
        fget=lambda self: self.my_settings.max_iterations,
        fset=lambda self, value: self.set_settings(
            IterativeSolverSettings(
                method=self.method,
                tolerance=self.tolerance,
                max_iterations=value,
                restart_size=self.restart_size,
            )
        ),
        doc="Maximum number of iterations accepted before convergence failure",
    )
    restart_size = property(
        fget=lambda self: self.my_settings.restart_size,
        fset=lambda self, value: self.set_settings(
            IterativeSolverSettings(
                method=self.method,
                tolerance=self.tolerance,
                max_iterations=self.max_iterations,
                restart_size=value,
            )
        ),
        doc="Number of iterations at which GMRES restarts",
    )

    def solve(self):

        self.ksp_reason = self.kernel.SolvePetsc_ksp_solve(self.linear_system.x)
        self.nit = self.kernel.SolvePetsc_KspSolveIterationNumber()

        if self.ksp_reason < 0:
            self.number_of_unsuccessful_iterations += self.nit
            raise IterativeSolverFailure(self.ksp_reason, self.nit)
        else:
            self.number_of_successful_iterations += self.nit

        return self.linear_system.x, self.nit

    def __str__(self):

        cpramg_description = (
            "activated" if self.activate_cpramg == True else "not activated"
        )
        return f"{super().__str__()}\n   Legacy Fortran90 implementation\n   Settings : {self.settings}\n   Preconditioner : CPR-AMG {cpramg_description}"


class LegacyDirectSolver(DirectSolver):
    """
    A structure that holds a direct Legacy KSP Object to solve the linear system
    """

    def __init__(
        self,
        linear_system,
        comm=PETSc.COMM_WORLD,
    ):
        super().__init__(linear_system)
        self.linear_system = linear_system
        self.kernel = get_kernel()
        self.activate_cpramg = False
        self.kernel.SolvePetsc_Init(0, 0.0, False, True)

    def solve(self):

        self.ksp_reason = self.kernel.SolvePetsc_ksp_solve(self.linear_system.x)
        self.nit = 1

        if self.ksp_reason < 0:
            raise DirectSolverFailure(
                f"Legacy KSP object returned error code: {self.ksp_reason}"
            )

        return self.linear_system.x, self.nit

    def __str__(self):
        return f"{super().__str__()}\n   Legacy Fortran90 implementation"
