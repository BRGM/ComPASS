#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
from .. import mpi
from .legacy_linear_solver import (
    LegacyDirectSolver,
    LegacyIterativeSolver,
    LegacyLinearSystem,
)
from .petsc_linear_solver import (
    PetscIterativeSolver,
    PetscDirectSolver,
    PetscLinearSystem,
)
from .solver import IterativeSolverSettings


def linear_solver(
    simulation,
    legacy=True,
    direct=False,
    activate_cpramg=None,
    tolerance=None,
    max_iterations=None,
    restart_size=None,
):
    """
    A function that manages linear solver instanciation from keyword parameters
    """
    if direct:
        if any((activate_cpramg, tolerance, max_iterations, restart_size)):
            mpi.master_print(
                "Invalid parameter(s) passed to linear_solver()\nIterative solver parameter(s) will be ignored in direct solver instanciation"
            )
        if legacy:
            return LegacyDirectSolver(LegacyLinearSystem(simulation))
        else:
            return PetscDirectSolver(PetscLinearSystem(simulation))
    else:
        # Default settings if not provided
        settings = IterativeSolverSettings(
            tolerance or 1e-6, max_iterations or 150, restart_size or 30,
        )
        # activate_cpramg defaults to True if not provided
        activate_cpramg = activate_cpramg if activate_cpramg is not None else True
        if legacy:
            return LegacyIterativeSolver(
                LegacyLinearSystem(simulation), settings, activate_cpramg
            )
        else:
            return PetscIterativeSolver(
                PetscLinearSystem(simulation), settings, activate_cpramg
            )
