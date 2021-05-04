#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
from .__init__ import *
from ..options import get, get_bool
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
from .preconditioners import CPRAMG, BlockJacobi
from .solver import IterativeSolverSettings


def linear_solver(
    simulation,
    legacy=True,
    direct=False,
    activate_cpramg=None,
    cpr_amg_type=None,
    tolerance=None,
    max_iterations=None,
    restart_size=None,
    from_options=False,
):
    """
    A function that manages linear solver instanciation from keyword parameters
    """

    # Command line options override the function's arguments
    if from_options == True:
        legacy_opt = get("--linear_solver_version", "legacy")
        legacy = True if legacy_opt == "legacy" else False
        direct = get_bool("--direct_linear_solver")
        activate_cpramg = not get_bool("--disable_cpramg")
        cpr_amg_type = get("--cpr_amg_type")

    if direct:
        if any((tolerance, max_iterations, restart_size, cpr_amg_type,)):
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
            if cpr_amg_type not in (None, "hypre"):
                mpi.master_print(
                    "Legacy solver only supports Hypre AMG in CPR-AMG preconditioner"
                )
            return LegacyIterativeSolver(
                LegacyLinearSystem(simulation), settings, activate_cpramg
            )
        else:
            linear_system = PetscLinearSystem(simulation)
            if activate_cpramg:
                cpr_amg_type = cpr_amg_type if cpr_amg_type is not None else "hypre"
                pc = CPRAMG(linear_system, amg_type=cpr_amg_type)
            else:
                pc = BlockJacobi(linear_system)
            return PetscIterativeSolver(linear_system, settings, pc=pc)
