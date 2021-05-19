#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
from .__init__ import *
from .. import options
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
    iterative_method=None,
    activate_cpramg=None,
    cpr_amg_type=None,
    tolerance=None,
    max_iterations=None,
    restart_size=None,
    from_options=False,
):
    """
    A function that manages linear solver instanciation from keyword parameters and command line options

    :param simulation: an instanciated simulation object
    :param legacy: Switch between the Fortran version of the solver (True) or the new petsc4py implementation (False), defaults to True
    :param direct: Switch between a direct LU solver (True) and GMRES (False), defaults to False
    :param iterative_method: Iterative method used to solve linear systems (see https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html)
    :param activate_cpramg: Turn on (True) or off (False) the use of the CPR-AMG preconditioner, defaults to True (even if set to None)
    :param cpr_amg_type: Switch between Hypre BoomerAMG ("hypre") and PETSc's built-in ("gamg") AMG procedures, defaults to "hypre"
    :param tolerance: Relative decrease in residual required for convergence (iterative solvers only), defaults to 1e-6
    :param max_iterations: Maximum number of iterations accepted before divergence (iterative solvers only), defaults to 150
    :param restart_size: Number of iterations at which GMRES restarts (iterative solvers only), defaults to 30
    :param from_options: Turn on (True) or off (False) the use of command line options, defaults to False
    If set to True, runtime command line options will override the function's arguments
    """

    if from_options == True:
        legacy_opt = options.database["linear_solver_version"]
        legacy = (legacy_opt == "legacy") if legacy_opt is not None else legacy
        direct = options.database["direct_linear_solver"] or direct
        activate_cpramg = (
            False if options.database["disable_cpramg"] else activate_cpramg
        )
        cpr_amg_type = options.database["cpr_amg_type"] or cpr_amg_type

    if direct:
        if any(
            (
                tolerance,
                max_iterations,
                restart_size,
                cpr_amg_type,
                activate_cpramg in (False, True),
            )
        ):
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
            iterative_method or "gmres",
            tolerance or 1e-6,
            max_iterations or 150,
            restart_size or 30,
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
