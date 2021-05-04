#
# This file is part of kernel.
#
# kernel.is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from . import mpi
from .newton import IterationExhaustion, NewtonFailure, LinearSolverFailure
from .utils.units import time_string
from ._kernel import get_kernel

# FIXME: computation time spent is to measured at the caller site
# comptime_start = MPI_WTIME()
# comptime_timestep = MPI_WTIME() - comptime_start # this is to be measureed externally (on the caller side)
# comptime_total = comptime_total + comptime_timestep # this is to be measured externally (on the caller side)


class AllAttemptsFailed(Exception):
    def __init__(self, attempts):
        self.attempts = attempts


def try_timestep(
    deltat, newton, simulation_context, display_residual_contributions=False
):
    # CHECKME: do we need to retrieve kernel here???
    kernel = get_kernel()
    kernel.IncCV_SaveIncPreviousTimeStep()
    try:
        mpi.master_print("trying newton with timestep:", time_string(deltat))
        iterations = newton.loop(deltat, display_residual_contributions)
        mpi.master_print(iterations)
    except LinearSolverFailure as e:
        mpi.master_print(e)
        if simulation_context and simulation_context.dump_system_on_ksp_failure:
            newton.lsolver.linear_system.dump_ascii(basename="ksp_failure_")
        if simulation_context and simulation_context.abort_on_ksp_failure:
            mpi.master_print(
                "\nComPASS - Abortion requested on linear solver failure\n"
            )
            mpi.abort()
        return False
    except IterationExhaustion as e:
        if simulation_context and simulation_context.abort_on_newton_failure:
            mpi.master_print("\nComPASS - Abortion requested on newton failure\n")
            mpi.abort()
        return False
    return True


def _flash_wells(simulation):
    # CHECKME: do we need to retrieve kernel here???
    kernel = get_kernel()
    # FIXME: we should be able to have no well model at all
    #        (when there are no wells)
    if simulation.well_model == "single_phase":
        kernel.DefFlashWells_TimeFlash_producers_single_phase()
    else:
        assert simulation.well_model == "two_phases"
        kernel.DefFlashWells_TimeFlash_producers_two_phases()
    kernel.DefFlashWells_TimeFlash_injectors()


def make_one_timestep(
    simulation,
    newton,
    timesteps,
    simulation_context=None,
    display_residual_contributions=False,
):
    # CHECKME: do we need to retrieve kernel here???
    kernel = get_kernel()
    # pressure drop are set explicitely so we don't need to recompute
    # in case of newton failure
    if simulation.unknown_producers_density:
        kernel.IncCVWells_estimate_producers_density(False)
        simulation.unknown_producers_density = False
    kernel.IncCVWells_UpdatePressureDrop()
    attempts = []
    for deltat in timesteps:
        if try_timestep(
            deltat, newton, simulation_context, display_residual_contributions
        ):
            mpi.master_print("Success with timestep: ", time_string(deltat))
            break
        attempts.append(deltat)
        mpi.master_print("Failure with timestep: ", time_string(deltat))
        newton.failures += 1
        newton.reset_loop()
    else:
        raise AllAttemptsFailed(attempts)
    _flash_wells(simulation)
    return deltat
