#
# This file is part of kernel.
#
# kernel.is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.newton import KspFailure, IterationExhaustion, NewtonFailure
from ComPASS.utils.units import day, year

# FIXME: computation time spent is to measured at the caller site
# comptime_start = MPI_WTIME()
#comptime_timestep = MPI_WTIME() - comptime_start # this is to be measureed externally (on the caller side)
#comptime_total = comptime_total + comptime_timestep # this is to be measured externally (on the caller side)

class AllAttemptsFailed(Exception):
    def __init__(self, attempts):
        self.attempts = attempts


def try_timestep(
    deltat, newton, simulation_context=None,
):
    # CHECKME: do we need to retrieve kernel here???
    kernel = ComPASS.kernel
    assert ComPASS.kernel
    kernel.IncCV_SaveIncPreviousTimeStep()
    kernel.IncCVWells_PressureDrop()
    try:
        mpi.master_print('trying newton with timestep:', deltat)
        iterations = newton.loop(deltat)
        mpi.master_print(iterations)
    except KspFailure as e:        
        mpi.master_print(
            'KSP failure - with reason', e.reason, 'after',
            kernel.SolvePetsc_KspSolveIterationNumber(), 'iterations'
        )
        if simulation_context and simulation_context.dump_system_on_ksp_failure:
            kernel.SolvePetsc_dump_system('ksp_failure_%03d' % newton.lsolver.failures)
        if simulation_context and simulation_context.abort_on_ksp_failure:
            mpi.master_print('\nComPASS - Abortion requested on linear solver failure\n')
            mpi.abort()
        return False 
    except IterationExhaustion as e:
        if simulation_context and simulation_context.abort_on_newton_failure:
            mpi.master_print('\nComPASS - Abortion requested on newton failure\n')
            mpi.abort()
        return False
    return True

def make_one_timestep(
    newton, timesteps, simulation_context=None,
):
    attempts = []
    for deltat in timesteps:
        if try_timestep(deltat, newton, simulation_context):
            break
        attempts.append(deltat)
        newton.failures+= 1   
        newton.reset_loop()
    else:
        raise AllAttemptsFailed(attempts)
    # CHECKME: do we need to retrieve kernel here???
    kernel = ComPASS.kernel
    assert ComPASS.kernel
    kernel.DefFlashWells_TimeFlash()
    return deltat
