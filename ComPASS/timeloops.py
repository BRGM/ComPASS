#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

from .utils.units import day, year
from .timestep_management import FixedTimeStep, TimeStepManager
from .simulation_context import SimulationContext
from .utils.units import time_string
from . import timestep
from . import mpi
from .dumps import Dumper
from ._kernel import get_kernel


def check_well_pressure(simulation):
    p_min, p_max = float('Inf'), -float('Inf')
    for state in [simulation.node_states(), simulation.fracture_states(), simulation.cell_states()]:
        if state.p.shape[0]>0:
            p_min = min(p_min, state.p.min())
            p_max = max(p_max, state.p.max())
    # whp = well head pressure
    simulation.production_whp()[:] = p_min - 1.
    simulation.injection_whp()[:] = p_max + 1.

class Snapshooter:

    def __init__(self, dumper):
        dumper.start_simulation()
        self.dumper = dumper
        self.nb_snaphots = 0
        self.latest_snapshot_time = None
        self.filename = self.dumper.to_output_directory('snapshots')
        if mpi.is_on_master_proc:
            with open(self.filename, 'w') as f:
                pass

    def new_tag(self):
        tag = '%05d' % self.nb_snaphots
        self.nb_snaphots+= 1
        return tag

    def shoot(self, t):
        tag = self.new_tag()
        assert self.latest_snapshot_time is None or t > self.latest_snapshot_time
        self.latest_snapshot_time = t
        if mpi.is_on_master_proc:
            with open(self.filename, 'a') as f:
                print(tag, '%.12g'%t, file=f)
        self.dumper.dump_states(tag)

n = 0 # iteration counter
shooter = None

def standard_loop(simulation,
                  initial_time=None, final_time=None,
                  initial_timestep=None, fixed_timestep=None,
                  output_period = None, output_every = None,
                  nb_output = None, nitermax = None, dumper=None,
                  iteration_callbacks = None, output_callbacks = None,
                  specific_outputs = None,
                  newton = None, context = None,
                  time_step_manager = None
                 ):
    """
    Performs a standard timeloop.

    :param simulation: The simulation object that the time loop is acting on.    
    :param initial_time: Starting time of the simulation. Defaults to 0 if not specified.
    :param final_time: Final time of the simulation (if any).
    :param initial_timestep: Initial timestep in seconds, will activate automatic timestep management.
            (cf. :mod:`ComPASS.timestep_management` module).
            Cannot be specified along with ``fixed_timestep``.
    :param fixed_timestep: Fixed timestep in seconds, will disable automatic timestep management
            (cf. :mod:`ComPASS.timestep_management` module).
            Cannot be specified along with ``initial_timestep``.
    :param output_period: Will dump simulation information every ``output_period`` seconds.
    :param output_every: Will dump simulation information every ``output_every`` iterations.
    :param nb_output: Will compute ``output_period`` from ``final_time`` so that there are
            ``nb_output`` dumps of simulation info. This parameter will not have effect if 
            ``output_period`` is defined.
    :param nitermax: Maximum number of iterations.    
    :param dumper: The object used to dump simulation (snaphots).    
    :param iteration_callbacks: A sequence that holds callbacks that will be called after each iteration.
        The callback signature must be `f(n, t)` where `n` is the iteration number and `t` is 
        the current time.    
    :param output_callbacks: A sequence that holds callbacks that will be called before each simulation ouput
        (cf. ``ouput_period`` and ``output_every``).
        The callback signature must be `f(n, t)` where `n` is the iteration number and `t` is 
        the current time.
    :param specific_outputs: .    
    :param newton: A :class:`ComPASS.newton.Newton` object. If not provided a default one will be created 
        by the :func:`~ComPASS.simulation.base.default_Newton` function.
    :param context: A :class:`ComPASS.simulation_context.SimulationContext` object that is used to
        control generic option. If not provided a default one will be created.
    :param time_step_manager: A specfic time manager that will override the parameters ``fixed_timestep``
        or  ``initial_timestep`` (cf. :mod:`ComPASS.timestep_management` module).
    :return: The time at the end of the time loop.
    """
    assert not (final_time is None and nitermax is None)
    if newton is None:
        newton = simulation.default_Newton()
    if context is None:
        context = SimulationContext()
    global n
    global shooter
    if output_period is None:
        if nb_output is None:
            output_period = final_time    
        else:
            nb_output = max(2, nb_output)
            if final_time:
                output_period = (max(tstart, final_time) - tstart) / (nb_output - 1)
    else:
        if nb_output is not None:
            print('WARNING: output_period is overriding nb_output in standard_loop.')
    assert not(output_period is None or output_period<=0)
    if time_step_manager:
        assert initial_timestep is None and fixed_timestep is None
        ts_manager = time_step_manager
    else:
        assert initial_timestep is None or fixed_timestep is None
        assert initial_timestep or fixed_timestep
        if fixed_timestep:
            ts_manager = FixedTimeStep(fixed_timestep)
        else:
            ts_manager = TimeStepManager(
                    initial_timestep, max(initial_timestep, output_period),
            )
    if iteration_callbacks is None:
        iteration_callbacks = tuple()
    if output_callbacks is None:
        output_callbacks = tuple()
    if specific_outputs is None:
        specific_outputs = []
    else:
        specific_outputs = list(specific_outputs)
        specific_outputs.sort()
    # this is necessary for well operating on pressures
    check_well_pressure(simulation)
    # InitPressureDrop
    kernel = get_kernel()
    kernel.IncCVWells_InitPressureDrop()
   
    #FIXME: t = ComPASS.get_current_time()
    t = initial_time if initial_time is not None else 0
    t_output = t
    #if initial_time:
    #    #FIXME: ComPASS.set_current_time(initial_time)
    #    t = initial_time
    if shooter is None:
        if dumper is None:
            dumper = Dumper(simulation)
        shooter = Snapshooter(dumper)
    else:
        t_output = t_output + output_period
    @mpi.on_master_proc
    def print_iteration_info():
        print()
        print('** Time Step (iteration):', n, '*'*50)
        if final_time:
            final_time_info = '-> ' + time_string(final_time)
        else:
            final_time_info = '-> NO final time'
        print('Current time:', time_string(t), final_time_info)
        # print('Timestep: %.3g s = %.3f d = %.3f y' % (
            # ts_manager.current_step, ts_manager.current_step / day, ts_manager.current_step / year))
    pcsp = np.copy(simulation.cell_states().p)
    pcsT = np.copy(simulation.cell_states().T)
    while (final_time is None or t < final_time) and (nitermax is None or n < nitermax):
        if ( t < t_output and specific_outputs and 
             t + ts_manager.current_step > specific_outputs[0] ):
            # CHECKME: broken if fixed time step
            ts_manager.current_step = specific_outputs[0] - t
            t_output = specific_outputs[0]
            del specific_outputs[0]
        if t >= t_output or not(output_every is None or n%output_every>0):
            for callback in output_callbacks:
                callback(n, t)
            shooter.shoot(t)
            # WARNING / CHECKME we may loose some outputs
            while (t_output < t):
                t_output = t_output + output_period
        n+= 1
        print_iteration_info()
        # --
        timestep.make_one_timestep(
                newton, ts_manager.steps(),
                simulation_context=context,
        )
        t += ts_manager.current_step
        mpi.master_print('max p variation', np.abs(simulation.cell_states().p-pcsp).max())
        mpi.master_print('max T variation', np.abs(simulation.cell_states().T-pcsT).max())
        # --
        #ComPASS.make_timestep(timestep)
        #t = ComPASS.get_current_time()
        #if fixed_timestep is None:
        #    timestep = ComPASS.get_timestep()
        #ComPASS.timestep_summary()
        # --
        for callback in iteration_callbacks:
            callback(n, t)
        pcsp = np.copy(simulation.cell_states().p)
        pcsT = np.copy(simulation.cell_states().T)
    # Output final time
    if shooter.latest_snapshot_time is None or shooter.latest_snapshot_time < t:
        for callback in output_callbacks:
            callback(n, t)
        shooter.shoot(t)
    if mpi.is_on_master_proc:
        if specific_outputs:
            mpi.master_print('WARNING')
            mpi.master_print('WARNING: Specific outputs were not reached:', specific_outputs)
            mpi.master_print('WARNING')
    mpi.synchronize()
    return t
