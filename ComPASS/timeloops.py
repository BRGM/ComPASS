#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
import logging

import numpy as np

from sortedcontainers import SortedKeyList

from .utils.units import day, year
from .timestep_management import FixedTimeStep, TimeStepManager
from .simulation_context import SimulationContext
from .utils.units import time_string
from . import timestep
from . import mpi
from .dumps import Dumper
from .options import get_callbacks_from_options
from ._kernel import get_kernel
from .utils.units import bar, year


Event = namedtuple("Event", ["time", "actions"])
LoopTick = namedtuple(
    "LoopTick", ["time", "iteration", "last_timestep"], defaults=(None,) * 2
)


def _make_event_list(events=None):
    if events is None:
        events = []
    return SortedKeyList(events, key=lambda event: event.time)


def check_well_pressure(simulation, pressure_offset=1 * bar):
    """
    Set all well head pressures:
        - to minimum reservoir pressure - pressure_offset for producers
        - to maximum reservoir pressure + pressure_offset for injectors

    :param simulation: the current simulation
    :param pressure_offset: the pressure offset to be applied (defaults to 1 bar)
    """
    p_min = simulation.all_states().p.min()
    p_max = simulation.all_states().p.max()
    # whp = well head pressure
    whp = simulation.production_whp()
    for wk, data in enumerate(simulation.producers_data()):
        if data.operating_code == "f" and data.imposed_flowrate == 0:
            logging.warning(
                f"Producer well {data.id} has no flow rate but is not marked as closed."
            )
            whp[wk] = p_max + pressure_offset
        else:
            whp[wk] = p_min - pressure_offset
    whp = simulation.injection_whp()
    for wk, data in enumerate(simulation.injectors_data()):
        if data.operating_code == "f" and data.imposed_flowrate == 0:
            logging.warning(
                f"Injector well {data.id} has no flow rate but is not marked as closed."
            )
            whp[wk] = p_min - pressure_offset
        else:
            whp[wk] = p_max + pressure_offset


class Snapshooter:
    def __init__(self, dumper):
        dumper.start_simulation()
        self.dumper = dumper
        self.nb_snaphots = 0
        self.latest_snapshot_time = None
        self.filename = self.dumper.to_output_directory("snapshots")
        if mpi.is_on_master_proc:
            with open(self.filename, "w") as f:
                pass

    def new_tag(self):
        tag = "%05d" % self.nb_snaphots
        self.nb_snaphots += 1
        return tag

    def shoot(self, t):
        tag = self.new_tag()
        assert self.latest_snapshot_time is None or t >= self.latest_snapshot_time
        if t == self.latest_snapshot_time:
            return
        self.latest_snapshot_time = t
        if mpi.is_on_master_proc:
            with open(self.filename, "a") as f:
                print(tag, "%.12g" % t, file=f)
                print(f"[Snapshooter] Saving output at time: {t} s = {t/year} y")

        self.dumper.dump_states(tag)


n = 0  # iteration counter FIXME: global variable
shooter = None


def standard_loop(
    simulation,
    initial_time=None,
    final_time=None,
    initial_timestep=None,
    fixed_timestep=None,
    output_period=None,
    output_every=None,
    nb_output=None,
    nitermax=None,
    dumper=None,
    iteration_callbacks=None,
    output_callbacks=None,
    specific_outputs=None,
    newton=None,
    context=None,
    time_step_manager=None,
    well_pressure_offset=1 * bar,
    events=None,
    output_before_start=True,
    output_after_loop=True,
    well_connections=None,
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
        The callback signature must be `f(tick)` where `tick` is compliant with the :py:class:`LoopTick`.
    :param output_callbacks: A sequence that holds callbacks that will be called before each simulation ouput
        (cf. ``ouput_period`` and ``output_every``).
        The callback signature must be `f(tick)` where `tick` is compliant with the :py:class:`LoopTick`.
    :param specific_outputs: .
    :param newton: A :class:`ComPASS.newton.Newton` object. If not provided a default one will be created
        by the :func:`~ComPASS.simulation.base.default_Newton` function.
    :param context: A :class:`ComPASS.simulation_context.SimulationContext` object that is used to
        control generic option. If not provided a default one will be created.
    :param time_step_manager: A specfic time manager that will override the parameters ``fixed_timestep``
        or  ``initial_timestep`` (cf. :mod:`ComPASS.timestep_management` module).
    :param well_pressure_offset: if given the corresponfing pressure offset will be apply to wells
        (this may be usefull to init wells operating on pressure) - if None, nothing is done. defaults to 1 bar.
    :param events: Is a sequence of events with actions to be process. The events must be compliant with
        the ComPASS.timeloops.Event namedtuple object.
    :param output_before_start: A boolean to specify if an output should be dumped before starting loop (default is True).
    :param output_after_loop: A boolean to specify if an output should be dumped when the loop is finished (default is True).
    :param well_connections: Well connections to be synchronized at the end of each time loop.
                             It must be a :py:class:`WellDataConnections` (defaults to simulation.well_connections).
    :return: The time at the end of the time loop.
    """
    assert not (final_time is None and nitermax is None)
    if newton is None:
        newton = simulation.default_Newton()
    if context is None:
        context = SimulationContext()
    if well_connections is None:
        well_connections = simulation.well_connections
    global n
    global shooter
    if output_period is None:
        if nb_output is not None:
            nb_output = max(2, nb_output)
            if final_time:
                output_period = (max(tstart, final_time) - tstart) / (nb_output - 1)
    else:
        if nb_output is not None:
            print("WARNING: output_period is overriding nb_output in standard_loop.")
    assert output_period is None or output_period > 0
    if time_step_manager:
        assert initial_timestep is None and fixed_timestep is None
        ts_manager = time_step_manager
    else:
        assert initial_timestep is None or fixed_timestep is None
        assert initial_timestep or fixed_timestep
        if fixed_timestep:
            ts_manager = FixedTimeStep(fixed_timestep)
        else:
            ts_manager = TimeStepManager(initial_timestep)
    if iteration_callbacks is None:
        iteration_callbacks = get_callbacks_from_options()
    else:
        iteration_callbacks = (
            tuple(callback for callback in iteration_callbacks)
            + get_callbacks_from_options()
        )
    if output_callbacks is None:
        output_callbacks = tuple()
    if specific_outputs is None:
        specific_outputs = []
    events = _make_event_list(events)
    if well_pressure_offset is not None:
        check_well_pressure(simulation, well_pressure_offset)
    # InitPressureDrop
    kernel = get_kernel()
    kernel.IncCVWells_InitPressureDrop()
    t = initial_time if initial_time is not None else 0
    while len(events) > 0 and events[0].time < t:
        mpi.master_print(f"WARNING: Event at time {events[0].time} is forgotten.")
        events.pop(0)
    if shooter is None:
        if dumper is None:
            dumper = Dumper(simulation)
        shooter = Snapshooter(dumper)

    def output_actions(tick):
        t, n, _ = tick
        for callback in output_callbacks:
            callback(n, t)
        shooter.shoot(t)

    events.update([Event(tout, [output_actions,]) for tout in specific_outputs])

    def add_output_event(tout):
        events.add(Event(tout, [output_actions,]))

    @mpi.on_master_proc
    def print_iteration_info():
        print()
        print("** Time Step (iteration):", n, "*" * 50)
        if final_time:
            final_time_info = "-> " + time_string(final_time)
        else:
            final_time_info = "-> NO final time"
        print("Current time:", time_string(t), final_time_info)
        # print('Timestep: %.3g s = %.3f d = %.3f y' % (
        # ts_manager.current_step, ts_manager.current_step / day, ts_manager.current_step / year))

    pcsp = np.copy(simulation.cell_states().p)
    pcsT = np.copy(simulation.cell_states().T)
    dt = None

    def process_events(tick):
        while len(events) > 0 and events[0].time <= tick.time:
            for action in events[0].actions:
                action(tick)
            events.pop(0)

    if output_period is not None:

        def push_reccuring_output_event(tick):
            t = tick.time
            if final_time is None or t + output_period <= final_time:
                events.add(
                    Event(
                        t + output_period, [output_actions, push_reccuring_output_event]
                    )
                )

        push_reccuring_output_event(LoopTick(time=t + output_period))
    if output_before_start:
        output_actions(LoopTick(time=t, iteration=n))
    process_events(LoopTick(time=t, iteration=n))
    while (final_time is None or t < final_time) and (nitermax is None or n < nitermax):
        dt_to_next_event = None
        if len(events) > 0:
            dt_to_next_event = events[0].time - t
        if final_time is not None:
            if dt_to_next_event is None:
                dt_to_next_event = final_time - t
            else:
                dt_to_next_event = min(dt_to_next_event, final_time - t)
        assert (
            dt_to_next_event is None or dt_to_next_event > 0
        ), f"dt to next event: {dt_to_next_event}"
        n += 1
        print_iteration_info()
        # --
        dt = timestep.make_one_timestep(
            newton,
            ts_manager.steps(upper_bound=dt_to_next_event),
            simulation_context=context,
        )
        well_connections.synchronize()
        assert (
            dt == ts_manager.current_step
        ), f"Timesteps differ: {dt} vs {ts_manager.current_step}"
        t += dt
        tick = LoopTick(time=t, iteration=n, last_timestep=dt)
        mpi.master_print(
            "max p variation", np.fabs(simulation.cell_states().p - pcsp).max()
        )
        mpi.master_print(
            "max T variation", np.fabs(simulation.cell_states().T - pcsT).max()
        )
        if output_every is not None and n % output_every == 0:
            add_output_event(t)
        process_events(tick)
        for callback in iteration_callbacks:
            try:
                callback(tick)
            except TypeError:
                mpi.master_print(
                    "WARNING: use the LoopTick API for iteration callbacks"
                )
                callback(n, t)
        pcsp[:] = simulation.cell_states().p
        pcsT[:] = simulation.cell_states().T
    if output_after_loop:
        output_actions(tick)
    # Check if some events were left unprocessed
    if mpi.is_on_master_proc:
        for event in events:
            mpi.master_print(
                f"WARNING: Event at time {events[0].time} = {events[0].time / year} y was not reached."
            )
    mpi.synchronize()  # is this usefull ?
    return t
