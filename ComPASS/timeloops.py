#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
import logging

from time import process_time

import numpy as np

from sortedcontainers import SortedKeyList

from .utils.units import day, year
from .timestep_management import FixedTimeStep, TimeStepManager
from .utils.units import time_string
from . import timestep
from . import mpi
from .newton import Newton, default_Newton
from .dumps import Dumper
from .callbacks import get_callbacks_from_options
from .utils.units import bar, year
from . import messages
from ComPASS._kernel import get_kernel

Event = namedtuple("Event", ["time", "actions"])
TimeloopTick = namedtuple(
    "TimeloopTick", ["time", "iteration", "latest_timestep"], defaults=(None,) * 2
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
    allreduce = mpi.communicator().allreduce
    lp_min = simulation.all_states().p.min()
    p_min = allreduce(lp_min, mpi.MPI.MIN)
    lp_max = simulation.all_states().p.max()
    p_max = allreduce(lp_max, mpi.MPI.MAX)
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
shooter = None  # FIXME: global variable

# FIXME: temporary workaround for bug #280
# https://gitlab.inria.fr/compass/v4/ComPASS/-/issues/280
default_newton_instance = None


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
    reset_iteration_counter=False,
    dumper=None,
    iteration_callbacks=None,
    output_callbacks=None,
    newton_iteration_callbacks=None,
    specific_outputs=None,
    newton=None,
    time_step_manager=None,
    well_pressure_offset=1 * bar,
    events=None,
    output_before_start=True,
    output_after_loop=True,
    well_connections=None,
    no_output=False,
    respect_final_time=True,
    display_residual_contributions=False,
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
    :param reset_iteration_counter: Reset iteration counter if True (default to false)
    :param dumper: The object used to dump simulation (snaphots).
    :param iteration_callbacks: A sequence that holds callbacks that will be called after each iteration.
        The callback signature must be `f(tick)` where `tick` is compliant with the :py:class:`TimeloopTick`.
    :param output_callbacks: A sequence that holds callbacks that will be called before each simulation ouput
        (cf. ``ouput_period`` and ``output_every``).
        The callback signature must be `f(tick)` where `tick` is compliant with the :py:class:`TimeloopTick`.
    :param newton_iteration_callbacks: A sequence that holds callbacks that will be called after each newton iteration.
        The callback must have an argument which is supposed to be a NewtonLoopTick.
    :param specific_outputs: A sequence of additional output times.
    :param newton: A :class:`ComPASS.newton.Newton` object. If not provided a default one will be created
        by the :func:`~ComPASS.simulation.base.default_Newton` function.
    :param time_step_manager: A specfic time manager that will override the parameters ``fixed_timestep``
        or  ``initial_timestep`` (cf. :mod:`ComPASS.timestep_management` module).
    :param well_pressure_offset: if given the corresponfing pressure offset will be apply to wells
        (this may be useful to init wells operating on pressure) - if None, nothing is done. defaults to 1 bar.
    :param events: Is a sequence of events with actions to be process. The events must be compliant with
        the ComPASS.timeloops.Event namedtuple object.
    :param output_before_start: A boolean to specify if an output should be dumped before starting loop (default is True).
    :param output_after_loop: A boolean to specify if an output should be dumped when the loop is finished (default is True).
    :param well_connections: Well connections to be synchronized at the end of each time loop.
                             It must be a :py:class:`WellDataConnections` (defaults to simulation.well_connections).
    :param no_output: Flag that will prevent any output (defaults to False)
    :param respect_final_time: Flag to adapt the final timestep to respect final time (defaults to False).
    :param display_residual_contributions: Detail contributions to the residual norm during
                                           newton convergence (defaults to False).
    :return: The time at the end of the time loop.
    """
    # FIXME: horrible global variables... to be removed... using OOP?
    global n, shooter, default_newton_instance
    t0 = initial_time or 0
    total_time = final_time
    if total_time:
        total_time -= t0
    if reset_iteration_counter:
        n = 0
    assert not (final_time is None and nitermax is None)
    if newton is None:
        if default_newton_instance is None:
            default_newton_instance = default_Newton(simulation)
        newton = default_newton_instance
        assert (
            newton.simulation == simulation
        ), "Unconsistent simulation object for default Newton instance"
    if well_connections is None:
        well_connections = simulation.well_connections
    if output_period is None:
        if nb_output is not None:
            nb_output = max(2, nb_output)
            if total_time is not None:
                output_period = total_time / (nb_output - 1)
            else:
                messages.warning(
                    "nb_output has no impact because final_time is not set in standard_loop."
                )
    else:
        if nb_output is not None:
            messages.warning("output_period is overriding nb_output in standard_loop.")
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
    tick0 = TimeloopTick(time=t0, iteration=n)
    if iteration_callbacks is not None:
        iteration_callbacks = tuple(cb for cb in iteration_callbacks)
    else:
        iteration_callbacks = ()
    iteration_callbacks += get_callbacks_from_options(newton, tick0, no_output)
    if output_callbacks is None:
        output_callbacks = tuple()
    if specific_outputs is None:
        specific_outputs = []
    # FIXME: cf. #411 we should use list
    if newton_iteration_callbacks is not None:
        if newton.iteration_callbacks is None:
            newton.iteration_callbacks = tuple(newton_iteration_callbacks)
        else:
            newton.iteration_callbacks = newton.iteration_callbacks + tuple(
                newton_iteration_callbacks
            )
    events = _make_event_list(events)
    while len(events) > 0 and events[0].time < t0:
        mpi.master_print(f"WARNING: Event at time {events[0].time} is forgotten.")
        events.pop(0)
    if shooter is None:
        if dumper is None:
            dumper = Dumper(simulation)
            # as dumper will create directories slave must wait for directory creation by master
            mpi.synchronize()
        shooter = Snapshooter(dumper)

    def output_actions(tick):
        if no_output:
            return
        t, n, _ = tick
        for callback in output_callbacks:
            # FIXME: why not callback(tick)?
            callback(n, t)
        shooter.shoot(t)

    events.update(
        [
            Event(
                tout,
                [
                    output_actions,
                ],
            )
            for tout in specific_outputs
        ]
    )

    def add_output_event(tout):
        events.add(
            Event(
                tout,
                [
                    output_actions,
                ],
            )
        )

    # FIXME: use tick as argument not (t, n)
    @mpi.on_master_proc
    def print_iteration_info(t, n):
        print()
        print("** Time Step (iteration):", n, "*" * 50)
        if final_time:
            final_time_info = (
                "-> "
                + time_string(final_time)
                + f" ({100*(t - t0)/total_time:.2f}% done)"
            )
        else:
            final_time_info = "-> NO final time"
        print("Current time:", time_string(t), final_time_info)
        # print('Timestep: %.3g s = %.3f d = %.3f y' % (
        # ts_manager.current_step, ts_manager.current_step / day, ts_manager.current_step / year))

    pcsp = np.copy(simulation.cell_states().p)
    pcsT = np.copy(simulation.cell_states().T)
    pcsS = np.copy(simulation.cell_states().S)
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

        push_reccuring_output_event(tick0)
    if output_before_start:
        output_actions(tick0)
    process_events(tick0)

    if well_pressure_offset is not None:
        check_well_pressure(simulation, well_pressure_offset)
    t = t0
    tick = tick0
    cpu_start_time = process_time()
    while (final_time is None or t < final_time) and (nitermax is None or n < nitermax):
        dt_to_next_event = None
        if len(events) > 0:
            dt_to_next_event = events[0].time - t
        if final_time is not None and respect_final_time:
            if dt_to_next_event is None:
                dt_to_next_event = final_time - t
            else:
                dt_to_next_event = min(dt_to_next_event, final_time - t)
        assert (
            dt_to_next_event is None or dt_to_next_event > 0
        ), f"dt to next event: {dt_to_next_event}"
        n += 1
        newton.iterations = n
        print_iteration_info(t, n)
        # --
        mpi.synchronize()
        dt = timestep.make_one_timestep(
            simulation,
            newton,
            TimeloopTick(t, n),
            ts_manager.steps(upper_bound=dt_to_next_event),
            display_residual_contributions=display_residual_contributions,
        )
        well_connections.synchronize()
        mpi.synchronize()

        assert (
            dt == ts_manager.current_step
        ), f"Timesteps differ: {dt} vs {ts_manager.current_step}"
        t += dt
        tick = TimeloopTick(time=t, iteration=n, latest_timestep=dt)

        allreduce = mpi.communicator().allreduce
        mpi.master_print(
            "max p variation",
            allreduce(np.fabs(simulation.cell_states().p - pcsp).max(), mpi.MPI.MAX),
        )
        mpi.master_print(
            "max T variation",
            allreduce(np.fabs(simulation.cell_states().T - pcsT).max(), mpi.MPI.MAX),
        )
        mpi.master_print(
            "max S variation",
            allreduce(np.fabs(simulation.cell_states().S - pcsS).max(), mpi.MPI.MAX),
        )
        if output_every is not None and n % output_every == 0:
            add_output_event(t)
        process_events(tick)
        for callback in iteration_callbacks:
            try:
                callback(tick)
            except TypeError:
                mpi.master_print(
                    "WARNING: use the TimeloopTick API for iteration callbacks"
                )
                callback(n, t)
        pcsp[:] = simulation.cell_states().p
        pcsT[:] = simulation.cell_states().T
        #################################################################################################################
        # Print mswell info every time iteration
        kernel = get_kernel()
        kernel.IncCVMSWells_print_info_to_file()  # directive _DEBUG_INCCVM_MSWELLS_ needs to be set to 1 in fortran file
        kernel.JacobianMSWells_print_IP_info_to_file(
            -1
        )  # directive _DEBUG_JAC_IP_MSWELLS_ needs to be set to 1 in fortran file
        #################################################################################################################
    if output_after_loop:
        output_actions(tick)
    # Check if some events were left unprocessed
    if mpi.is_on_master_proc:
        mpi.master_print(
            f"Elapsed computing time (including output): {process_time() - cpu_start_time}"
        )
        for event in events:
            mpi.master_print(
                f"WARNING: Event at time {event.time} = {event.time / year} y was not reached."
            )
    mpi.synchronize()  # is this useful ?
    return t
