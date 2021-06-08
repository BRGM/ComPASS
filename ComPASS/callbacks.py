from .options import compass_config
from .runtime import to_output_directory
import os
from . import mpi


class Flag:
    def __init__(self):
        self.is_on = False

    def __call__(self, *args, **kwargs):
        pass


class AlwaysOnFlag(Flag):
    def __init__(self):
        super().__init__()
        self.is_on = True


class TimestepFlag(Flag):
    """
    A flag that will turn on at the first timestep where tick.time > self.t
    and turn off for the rest of the timeloop
    """

    def __init__(self, t):
        super().__init__()
        self.has_been_triggered = False
        self.t = t

    def __call__(self, tick):

        if not self.has_been_triggered:
            if tick.time > self.t:
                if not self.is_on:
                    self.is_on = True
                    self.t = tick.time
                elif self.is_on:
                    self.is_on = False
                    self.has_been_triggered = True


class NewtonDumper:
    def __init__(self, dump_functions, basedir, flag=None):

        self.dump_functions = dump_functions
        self.basedir = to_output_directory(basedir)
        self.flag = flag or AlwaysOnFlag()

    def __call__(self, newton_tick):

        time_tick = newton_tick.timeloop_tick
        self.flag(time_tick)

        if self.flag.is_on:
            tit = time_tick.iteration
            time = time_tick.time
            dt = newton_tick.current_dt
            nit = newton_tick.iteration
            dump_dirname = (
                f"{self.basedir}/it={tit}_t={time:.4e}/dt={dt:.4e}/newton_it={nit}"
            )

            try:
                os.makedirs(dump_dirname)
            except FileExistsError:
                pass

            for function in self.dump_functions:
                function(basename=dump_dirname)


class InterruptTrigger:
    def __init__(self, message=None):
        self.message = message

    def __call__(self, *args, **kwargs):
        mpi.master_print(f"\n ComPASS - Abortion requested {self.message}\n")
        mpi.abort()


class TimestepInterruptTrigger(InterruptTrigger):
    def __init__(self, timeflag, message=None):
        self.flag = timeflag
        self.message = message or f" at time t > {self.flag.t:.4e} s\n"

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.is_on:
            super().__call__()


class NewtonLogCallback:
    def __init__(self, filename, newton):
        self.filename = to_output_directory(filename)
        self.newton = newton
        with open(self.filename, "w+") as f:
            f.write(
                "Iteration number, Time, Successful timestep, Newton iterations, Linear iterations\n\n"
            )

    def __call__(self, tick):
        def writer(tick):
            with open(self.filename, "a") as f:
                f.write(
                    f"{tick.iteration},   {tick.time:.13e}, {tick.latest_timestep:.13e}, {len(self.newton.lsolver_iterations)}, {self.newton.lsolver_iterations}\n"
                )

        mpi.on_master_proc(writer)(tick)


def get_callbacks_from_options(newton, tick0):

    timestep_callbacks = []
    newton_iteration_callbacks = []
    newton_failure_callbacks = []
    linear_failure_callbacks = []
    linear_system = newton.lsolver.linear_system

    if compass_config.get("callbacks.dump_system_on_linear_failure"):
        dump_trigger = NewtonDumper(
            (linear_system.dump_binary, newton.lsolver.write_history),
            "linear_systems/",
        )
        linear_failure_callbacks.extend([dump_trigger])
    if compass_config.get("callbacks.abort_on_linear_failure"):
        linear_failure_callbacks.append(InterruptTrigger("on linear failure"))
    if compass_config.get("callbacks.abort_on_newton_failure"):
        newton_failure_callbacks.append(InterruptTrigger("on Newton failure"))

    if compass_config.get("callbacks.linear_system_dump"):
        t_dump_raw = compass_config["callbacks.linear_system_dump"]
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            dump_trigger = NewtonDumper(
                (linear_system.dump_ascii, newton.lsolver.write_history),
                "linear_systems/",
                flag=TimestepFlag(float(t_dump)),
            )
            newton_iteration_callbacks.extend([dump_trigger])

    if compass_config.get("callbacks.linear_system_binary_dump"):
        t_dump_raw = compass_config["callbacks.linear_system_binary_dump"]
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            dump_trigger = NewtonDumper(
                (linear_system.dump_binary, newton.lsolver.write_history),
                "linear_systems/",
                flag=TimestepFlag(float(t_dump)),
            )
            newton_iteration_callbacks.extend([dump_trigger])

    if compass_config.get("callbacks.newton_log"):
        newton_log_filename = compass_config["callbacks.newton_log"]
        newton_log_callback = NewtonLogCallback(newton_log_filename, newton)
        timestep_callbacks.append(newton_log_callback)

    if compass_config.get("callbacks.abort") is not None:
        t_kill = compass_config["callbacks.abort"]
        t_kill = float(t_kill)
        kill_flag = TimestepFlag(t_kill)
        kill_trigger = TimestepInterruptTrigger(kill_flag)
        timestep_callbacks.append(kill_trigger)

    newton.iteration_callbacks += tuple(newton_iteration_callbacks)
    newton.failure_callbacks += tuple(newton_failure_callbacks)
    newton.lsolver.failure_callbacks += tuple(linear_failure_callbacks)
    return tuple(timestep_callbacks)
