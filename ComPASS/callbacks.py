from . import options
from .runtime import to_output_directory
import os
from . import mpi


class Flag:
    def __init__(self):
        self.is_on = False


class AlwaysOnFlag(Flag):
    def __init__(self):
        super().__init__()
        self.is_on = True


class TimestepFlag(Flag):
    def __init__(self, t):
        super().__init__()
        self.has_been_triggered = False
        self.t = t

    def __call__(self, tick):

        if not self.has_been_triggered:
            self.tick = tick
            if tick.time >= self.t:
                if not self.is_on:
                    self.is_on = True
                elif self.is_on:
                    self.is_on = False
                    self.has_been_triggered = True


class DumpLinearSystemTrigger:
    def __init__(self, dump_function, flag):
        self.flag = flag
        self.dump_function = dump_function
        self.dirname = to_output_directory("linear_systems")
        try:
            os.mkdir(self.dirname)
        except FileExistsError:
            pass

    def __call__(self, dt, newton_it):

        tick = self.flag.tick
        if self.flag.is_on:
            timestep_dirname = (
                self.dirname + f"/it={tick.iteration}_t={tick.time:.3e}_dt={dt:.3e}/"
            )
            dump_dirname = timestep_dirname + f"Newton_{newton_it}/"
            try:
                os.mkdir(timestep_dirname)
            except FileExistsError:
                pass
            try:
                os.mkdir(dump_dirname)
            except FileExistsError:
                pass
            self.dump_function(basename=dump_dirname)


class InterruptTrigger:
    def __init__(self, flag):
        self.flag = flag

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.is_on:
            mpi.master_print(
                f"\nComPASS - Abortion requested at time t > {self.flag.t:.4e} s\n"
            )
            mpi.abort()


class NewtonLogCallback:
    def __init__(self, filename, newton):
        self.filename = to_output_directory(filename)
        self.newton = newton
        with open(self.filename, "w+") as f:
            f.write(
                "Iteration number, Time, Successful timestep, Newton iterations, Linear iterations\n\n"
            )

    def __call__(self, tick):
        with open(self.filename, "a") as f:
            # FIXME : Opening file at every timestep isn't performance-wise ideal
            f.write(
                f"{tick.iteration},   {tick.time:.13e}, {tick.latest_timestep:.13e}, {len(self.newton.lsolver_iterations)}, {self.newton.lsolver_iterations}\n"
            )


def get_callbacks_from_options(newton, tick0):
    """ Possible options : --dump_ls <comma_separated_times>
                                  --> writes linear systems in file in ASCII mode
                           --dump_ls_binary <comma_separated_times>
                                  --> writes linear systems in file in binary mode (not available with Legacy implementation)
                           --kill <time>
                                  --> kills execution
                           --newton_log <filename>
                                  --> Writes a history of the simulation timesteps, newton iterations
                                      and linear iterations in file filename
    example : --dump_ls 0.0,1.5e6 --kill 1.5e6 """

    callbacks = []
    newton_callbacks = []
    linear_system = newton.lsolver.linear_system

    t_dump_raw = options.database["dump_ls"]
    if t_dump_raw is not None:
        t_dump_list = t_dump_raw.split(",")
        for t_dump in t_dump_list:
            t_dump = float(t_dump)
            dump_flag = TimestepFlag(t_dump)
            dump_flag(tick0)
            dump_trigger = DumpLinearSystemTrigger(linear_system.dump_ascii, dump_flag)
            callbacks.append(dump_flag)
            newton_callbacks.append(dump_trigger)

    t_dump_b_raw = options.database["dump_ls_binary"]
    if t_dump_b_raw is not None:
        t_dump_b_list = t_dump_b_raw.split(",")
        for t_dump_b in t_dump_b_list:
            t_dump_b = float(t_dump_b)
            dump_flag_b = TimestepFlag(t_dump_b)
            dump_flag_b(tick0)
            dump_trigger = DumpLinearSystemTrigger(
                linear_system.dump_binary, dump_flag_b
            )
            callbacks.append(dump_flag_b)
            newton_callbacks.append(dump_trigger)

    newton_log_filename = options.database["newton_log"]
    if newton_log_filename is not None:
        newton_log_callback = NewtonLogCallback(newton_log_filename, newton)
        callbacks.append(newton_log_callback)

    t_kill = options.database["kill"]
    if t_kill is not None:
        t_kill = float(t_kill)
        kill_flag = TimestepFlag(t_kill)
        kill_trigger = InterruptTrigger(kill_flag)
        callbacks.append(kill_trigger)

    newton.callbacks += tuple(newton_callbacks)
    return tuple(callbacks)
