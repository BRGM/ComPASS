import sys
import os
from collections import namedtuple
import ComPASS
from ._kernel import get_kernel
from . import mpi
from .runtime import to_output_directory


def get(name, default=None):
    _, *args = sys.argv
    if name not in args:
        return default
    idx = args.index(name) + 1
    assert idx < len(args), "no value provided after option"
    return args[idx]


def get_bool(name, default=False):
    _, *args = sys.argv
    if name not in args:
        return default
    else:
        return True


class Flag:
    def __init__(self):
        self._status = False

    @property
    def on(self):
        return self._status


class AlwaysOnFlag(Flag):
    def __init__(self):
        super().__init__()
        self._status = True


class TimestepFlag(Flag):
    def __init__(self, t):
        super().__init__()
        self.t = t

    def __call__(self, tick):
        if tick.time >= self.t:
            self._status = True


class DumpLinearSystemTrigger:
    def __init__(self, dump_function, flag):
        self.flag = flag
        self.dump_function = dump_function
        self.dirname = to_output_directory("linear_systems")
        try:
            os.mkdir(self.dirname)
        except FileExistsError:
            pass

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.on:
            dump_dirname = self.dirname + f"/it={tick.iteration}_t={tick.time:.3e}/"
            try:
                os.mkdir(dump_dirname)
            except FileExistsError:
                pass
            print(">" * 30, "Dump linear system")
            self.dump_function(basename=dump_dirname)


class InterruptTrigger:
    def __init__(self, flag):
        self.flag = flag

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.on:
            print(">" * 30, "Kill execution")
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


def get_callbacks_from_options(newton):
    """ Possible options : --dump_ls <time>
                                  --> writes linear system in file in ASCII mode
                           --dump_ls_binary <time>
                                  --> writes linear system in file in binary mode (not available with Legacy implementation)
                           --kill <time>
                                  --> kills execution
                           --newton_log <filename>
                                  --> Writes a history of the simulation timesteps, newton iterations
                                      and linear iterations in file filename
    example : --dump_ls 1.5e6 --kill 1.5e6 """

    callbacks = []
    linear_system = newton.lsolver.linear_system
    t_dump = get("--dump_ls")
    if t_dump is not None:
        t_dump = float(t_dump)
        dump_flag = TimestepFlag(t_dump)
        dump_trigger = DumpLinearSystemTrigger(linear_system.dump_ascii, dump_flag)
        callbacks.append(dump_trigger)

    t_dump_b = get("--dump_ls_binary")
    if t_dump_b is not None:
        t_dump_b = float(t_dump_b)
        dump_flag_b = TimestepFlag(t_dump_b)
        dump_trigger = DumpLinearSystemTrigger(linear_system.dump_binary, dump_flag_b)
        callbacks.append(dump_trigger)

    newton_log_filename = get("--newton_log")
    if newton_log_filename is not None:
        newton_log_callback = NewtonLogCallback(newton_log_filename, newton)
        callbacks.append(newton_log_callback)

    t_kill = get("--kill")
    if t_kill is not None:
        t_kill = float(t_kill)
        kill_flag = TimestepFlag(t_kill)
        kill_trigger = InterruptTrigger(kill_flag)
        callbacks.append(kill_trigger)

    return tuple(callbacks)


if __name__ == "__main__":
    print("toto=", get("--toto"))
    print("tutu=", get("--tutu", "pas de tutu"))
    print("titi=", get("--titi", "pas de titi"))
