import sys
import os
from collections import namedtuple
import ComPASS
from ._kernel import get_kernel


def get(name, default=None):
    _, *args = sys.argv
    if name not in args:
        return default
    idx = args.index(name) + 1
    assert idx < len(args), "no value provided after option"
    return args[idx]


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
    def __init__(self, dump_function, flag, binary=False):
        self.flag = flag
        self.dump_function = dump_function
        self.binary = binary

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.on:
            print(">" * 30, "Dump linear system")
            self.dump_function(basename="t=%e_" % tick.time, binary=self.binary)


class InterruptTrigger:
    def __init__(self, flag):
        self.flag = flag

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.on:
            print(">" * 30, "Kill execution")
            1 / 0


def get_callbacks_from_options(newton):
    """ Possible options : --dump_ls <time>
                                  --> writes linear system in file in ASCII mode
                           --dump_ls_binary <time>
                                  --> writes linear system in file in binary mode (not available with Legacy implementation)
                           --kill <time>
                                  --> kills execution
    example : --dump_ls 1.5e6 --kill 1.5e6 """

    callbacks = []
    lsolver = newton.lsolver
    dump_ls = get("--dump_ls")
    if dump_ls is not None:
        tdump = float(dump_ls)
        timeFlag = TimestepFlag(tdump)
        dumpTrigger = DumpLinearSystemTrigger(lsolver.dump_system, timeFlag)
        callbacks.append(dumpTrigger)

    dump_ls = get("--dump_ls_binary")
    if dump_ls is not None:
        tdump = float(dump_ls)
        timeFlag = TimestepFlag(tdump)
        dumpTrigger = DumpLinearSystemTrigger(
            lsolver.dump_system, timeFlag, binary=True
        )
        callbacks.append(dumpTrigger)

    kill = get("--kill")
    if kill is not None:
        tkill = float(kill)
        timeFlag = TimestepFlag(tkill)
        kill_trigger = InterruptTrigger(timeFlag)
        callbacks.append(kill_trigger)
    return tuple(callbacks)


if __name__ == "__main__":
    print("toto=", get("--toto"))
    print("tutu=", get("--tutu", "pas de tutu"))
    print("titi=", get("--titi", "pas de titi"))
