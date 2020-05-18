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
    def __init__(self, flag, petsc):
        self.flag = flag
        self.petsc = petsc

    def __call__(self, tick):
        # Update flag
        self.flag(tick)
        if self.flag.on:
            print(">" * 30, "Dump linear system")
            # print(n,t, self.flag.n, self.flag.parent_flag.t)
            assert self.petsc is not None
            self.petsc.dump_LinearSystem()


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
    """ Possible options : --dump_ls <time>,<newton_iteration>
                                  --> writes linear system in file in ASCII mode (Matrix and RHS)
                           --kill <time>,<newton_iteration>
                                  --> kills execution
    example : --dump_ls 1.5e6,2 --kill 1.5e6,2 """

    callbacks = []
    petsc = newton.solver_fmk
    dump_ls = get("--dump_ls")
    if dump_ls is not None:
        tdump = float(dump_ls)
        timeFlag = TimestepFlag(tdump)
        dumpTrigger = DumpLinearSystemTrigger(timeFlag, petsc)
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
