
from . import base
from . import data
from . import utils
from . import state

from .._kernel import simulation_wrapper


_not_available = object()


class Simulation:
    def __setattr__(self, name, value):
        assert False, "no modification for now !"

    def __getattr__(self, name):
        for src in (state, base, data, utils, simulation_wrapper):
            value = getattr(src, name, _not_available)
            if value is not _not_available:
                return value
        raise AttributeError(f"'Simulation' object has no attribute {name!r}")

    def __dir__(self):
        res = super().__dir__()
        for src in (state, base, data, utils, simulation_wrapper):
            res.extend(dir(src))
        return res


self = Simulation()
