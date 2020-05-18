from functools import partial
from . import fake_methods
from . import base
from . import data
from . import utils
from . import state

from .._kernel import simulation_wrapper
from ..wells.wells import get_well_data
from ..wells.connections import WellDataConnections, add_well_connections

_not_available = object()
_fake_members = state, base, data, utils, simulation_wrapper


class SimmulationBase:
    """
        A temporary object to group true members for simulation objects.
    """

    def __init__(self, well_data_provider):
        self.well_data_provider = well_data_provider
        self.well_connections = WellDataConnections()

    def add_well_connections(self, well_pairs=None, proc_requests=None):
        """
            Given a sequence of well ids pairs `(source, target)`,
            create connections between simulations domains so that domains
            where the target exists can acces flowrate informations from the source
            well (wherever they are).

            The well informations are explicitely synchronised by calling
            `self.well_connections.synchronize()`.

            The source information can be access using the `__getitem__` method
            on self.well_connections : e.g. `self.well_connections[source].mass_flowrate`.

            :param well_pairs: a sequence of well ids pairs `(source, target)`
            :param proc_requests: a sequence of pair (proc, list of wells to make available)
        """
        add_well_connections(
            self.well_connections,
            self.well_data_provider,
            well_pairs=well_pairs,
            proc_requests=proc_requests,
        )


class Simulation:
    def __init__(self):
        # self.base = SimulationBase() # is not allowed because __setattr__ is prohibited
        well_data_provider = partial(get_well_data, self)
        self.__dict__["base"] = SimmulationBase(well_data_provider)

    def __setattr__(self, name, value):
        assert False, "no modification for now !"

    def __getattr__(self, name):
        try:
            return getattr(self.__dict__["base"], name)
        except AttributeError:
            pass
        value = getattr(fake_methods, name, _not_available)
        if value is not _not_available:
            return partial(value, self)
        for src in _fake_members:
            value = getattr(src, name, _not_available)
            if value is not _not_available:
                return value
        raise AttributeError(f"'Simulation' object has no attribute {name!r}")

    def __dir__(self):
        res = super().__dir__()
        res.extend(dir(fake_methods))
        for src in _fake_members:
            res.extend(dir(src))
        return res


self = Simulation()
