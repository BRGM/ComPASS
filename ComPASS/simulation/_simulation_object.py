from functools import partial, update_wrapper
from . import AlignmentMethod
from . import fake_methods
from . import base
from . import data
from . import utils
from .. import mpi

from .._kernel import simulation_wrapper
from ..wells.wells import get_wellhead
from ..wells.connections import WellDataConnections, add_well_connections
from ..properties.physical_properties import FluidMixtureProperties

_not_available = object()
_fake_members = base, data, utils, simulation_wrapper


class SimulationInfo:
    def __init__(self):
        self.system = None
        self.ghosts_synchronizer = None


class SimmulationBase:
    """
    A temporary object to group true members for simulation objects.
    """

    def __init__(self, model, well_data_provider):
        self.info = SimulationInfo()
        # FIXME: should be encapsulated elsewhere
        self.alignment = AlignmentMethod.inverse_diagonal
        self.initialized = False
        self.mesh_is_local = False
        self.well_data_provider = well_data_provider
        self.well_connections = WellDataConnections()
        self.well_model = None
        self.unknown_producers_density = True
        self.scheme = None
        # FIXME: can number of components be different in each phase?
        np, nc = model  # number of components, number of phases
        self.fluid_mixture = FluidMixtureProperties(np, nc)

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

    @property
    def is_sequential(self):
        return mpi.communicator().size == 1


class Simulation:
    already_instanciated = False

    def __init__(self, kernel=None):
        assert not Simulation.already_instanciated, "Simulation must be a singleton"
        Simulation.already_instanciated = True
        # a copy of the python object is made on the C++ side
        # so that the simulation will not be garbage collected
        assert kernel is not None
        kernel.register_simulation(self)

        # self.base = SimulationBase() # is not allowed because __setattr__ is prohibited
        well_data_provider = partial(get_wellhead, self)
        model = (kernel.number_of_phases(), kernel.number_of_components())
        self.__dict__["base"] = SimmulationBase(model, well_data_provider)
        self.set_viscosity_functions()
        self.set_molar_density_functions(
            _update_volumetric_mass_density_functions=False
        )
        # also updates the volumetric mass density functions
        self.set_components_molar_mass()
        self.set_molar_enthalpy_functions()

    def __setattr__(self, name, value):
        setattr(self.__dict__["base"], name, value)

    def __getattr__(self, name):
        value = getattr(fake_methods, name, _not_available)
        if value is not _not_available:
            return update_wrapper(partial(value, self), value)
        for src in _fake_members:
            value = getattr(src, name, _not_available)
            if value is not _not_available:
                return value
        try:
            return getattr(self.__dict__["base"], name)
        except AttributeError:
            pass
        raise AttributeError(f"'Simulation' object has no attribute {name!r}")

    def __dir__(self):
        res = super().__dir__()
        res.extend(dir(fake_methods))
        for src in _fake_members:
            res.extend(dir(src))
        return res
