#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .. import messages
from .viscosities import *
from .physical_properties import PhaseProperty
from .physical_properties_wrappers import register_property


def set_viscosity_functions(simulation, functions=None):
    """
    Set the viscosity functions (with and without the derivatives)
    for all phases.
    If no viscosity function is given, use default values.
    To use another function, give the function
    :param functions: (optional) list of PhaseProperty (one by phase)
    """
    if functions is not None:
        messages.warning("overriding default viscosity functions")
    else:
        functions = get_default_viscosity_functions(simulation)

    register_c_property_module = PhaseProperty(
        with_derivatives=simulation.register_c_viscosities_with_derivatives,
        without_derivatives=simulation.register_c_viscosities_without_derivatives,
    )
    register_property(
        simulation,
        "dynamic_viscosity",
        functions,
        register_c_property_module,
    )


def get_default_viscosity_functions(simulation):
    eos_name = simulation.eos_name()

    if eos_name == "diphasic" or eos_name == "immiscible2ph":
        assert_diphasic_phase_indexes(simulation)
        return [
            gas_diphasic_viscosities,
            liquid_diphasic_viscosities,
        ]

    elif eos_name == "water2ph":
        assert_diphasic_phase_indexes(simulation)
        return [gas_water2ph_viscosities, liquid_water2ph_viscosities]

    elif eos_name == "linear_water":
        return pure_phase_viscosities

    elif eos_name == "brine":
        assert_salt_component_index(simulation)
        return brine_viscosities

    else:
        raise "viscosity not implemented for this eos"


def assert_diphasic_phase_indexes(simulation):
    gas_index = simulation.phase_index(simulation.Phase.gas)
    liquid_index = simulation.phase_index(simulation.Phase.liquid)
    # assumed in diphasic viscosities definition
    assert gas_index == 0
    assert liquid_index == 1
