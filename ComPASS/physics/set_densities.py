#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .. import messages
from .densities import *
from .physical_properties import PhaseProperty
from .physical_properties_wrappers import register_property


def set_molar_density_functions(simulation, functions=None, *, check_derivatives=True):
    """
    Set the molar density functions (with and without the derivatives)
    for all phases.
    If no molar density function is given, use default values.
    To use another function, give the function
    :param functions: (optional) list of PhaseProperty (one by phase)
    :param check_derivatives: (True by default) check the derivatives of the property
    """
    if functions is not None:
        messages.warning("overriding default molar density functions")
    else:
        functions = get_default_molar_density_functions(simulation)

    register_c_property_module = PhaseProperty(
        with_derivatives=simulation.register_c_molar_densities_with_derivatives,
        without_derivatives=simulation.register_c_molar_densities_without_derivatives,
    )
    register_property(
        simulation,
        "molar_density",
        functions,
        register_c_property_module,
        check_derivatives,
    )


def get_default_molar_density_functions(simulation):
    eos_name = simulation.eos_name()

    if eos_name == "diphasic" or eos_name == "immiscible2ph":
        assert_diphasic_phase_indexes(simulation)
        return [
            gas_diphasic_molar_densities,
            liquid_diphasic_molar_densities,
        ]

    elif eos_name == "water2ph":
        assert_diphasic_phase_indexes(simulation)
        # in ComPASS water2ph eos, f_MolarDensity is mass density, see #348
        return [
            gas_water2ph_volumetric_mass_densities,
            liquid_water2ph_volumetric_mass_densities,
        ]

    elif eos_name == "linear_water":
        # in ComPASS linear_water eos, f_MolarDensity is mass density, see #51
        return pure_phase_volumetric_mass_densities

    elif eos_name == "brine":
        # in ComPASS brine eos, f_MolarDensity is mass density, see #51
        assert_salt_component_index(simulation)
        return brine_volumetric_mass_densities

    else:
        raise "molar density not implemented for this eos"


def assert_diphasic_phase_indexes(simulation):
    gas_index = simulation.phase_index(simulation.Phase.gas)
    liquid_index = simulation.phase_index(simulation.Phase.liquid)
    # assumed in diphasic viscosities definition
    assert gas_index == 0
    assert liquid_index == 1
