#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .. import messages
from .enthalpies import *
from .physical_properties import PhaseProperty
from .physical_properties_wrappers import register_property


def set_molar_enthalpy_functions(simulation, functions=None, *, check_derivatives=True):
    """
    Set the molar enthalpy functions (with and without the derivatives)
    for all phases.
    If no molar enthalpy function is given, use default values.
    To use another function, give the function
    :param functions: (optional) list of PhaseProperty (one by phase)
    :param check_derivatives: (True by default) check the derivatives of the property
    """
    if functions is not None:
        messages.warning("overriding default molar enthalpy functions")
    else:
        functions = get_default_molar_enthalpy_functions(simulation)

    register_c_property_module = PhaseProperty(
        with_derivatives=simulation.register_c_molar_enthalpies_with_derivatives,
        without_derivatives=simulation.register_c_molar_enthalpies_without_derivatives,
    )
    register_property(
        simulation,
        "molar_enthalpy",
        functions,
        register_c_property_module,
        check_derivatives,
    )


def get_default_molar_enthalpy_functions(simulation):
    physics_name = simulation.physics_name()

    if physics_name == "diphasic":
        assert_diphasic_phase_indexes(simulation)
        assert_diphasic_components_indexes(simulation)
        return [
            gas_diphasic_molar_enthalpies,
            liquid_diphasic_molar_enthalpies,
        ]

    elif physics_name == "immiscible2ph":
        assert_diphasic_phase_indexes(simulation)
        assert_immiscible2ph_components_indexes(simulation)
        return [
            gas_immiscible2ph_molar_enthalpies,
            liquid_diphasic_molar_enthalpies,
        ]

    elif physics_name == "water2ph":
        assert_diphasic_phase_indexes(simulation)
        return [gas_water2ph_enthalpies, liquid_water2ph_enthalpies]

    elif physics_name == "linear_water":
        return pure_phase_enthalpies

    elif physics_name == "brine":
        return brine_enthalpies

    else:
        raise "molar enthalpy not implemented for this physics"


def assert_diphasic_phase_indexes(simulation):
    gas_index = simulation.phase_index(simulation.Phase.gas)
    liquid_index = simulation.phase_index(simulation.Phase.liquid)
    # assumed in diphasic enthalpies definition
    assert gas_index == 0
    assert liquid_index == 1
