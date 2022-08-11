#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from numba import jit
from .. import messages
from .densities import *
from .physical_properties import PhaseProperty
from .physical_properties_wrappers import register_property


def set_molar_density_functions(
    simulation,
    functions=None,
    *,
    check_derivatives=True,
    _update_volumetric_mass_density_functions=True,
):
    """
    Set the molar density functions (with and without the derivatives)
    for all phases.
    If no molar density function is given, use default values.
    To use another function, give the function
    :param functions: (optional) list of PhaseProperty (one by phase)
    :param check_derivatives: (True by default) check the derivatives of the property
    :param _update_volumetric_mass_density_functions: (True by default) update automatically the volumetric mass density functions
    """
    if functions is not None:
        if _update_volumetric_mass_density_functions:
            messages.warning(
                "overriding default molar and volumetric mass density functions"
            )
        else:
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
    if _update_volumetric_mass_density_functions:
        _set_volumetric_mass_density_functions(
            simulation,
            check_derivatives=check_derivatives,
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


def _set_volumetric_mass_density_functions(
    simulation,
    *,
    check_derivatives=True,
):
    """
    Set the volumetric mass density functions (with and without the derivatives)
    for all phases with volumetric_mass_density = (Sum Ci Mi)*molar_density
    :param check_derivatives: (True by default) check the derivatives of the property
    """
    functions = build_volumetric_mass_density_from_molar_density(
        simulation,
        simulation.fluid_mixture.components_molar_mass,
    )

    register_c_property_module = PhaseProperty(
        with_derivatives=simulation.register_c_volumetric_mass_densities_with_derivatives,
        without_derivatives=simulation.register_c_volumetric_mass_densities_without_derivatives,
    )
    register_property(
        simulation,
        "volumetric_mass_density",
        functions,
        register_c_property_module,
        check_derivatives,
    )


def get_default_components_molar_mass(simulation):
    eos_name = simulation.eos_name()

    if eos_name == "diphasic" or eos_name == "immiscible2ph":
        # components indexes guessed in diphasic_components_molar_mass
        assert_diphasic_components_indexes(simulation)
        return diphasic_components_molar_mass

    elif eos_name == "water2ph":
        # in ComPASS water2ph eos, f_MolarDensity is mass density, see #348
        return 1.0

    elif eos_name == "linear_water":
        # in ComPASS linear_water eos, f_MolarDensity is mass density, see #51
        return 1.0

    elif eos_name == "brine":
        # in ComPASS brine eos, f_MolarDensity is mass density, see #51
        return [1.0, 1.0]

    else:
        raise "no default components molar mass for this eos"


def set_components_molar_mass(simulation, components_molar_mass=None):
    if components_molar_mass is None:
        components_molar_mass = get_default_components_molar_mass(simulation)
    else:
        messages.warning(
            "overriding default components molar mass and mass density functions"
        )

    simulation.fluid_mixture.set_components_molar_mass(components_molar_mass)
    _set_volumetric_mass_density_functions(simulation)


def build_volumetric_mass_density_from_molar_density(
    simulation,
    components_molar_mass,  # np.array of size n_comp
):
    volumetric_mass_density_functions = []
    for phase_property in simulation.fluid_mixture.properties["molar_density"]:
        volumetric_mass_density_functions.append(
            build_phase_volumetric_mass_density_from_molar_density(
                components_molar_mass,
                PhaseProperty(
                    with_derivatives=phase_property.py_with_derivatives,
                    without_derivatives=phase_property.py_without_derivatives,
                ),
            )
        )

    return volumetric_mass_density_functions


def build_phase_volumetric_mass_density_from_molar_density(
    components_molar_mass,  # np.array of size n_comp
    phase_molar_density_functions,  # PhaseProperty
):
    # is it necessary to compile molar_density first ?
    compiled_phase_molar_densities_with_derivatives = jit(
        phase_molar_density_functions.with_derivatives, nopython=True
    )
    compiled_phase_molar_densities_without_derivatives = jit(
        phase_molar_density_functions.without_derivatives, nopython=True
    )

    def volumetric_mass_density_with_derivatives(X, dfdX):
        phase_molar_density = compiled_phase_molar_densities_with_derivatives(X, dfdX)
        # I removed MCP from the formula
        # MCP(AIR_COMP, iph)*M_air*C(AIR_COMP) + MCP(WATER_COMP, iph)*M_H2O*C(WATER_COMP)
        weighted_components_molar_mass = X.molar_fractions.dot(components_molar_mass)
        dfdX.pressure *= weighted_components_molar_mass
        dfdX.temperature *= weighted_components_molar_mass
        # is isscalar quicker than empty list + append + conversion to np.array ?
        if np.isscalar(phase_molar_density):
            weighted_phase_molar_density = phase_molar_density * components_molar_mass
        else:
            weighted_phase_molar_density = []
            for Mcomp in components_molar_mass:
                weighted_phase_molar_density.append(Mcomp * phase_molar_density)
            weighted_phase_molar_density = np.array(weighted_phase_molar_density).T
        dfdX.molar_fractions *= weighted_components_molar_mass
        dfdX.molar_fractions += weighted_phase_molar_density
        return weighted_components_molar_mass * phase_molar_density

    def volumetric_mass_density_without_derivatives(X):
        phase_molar_density = compiled_phase_molar_densities_without_derivatives(X)
        # I removed MCP from the formula
        # MCP(AIR_COMP, iph)*M_air*C(AIR_COMP) + MCP(WATER_COMP, iph)*M_H2O*C(WATER_COMP)
        weighted_components_molar_mass = X.molar_fractions.dot(components_molar_mass)
        return weighted_components_molar_mass * phase_molar_density

    return PhaseProperty(
        with_derivatives=volumetric_mass_density_with_derivatives,
        without_derivatives=volumetric_mass_density_without_derivatives,
    )


def assert_diphasic_phase_indexes(simulation):
    gas_index = simulation.phase_index(simulation.Phase.gas)
    liquid_index = simulation.phase_index(simulation.Phase.liquid)
    # assumed in diphasic viscosities definition
    assert gas_index == 0
    assert liquid_index == 1
