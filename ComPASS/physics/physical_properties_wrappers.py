#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
import numpy as np
from .utils import property_sanity_check


# register the phase properties in Python and Fortran
# property_name = 'dynamic_viscosity' or 'molar_density'
def register_property(
    simulation,
    property_name,
    property_functions,
    register_c_property_module,
    check_derivatives,
):

    n_phases = simulation.number_of_phases()
    if n_phases == 1 and not isinstance(property_functions, (list, np.ndarray)):
        #  it is possible that property_functions is not a list
        # if there is only one phase in the eos
        property_functions = [property_functions]
    assert (
        len(property_functions) == n_phases
    ), "You must give one property function per phase"

    for i, prop in enumerate(property_functions):
        simulation.fluid_mixture.register(property_name, i, prop)

    if check_derivatives:
        # use the vectorized functions
        property_sanity_check(simulation, property_name)

    register_c_property_module.with_derivatives(
        *(
            p.c_with_derivatives.address
            for p in simulation.fluid_mixture.properties[property_name]
        )
    )
    # necessary if there is a well
    register_c_property_module.without_derivatives(
        *(
            p.c_without_derivatives.address
            for p in simulation.fluid_mixture.properties[property_name]
        )
    )


def build_molar_fractions(simulation, salt_molar_fraction):
    eos_name = simulation.eos_name()
    if salt_molar_fraction is None:
        assert eos_name == "linear_water", "Molar fractions missing in brine eos"
        molar_fractions = None
    else:
        assert eos_name == "brine", "Do not give molar fractions in linear_water eos"
        # in brine eos, only Cs is given amoung the two components
        assert int(simulation.Component.salt) - 1 == 0
        # transpose is necessary when salt_molar_fraction is an array
        molar_fractions = np.array([salt_molar_fraction, 1.0 - salt_molar_fraction]).T
    return molar_fractions


# common function to compute phase property
# property can be "dynamic_viscosity" or "molar_density"
def phase_property(
    simulation,
    property_name,
    phase_id,
    X,
):
    return np.asarray(
        simulation.fluid_mixture.properties[property_name][
            phase_id
        ].vect_without_derivatives(X)
    )


# brine : salt_molar_fraction is a scalar
# linear_water : salt_molar_fraction is absent
def molar_density(simulation, pressure, temperature, salt_molar_fraction=None):
    iph = simulation.phase_index(simulation.Phase.single_phase)
    molar_fractions = build_molar_fractions(simulation, salt_molar_fraction)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)

    return phase_property(simulation, "molar_density", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def gas_molar_density(simulation, pressure, temperature, molar_fractions=None):
    iph = simulation.phase_index(simulation.Phase.gas)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "molar_density", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def liquid_molar_density(simulation, pressure, temperature, molar_fractions=None):
    iph = simulation.phase_index(simulation.Phase.liquid)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "molar_density", iph, X)


# brine : salt_molar_fraction is a scalar
# linear_water : salt_molar_fraction is absent
def volumetric_mass_density(
    simulation, pressure, temperature, salt_molar_fraction=None
):
    iph = simulation.phase_index(simulation.Phase.single_phase)
    molar_fractions = build_molar_fractions(simulation, salt_molar_fraction)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "volumetric_mass_density", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def gas_volumetric_mass_density(
    simulation, pressure, temperature, molar_fractions=None
):
    iph = simulation.phase_index(simulation.Phase.gas)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "volumetric_mass_density", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def liquid_volumetric_mass_density(
    simulation, pressure, temperature, molar_fractions=None
):
    iph = simulation.phase_index(simulation.Phase.liquid)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "volumetric_mass_density", iph, X)


# brine : salt_molar_fraction is a scalar
# linear_water : salt_molar_fraction is missing
def dynamic_viscosity(simulation, pressure, temperature, salt_molar_fraction=None):
    iph = simulation.phase_index(simulation.Phase.single_phase)
    molar_fractions = build_molar_fractions(simulation, salt_molar_fraction)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)

    return phase_property(simulation, "dynamic_viscosity", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def gas_dynamic_viscosity(simulation, pressure, temperature, molar_fractions=None):
    iph = simulation.phase_index(simulation.Phase.gas)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "dynamic_viscosity", iph, X)


# diphasic, immiscible2ph : molar_fractions is a vector
# water2ph : molar_fractions is absent
def liquid_dynamic_viscosity(simulation, pressure, temperature, molar_fractions=None):
    iph = simulation.phase_index(simulation.Phase.liquid)
    X = simulation.Xalpha(pressure, temperature, molar_fractions)
    return phase_property(simulation, "dynamic_viscosity", iph, X)
