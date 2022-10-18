#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
from .physical_properties import PhaseProperty


def assert_diphasic_components_indexes(simulation):
    air_component = int(simulation.Component.air) - 1
    water_component = int(simulation.Component.water) - 1
    # assumed in diphasic enthalpy definition
    assert air_component == 0
    assert water_component == 1


# diphasic : this function can be called over scalar X only
def gas_diphasic_molar_enthalpy_with_derivatives(X, dfdX):
    M_air = 29.0e-3
    M_H2O = 18.0e-3
    a = 1990.89e3
    b = 190.16e3
    cc = -1.91264e3
    d = 0.2997e3
    CpGas = 1000.0

    Ts = X.temperature / 100.0
    ss = a + b * Ts + cc * Ts**2.0 + d * Ts**3.0
    dssdT = (b + 2.0 * cc * Ts + 3.0 * d * Ts**2.0) / 100.0
    beta_air = CpGas * M_air
    beta_water = M_H2O
    partial_molar_enthalpy = np.array([beta_air * X.temperature, beta_water * ss])
    dpartial_molar_enthalpydT = np.array([beta_air, beta_water * dssdT])

    dfdX.pressure = 0
    dfdX.molar_fractions[...] = partial_molar_enthalpy
    dfdX.temperature = X.molar_fractions.dot(dpartial_molar_enthalpydT)
    return X.molar_fractions.dot(partial_molar_enthalpy)


def gas_diphasic_molar_enthalpy_without_derivatives(X):
    M_air = 29.0e-3
    M_H2O = 18.0e-3
    a = 1990.89e3
    b = 190.16e3
    cc = -1.91264e3
    d = 0.2997e3
    CpGas = 1000.0

    Ts = X.temperature / 100.0
    ss = a + b * Ts + cc * Ts**2.0 + d * Ts**3.0
    beta_air = CpGas * M_air
    beta_water = M_H2O
    partial_molar_enthalpy = np.array([beta_air * X.temperature, beta_water * ss])
    return X.molar_fractions.dot(partial_molar_enthalpy)


gas_diphasic_molar_enthalpies = PhaseProperty(
    with_derivatives=gas_diphasic_molar_enthalpy_with_derivatives,
    without_derivatives=gas_diphasic_molar_enthalpy_without_derivatives,
)


def pure_liquid_water_molar_enthalpy_with_derivatives(X, dfdX):
    M_H2O = 18.0e-3
    a = -14.4319e3
    b = 4.70915e3
    cc = -4.87534
    d = 1.45008e-2
    T0 = 273.0

    TdegC = X.temperature - T0
    ss = a + b * TdegC + cc * (TdegC**2.0) + d * (TdegC**3.0)
    dssdT = b + 2.0 * cc * TdegC + 3.0 * d * (TdegC**2.0)
    #   FIXME: all components have the same contributions

    dfdX.pressure = 0
    dfdX.molar_fractions.fill(0)
    dfdX.temperature = dssdT * M_H2O
    return M_H2O * ss


def pure_liquid_water_molar_enthalpy_without_derivatives(X):
    M_H2O = 18.0e-3
    a = -14.4319e3
    b = 4.70915e3
    cc = -4.87534
    d = 1.45008e-2
    T0 = 273.0

    TdegC = X.temperature - T0
    ss = a + b * TdegC + cc * (TdegC**2.0) + d * (TdegC**3.0)
    #   FIXME: all components have the same contributions
    return M_H2O * ss


liquid_diphasic_molar_enthalpies = PhaseProperty(
    with_derivatives=pure_liquid_water_molar_enthalpy_with_derivatives,
    without_derivatives=pure_liquid_water_molar_enthalpy_without_derivatives,
)


def assert_immiscible2ph_components_indexes(simulation):
    air_component = int(simulation.Component.air) - 1
    water_component = int(simulation.Component.water) - 1
    # assumed in immiscible2ph enthalpy definition
    assert air_component == 0
    assert water_component == 1


# immiscible2ph : this function can be called over scalar X only
# immiscible2ph needs its own function (distinct with diphasic)
# because the partial derivatives wrt the molar fractions are distinct.
def gas_immiscible2ph_molar_enthalpy_with_derivatives(X, dfdX):
    M_air = 29.0e-3
    CpGas = 1000.0

    beta_air = CpGas * M_air
    # beta_water = 0  # the water component cannot be in the gas phase
    partial_molar_enthalpy = np.array([beta_air * X.temperature, 0.0])
    dpartial_molar_enthalpydT = np.array([beta_air, 0.0])

    dfdX.pressure = 0
    dfdX.molar_fractions[...] = partial_molar_enthalpy
    dfdX.temperature = X.molar_fractions.dot(dpartial_molar_enthalpydT)
    return X.molar_fractions.dot(partial_molar_enthalpy)


def gas_immiscible2ph_molar_enthalpy_without_derivatives(X):
    M_air = 29.0e-3
    CpGas = 1000.0

    beta_air = CpGas * M_air
    # beta_water = 0  # the water component cannot be in the gas phase
    partial_molar_enthalpy = np.array([beta_air * X.temperature, 0.0])
    return X.molar_fractions.dot(partial_molar_enthalpy)


gas_immiscible2ph_molar_enthalpies = PhaseProperty(
    with_derivatives=gas_immiscible2ph_molar_enthalpy_with_derivatives,
    without_derivatives=gas_immiscible2ph_molar_enthalpy_without_derivatives,
)


def gas_water2ph_enthalpy_with_derivatives(X, dfdX):
    dfdX.pressure = 0
    dfdX.molar_fractions.fill(0)
    dfdX.temperature = 190.16e1
    return 1990.89e3 + 190.16e1 * X.temperature


def gas_water2ph_enthalpy_without_derivatives(X):
    return 1990.89e3 + 190.16e1 * X.temperature


gas_water2ph_enthalpies = PhaseProperty(
    with_derivatives=gas_water2ph_enthalpy_with_derivatives,
    without_derivatives=gas_water2ph_enthalpy_without_derivatives,
)


# same as pure_liquid_water_enthalpy without *M_H2O
def pure_liquid_water_enthalpy_with_derivatives(X, dfdX):
    a = -14.4319e3
    b = 4.70915e3
    cc = -4.87534
    d = 1.45008e-2
    T0 = 273.0

    TdegC = X.temperature - T0
    ss = a + b * TdegC + cc * (TdegC**2.0) + d * (TdegC**3.0)
    dssdT = b + 2.0 * cc * TdegC + 3.0 * d * (TdegC**2.0)

    dfdX.pressure = 0
    dfdX.molar_fractions.fill(0)
    dfdX.temperature = dssdT
    return ss


def pure_liquid_water_enthalpy_without_derivatives(X):
    a = -14.4319e3
    b = 4.70915e3
    cc = -4.87534
    d = 1.45008e-2
    T0 = 273.0

    TdegC = X.temperature - T0
    return a + b * TdegC + cc * (TdegC**2) + d * (TdegC**3)


liquid_water2ph_enthalpies = PhaseProperty(
    with_derivatives=pure_liquid_water_enthalpy_with_derivatives,
    without_derivatives=pure_liquid_water_enthalpy_without_derivatives,
)


# def assert_salt_component_index(simulation):
#     salt_component = int(simulation.Component.salt) - 1
#     # assumed in brine definition
#     assert salt_component == 0

# FIXME: this should depend on salt concentration and use specific enthalpy
brine_enthalpies = PhaseProperty(
    with_derivatives=pure_liquid_water_enthalpy_with_derivatives,
    without_derivatives=pure_liquid_water_enthalpy_without_derivatives,
)


def build_pure_phase_enthalpy(
    volumetric_heat_capacity=1.0,
):
    def pure_phase_molar_enthalpy_with_derivatives(X, dfdX):
        dfdX.pressure = 0
        dfdX.temperature = volumetric_heat_capacity
        dfdX.molar_fractions.fill(0)
        return volumetric_heat_capacity * X.temperature

    def pure_phase_molar_enthalpy_without_derivatives(X):
        return volumetric_heat_capacity * X.temperature

    return PhaseProperty(
        with_derivatives=pure_phase_molar_enthalpy_with_derivatives,
        without_derivatives=pure_phase_molar_enthalpy_without_derivatives,
    )


# linear_water -> pure_phase
pure_phase_enthalpies = build_pure_phase_enthalpy()
