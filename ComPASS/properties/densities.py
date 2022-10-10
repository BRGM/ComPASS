#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
from .physical_properties import PhaseProperty
from .utils import constant_physical_property


def gas_diphasic_molar_density_with_derivatives(X, dfdX):
    Rgp = 8.314
    one_over_RgpT = 1.0 / (Rgp * X.temperature)
    dfdX.pressure = one_over_RgpT
    dfdX.temperature = -X.pressure * one_over_RgpT / X.temperature
    dfdX.molar_fractions.fill(0)
    return X.pressure * one_over_RgpT


def gas_diphasic_molar_density_without_derivatives(X):
    Rgp = 8.314
    return X.pressure / (Rgp * X.temperature)


gas_diphasic_molar_densities = PhaseProperty(
    with_derivatives=gas_diphasic_molar_density_with_derivatives,
    without_derivatives=gas_diphasic_molar_density_without_derivatives,
)
# M_H2O = 18.0e-3
liquid_diphasic_molar_densities = constant_physical_property(1000.0 / 18.0e-3)


def assert_diphasic_components_indexes(simulation):
    air_component = int(simulation.Component.air) - 1
    water_component = int(simulation.Component.water) - 1
    # assumed in diphasic definition
    assert air_component == 0
    assert water_component == 1


# M_H2O = 18.0e-3
# M_air = 29.0e-3
diphasic_components_molar_mass = [29.0e-3, 18.0e-3]


def gas_water2ph_volumetric_mass_density_with_derivatives(X, dfdX):
    dfdX.molar_fractions.fill(0)
    # Z = 1.0  # CHECKME: ?
    M_H2O = 0.018016
    R = 8.3145
    # M_over_RTZ = M_H2O / (R * X.temperature * Z)
    M_over_RTZ = M_H2O / (R * X.temperature)
    dfdX.pressure = M_over_RTZ
    dfdX.temperature = -X.pressure * M_over_RTZ / X.temperature
    return X.pressure * M_over_RTZ


def gas_water2ph_volumetric_mass_density_without_derivatives(X):
    # Z = 1.0  # CHECKME: ?
    M_H2O = 0.018016
    R = 8.3145
    # return X.pressure * M_H2O / (R * X.temperature * Z)
    return X.pressure * M_H2O / (R * X.temperature)


# in ComPASS water2ph physics, f_MolarDensity is mass density, see #348
gas_water2ph_volumetric_mass_densities = PhaseProperty(
    with_derivatives=gas_water2ph_volumetric_mass_density_with_derivatives,
    without_derivatives=gas_water2ph_volumetric_mass_density_without_derivatives,
)


def liquid_water2ph_volumetric_mass_density_with_derivatives(X, dfdX):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23

    Temperature_square = X.temperature**2
    ss = rho0 + a * X.temperature + b * Temperature_square
    dssdT = a + 2 * b * X.temperature
    cw = (
        a1
        + a2 * X.pressure
        + X.temperature * (b1 + b2 * X.pressure)
        + Temperature_square * (c1 + c2 * X.pressure)
    )
    dcwdp = a2 + b2 * X.temperature + c2 * Temperature_square
    dcwdT = b1 + b2 * X.pressure + 2 * X.temperature * (c1 + c2 * X.pressure)
    psat = 1.0e-3 * (X.temperature - 273.0) ** 4.0  # use extern function psat ???
    dpsatdT = 4.0e-3 * (X.temperature - 273.0) ** 3.0
    p_rel = X.pressure - psat
    dfdX.molar_fractions.fill(0)
    dfdX.pressure = ss * (dcwdp * p_rel + cw)
    dfdX.temperature = dssdT * (1.0 + cw * p_rel) + ss * (dcwdT * p_rel - cw * dpsatdT)
    return ss * (1.0 + cw * p_rel)


def liquid_water2ph_volumetric_mass_density_without_derivatives(X):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    Temperature_square = X.temperature**2
    ss = rho0 + a * X.temperature + b * Temperature_square
    cw = (
        a1
        + a2 * X.pressure
        + X.temperature * (b1 + b2 * X.pressure)
        + Temperature_square * (c1 + c2 * X.pressure)
    )
    psat = 1.0e-3 * (X.temperature - 273.0) ** 4.0  # use extern function psat ???
    p_rel = X.pressure - psat
    return ss * (1.0 + cw * p_rel)


# in ComPASS water2ph physics, f_MolarDensity is mass density, see #348
liquid_water2ph_volumetric_mass_densities = PhaseProperty(
    with_derivatives=liquid_water2ph_volumetric_mass_density_with_derivatives,
    without_derivatives=liquid_water2ph_volumetric_mass_density_without_derivatives,
)


def assert_salt_component_index(simulation):
    salt_component = int(simulation.Component.salt) - 1
    # assumed in brine definition
    assert salt_component == 0


def brine_volumetric_mass_density_with_derivatives(X, dfdX):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    # rs = 0.0 # CHECKME: What is this???? Link to capillary pressure?
    cwrs = 1.0  # + 5.e-2*rs
    cw = cwrs * (
        a1
        + a2 * X.pressure
        + X.temperature * (b1 + b2 * X.pressure)
        + X.temperature**2 * (c1 + c2 * X.pressure)
    )
    dcwdp = cwrs * (a2 + b2 * X.temperature + c2 * X.temperature**2)
    dcwdT = cwrs * (
        (b1 + b2 * X.pressure) + X.temperature * 2.0 * (c1 + c2 * X.pressure)
    )
    # Salt correction
    sc = (rho0 + a * X.temperature + b * X.temperature**2) * (
        1.0 + 6.51e-4 * X.molar_fractions[0]
    )  # 0 is salt ???
    dscdT = (a + 2.0 * b * X.temperature) * (
        1.0 + 6.51e-4 * X.molar_fractions[0]
    )  # 0 is salt ???
    dscdCs = 6.51e-4 * (rho0 + a * X.temperature + b * X.temperature**2)
    # psat of pure water
    psat = (X.temperature - 273.0) ** 4.0 / 1.0e3
    dpsatdT = 4.0 * (X.temperature - 273.0) ** 3.0 / 1.0e3
    p_rel = X.pressure - psat

    dfdX.pressure = sc * dcwdp * p_rel + sc * cw
    dfdX.temperature = dscdT * (1.0 + cw * p_rel) + sc * (dcwdT * p_rel - cw * dpsatdT)
    dfdCs = dscdCs * (1.0 + cw * p_rel)
    dfdX.molar_fractions[...] = np.array([dfdCs, 0.0])
    return sc * (1.0 + cw * p_rel)


def brine_volumetric_mass_density_without_derivatives(X):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    # rs = 0.0 # CHECKME: What is this???? Link to capillary pressure?
    cwrs = 1.0  # + 5.e-2*rs
    # Salt correction
    sc = (rho0 + a * X.temperature + b * X.temperature**2) * (
        1.0 + 6.51e-4 * X.molar_fractions[0]
    )  # 0 is salt ???
    cw = cwrs * (
        a1
        + a2 * X.pressure
        + X.temperature * (b1 + b2 * X.pressure)
        + X.temperature**2 * (c1 + c2 * X.pressure)
    )
    psat = (X.temperature - 273.0) ** 4.0 / 1.0e3
    p_rel = X.pressure - psat
    return sc * (1.0 + cw * p_rel)


# in ComPASS brine physics, f_MolarDensity is mass density, see #51
brine_volumetric_mass_densities = PhaseProperty(
    with_derivatives=brine_volumetric_mass_density_with_derivatives,
    without_derivatives=brine_volumetric_mass_density_without_derivatives,
)


def build_pure_phase_volumetric_mass_density(
    specific_mass=1.0,
    compressibility=0.0,
    thermal_expansivity=0.0,
    reference_pressure=1.0e5,  # 1 bar
    reference_temperature=293.15,  # 20°C,
):
    def pure_phase_volumetric_mass_density_with_derivatives(X, dfdX):
        val = (
            specific_mass
            * np.exp(compressibility * (X.pressure - reference_pressure))
            * np.exp(thermal_expansivity * (X.temperature - reference_temperature))
        )
        dfdX.pressure = compressibility * val
        dfdX.temperature = thermal_expansivity * val
        dfdX.molar_fractions.fill(0)
        return val

    def pure_phase_volumetric_mass_density_without_derivatives(X):
        return (
            specific_mass
            * np.exp(compressibility * (X.pressure - reference_pressure))
            * np.exp(thermal_expansivity * (X.temperature - reference_temperature))
        )

    return PhaseProperty(
        with_derivatives=pure_phase_volumetric_mass_density_with_derivatives,
        without_derivatives=pure_phase_volumetric_mass_density_without_derivatives,
    )


# in ComPASS linear_water physics, f_MolarDensity is mass density, see #51
# linear_water -> pure_phase
pure_phase_volumetric_mass_densities = build_pure_phase_volumetric_mass_density()
