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


# diphasic
gas_diphasic_viscosities = constant_physical_property(2.0e-5)
liquid_diphasic_viscosities = constant_physical_property(1.0e-3)
# linear_water -> pure_phase
pure_phase_viscosities = constant_physical_property(1.0)


def gas_water2ph_viscosity_with_derivatives(X, dfdX):
    dfdX.pressure = 0
    dfdX.molar_fractions.fill(0)
    dfdX.temperature = 0.361e-7
    return (0.361 * X.temperature - 10.2) * 1.0e-7


def gas_water2ph_viscosity_without_derivatives(X):
    return (0.361 * X.temperature - 10.2) * 1.0e-7


gas_water2ph_viscosities = PhaseProperty(
    with_derivatives=gas_water2ph_viscosity_with_derivatives,
    without_derivatives=gas_water2ph_viscosity_without_derivatives,
)


def liquid_water2ph_viscosity_with_derivatives(X, dfdX):
    dfdX.pressure = 0
    dfdX.molar_fractions.fill(0)
    Tref = X.temperature - 273 - 8.435
    b = np.sqrt(8078.4 + Tref**2)
    a = 0.021482 * (Tref + b) - 1.2
    da = 0.021482 * (1.0 + Tref / b)
    dfdX.temperature = -1.0e-3 * da / (a**2)
    return 1.0e-3 / a


def liquid_water2ph_viscosity_without_derivatives(X):
    Tref = X.temperature - 273 - 8.435
    b = np.sqrt(8078.4 + Tref**2)
    a = 0.021482 * (Tref + b) - 1.2
    return 1.0e-3 / a


liquid_water2ph_viscosities = PhaseProperty(
    with_derivatives=liquid_water2ph_viscosity_with_derivatives,
    without_derivatives=liquid_water2ph_viscosity_without_derivatives,
)


def assert_salt_component_index(simulation):
    salt_component = int(simulation.Component.salt) - 1
    # assumed in brine definition
    assert salt_component == 0


def brine_viscosity_with_derivatives(X, dfdX):
    Tref = X.temperature - 273 - 8.435
    b = np.sqrt(8078.4 + Tref**2)
    # ns : no salt
    ns = 0.021482 * (Tref + b) - 1.2
    dnsdT = 0.021482 * (1.0 + Tref / b)
    # sc : salt correction, C[0] is C[salt]
    sc = 1.0 + X.molar_fractions[0] * 1.34 + 6.12 * X.molar_fractions[0] ** 2
    dscdCs = 1.34 + 2 * 6.12 * X.molar_fractions[0]
    dmudCs = 1.0e-3 * dscdCs / ns

    dfdX.pressure = 0
    dfdX.temperature = -1.0e-3 * sc * dnsdT / (ns**2)
    dfdX.molar_fractions[...] = np.array([dmudCs, 0.0])
    return 1.0e-3 * sc / ns


def brine_viscosity_without_derivatives(X):
    Tref = X.temperature - 273 - 8.435
    b = np.sqrt(8078.4 + Tref**2)
    # ns : no salt
    ns = 0.021482 * (Tref + b) - 1.2
    # sc : salt correction, C[0] is C[salt]
    sc = 1.0 + X.molar_fractions[0] * 1.34 + 6.12 * X.molar_fractions[0] ** 2
    return 1.0e-3 * sc / ns


brine_viscosities = PhaseProperty(
    with_derivatives=brine_viscosity_with_derivatives,
    without_derivatives=brine_viscosity_without_derivatives,
)
