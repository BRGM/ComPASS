#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
import numpy as np
from .physical_properties import PhaseProperty


def Xalpha(simulation, pressure, temperature, molar_fractions=None):
    return simulation.fluid_mixture.phase_state_type.Xalpha(
        pressure,
        temperature,
        molar_fractions,
    )


def empty_Xalpha(simulation, len=None):
    return simulation.fluid_mixture.phase_state_type.empty_Xalpha(len)


def constant_physical_property(a):
    def cst_property_with_derivatives(X, dfdX):
        dfdX.pressure = 0
        dfdX.temperature = 0
        dfdX.molar_fractions.fill(0)
        return a  # necessary to init with np.double(a) ?

    def cst_property_without_derivatives(X):
        return a  # necessary to init with np.double(a) ?

    return PhaseProperty(
        with_derivatives=cst_property_with_derivatives,
        without_derivatives=cst_property_without_derivatives,
    )
