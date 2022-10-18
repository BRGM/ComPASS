#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
import numpy as np
from ComPASS.utils.units import bar, MPa
from ComPASS import messages
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


def check_phase_property_derivatives(
    simulation,
    properties,
    Xalpha,
    dh=1e-6,
    atol=1e-6,
    rtol=1e-5,
):
    ncomp = simulation.number_of_components()
    # properties contains one property by phase
    for property in properties:
        # property contains one function with the derivatives computation
        # and one function without the derivatives computation
        f_without_d = property.without_derivatives  # function val = f(X)
        f_with_d = property.with_derivatives  # function val = f(X, dfdX)

        for fieldname in Xalpha.dtype.names:
            # masks allows to check molar fractions (array) in the same way
            # as pressure and temperature
            masks = [1]
            if fieldname == "molar_fractions":
                masks = np.identity(ncomp)
            for mask in masks:
                Xplus = Xalpha.copy()
                Xplus[fieldname] += dh * mask  # might be out of function definition
                Xminus = Xalpha.copy()
                Xminus[fieldname] -= dh * mask  # might be out of function definition
                dfdX = simulation.empty_Xalpha(np.size(Xalpha))
                finite_diff = (f_without_d(Xplus) - f_without_d(Xminus)) / 2 / dh
                f_with_d(Xalpha, dfdX)  # fill dfdX
                # extract the good column
                function_derivatives = dfdX[fieldname].dot(mask)
                assert np.allclose(
                    finite_diff, function_derivatives, atol=atol, rtol=rtol
                ), "derivatives are wrong"
    return True


# check that the value computed with the derivatives val = property(X, dfdX)
# is equal to the one computed without the derivatives val = property(X)
def check_phase_property_coherent_values(
    simulation,
    properties,
    Xalpha,
    atol=1e-6,
    rtol=1e-5,
    output_maxerror=False,
):
    # properties contains one property by phase
    for property in properties:
        # property contains one function with derivatives computation
        # and one function without the derivatives computation
        f_with_d = property.with_derivatives  # function val = f(X, dfdX)
        f_without_d = property.without_derivatives  # function val = f(X)
        val_without_d = f_without_d(Xalpha)
        dfdX = simulation.empty_Xalpha(np.size(Xalpha))
        val_with_d = f_with_d(Xalpha, dfdX)

        if output_maxerror:
            print(np.fabs(val_with_d - val_without_d).max())
    return np.allclose(val_with_d, val_without_d, atol=atol, rtol=rtol)


def property_sanity_check(
    simulation,
    property_name,
    pressure=np.logspace(np.log10(bar), np.log10(50 * MPa), 100),
    temperature=np.linspace(273.15, 573.15, 100),
    molar_fractions=None,
):

    if molar_fractions is None:
        array_size = max(np.size(pressure), np.size(temperature))
        molar_fractions = default_linspace_molar_fractions(
            array_size,
            simulation.number_of_components(),
        )
    Xalpha = simulation.Xalpha(pressure, temperature, molar_fractions)

    fluid_mixture_properties = simulation.fluid_mixture.properties[property_name]

    # properties contains one PhaseProperty by phase, each phase contains
    # the vectorize property with and without the derivatives
    properties = [
        PhaseProperty(
            with_derivatives=fluid_mixture_properties[i].vect_with_derivatives,
            without_derivatives=fluid_mixture_properties[i].vect_without_derivatives,
        )
        for i in range(simulation.number_of_phases())
    ]

    assert check_phase_property_coherent_values(
        simulation, properties, Xalpha
    ), "Property values differ with and without the derivatives computation"
    assert check_phase_property_derivatives(simulation, properties, Xalpha)


def default_linspace_molar_fractions(npts, ncomp):
    if ncomp == 1:
        return np.array([np.linspace(0.0, 1.0, npts)]).T
    elif ncomp == 2:
        return np.array(
            [np.linspace(0.0, 1.0, npts), 1 - np.linspace(0.0, 1.0, npts)]
        ).T
    else:
        messages.error(
            "No default molar fractions in property_sanity_check when nc > 2"
        )
