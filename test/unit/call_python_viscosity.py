import numpy as np
from ComPASS.properties.viscosities import liquid_diphasic_viscosities
from ComPASS.utils.units import *


p = 1.0 * bar
T = 280.0
C = np.array([0.0, 1.0])

# Test the implementation outside of ComPASS (without
# loading an eos, then without the simulation object)
from ComPASS.properties.physical_properties import PhaseStateStruct

# define the type of Xalpha (needs the number of components)
phase_state_type = PhaseStateStruct(2)
X = phase_state_type.Xalpha(p, T, C)
dfdX = phase_state_type.empty_Xalpha()
val_without_derivatives = liquid_diphasic_viscosities.without_derivatives(X)
val_with_derivatives = liquid_diphasic_viscosities.with_derivatives(X, dfdX)
assert np.allclose(val_without_derivatives, val_with_derivatives)
print(
    "without loading an eos the viscosity is ",
    val_with_derivatives,
    " and the derivatives are ",
    dfdX,
)
# Create the object to use the vectorize function
from ComPASS.properties.physical_properties import CompiledPhaseProperty

compiled_liquid_diphasic_viscosities = CompiledPhaseProperty(
    phase_state_type,
    liquid_diphasic_viscosities.with_derivatives,
    liquid_diphasic_viscosities.without_derivatives,
)

npts = 5
Xarray = np.empty(npts, dtype=phase_state_type.Xt)
Xarray["pressure"] = np.logspace(np.log10(bar), np.log10(50 * MPa), npts)
Xarray["temperature"] = np.linspace(273.15, 573.15, npts)
Xarray["molar_fractions"] = np.array(
    [np.linspace(0.0, 1.0, npts), 1 - np.linspace(0.0, 1.0, npts)]
).T
print(
    "Using the vectorize property ",
    compiled_liquid_diphasic_viscosities.vect_without_derivatives(Xarray),
)

# Test the utilities of ComPASS
import ComPASS

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("diphasic")
# X = simulation.Xalpha(p, T, C)
# print("using the simulation facilities",
#     liquid_diphasic_viscosities.without_derivatives(X))

# Test the python wrapper
print("Using the simulation facilities ", simulation.liquid_dynamic_viscosity(p, T, C))

# Test the cpp wrapper
print("Test with the cpp function ", simulation.cpp_liquid_dynamic_viscosity(p, T, C))
