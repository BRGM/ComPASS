#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.properties.utils import constant_physical_property
from ComPASS.properties.densities import build_pure_phase_volumetric_mass_density
from ComPASS.properties.enthalpies import build_pure_phase_enthalpy


p0 = 1.0  # initial reservoir pressure
T0 = 1.0  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
phi = 0.15  # reservoir porosity
k = 1e-12  # reservoir permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K
rhof = 1.0
rhocpf = 0.8
rhocpr = 0.9

# A 10m cube (ie 1E3m volume) with center at (0,0,0)
Lx, Ly, Lz = 10, 10, 10
Ox, Oy, Oz = -5, -5, -5
nx, ny, nz = 10, 10, 10
volume = Lx * Ly * Lz

simulation = ComPASS.load_physics("linear_water")
simulation.set_gravity(0)

simulation.set_molar_enthalpy_functions(
    build_pure_phase_enthalpy(volumetric_heat_capacity=rhocpf),
)
simulation.set_molar_density_functions(
    build_pure_phase_volumetric_mass_density(
        specific_mass=rhof,
        compressibility=0.1,
        thermal_expansivity=0.1,
        reference_pressure=p0,
        reference_temperature=T0,
    ),
    # p0 and T0 are not scaled, out of range of the exp function
    # so not possible to check the densities values
    check_derivatives=False,
)
simulation.set_viscosity_functions(constant_physical_property(1.0))
simulation.set_rock_volumetric_heat_capacity(rhocpr)


grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

ComPASS.set_output_directory_and_logfile(__file__)

simulation.init(
    mesh=grid,
    cell_porosity=phi,
    cell_permeability=k,
    cell_thermal_conductivity=K,
)

# Init physical values
states = simulation.all_states()
states.p[:] = p0
states.T[:] = T0
states.C[:] = 1.0
states.S[:] = 1.0
states.context[:] = simulation.Context.single_context

# Initial state
mass, energy = simulation.total_accumulation()

assert np.allclose(mass, volume * phi * rhof)
assert np.allclose(energy, volume * (phi * rhocpf + (1 - phi) * rhocpr))
assert np.allclose(
    simulation.total_phase_volume(simulation.Phase.single_phase), volume * phi
)
