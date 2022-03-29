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

simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)
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
states.S[:, :] = 0.5
states.context[:] = simulation.Context.diphasic

assert np.allclose(
    simulation.total_phase_volume(simulation.Phase.gas), 0.5 * volume * phi
)
