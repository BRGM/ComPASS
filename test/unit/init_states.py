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
from ComPASS.timeloops import standard_loop


ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")

pres = 20.0 * MPa  # initial reservoir pressure
Tres = degC2K(
    70.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
omega_reservoir = 0.15  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K

Lx, Ly, Lz = 1000.0, 1000.0, 1000.0
Ox, Oy, Oz = 0, 0, 0
nx, ny, nz = 2, 2, 1

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

simulation.init(
    mesh=grid,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

# Create a specific state and set all degrees of freedom
initial_state = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(initial_state)

# You can iterate over states and print them
for state in simulation.cell_states():
    print(state)

# You can also iterate jointly over states and positions
# but this will not be very efficient
rho = 1E3
g = simulation.get_gravity()
for state, position in zip(simulation.all_states(), simulation.all_positions()):
    state.p = rho * g * position[2]

# Or rather use the direct accesor over the whole array
simulation.all_states().p[:] = rho * g * simulation.all_positions()[:, 2]

print(simulation.node_states().p)

# You can acces a specific state and change its values
S = simulation.all_states()[2]
S.p = 3
S.S[:]=0.5
print(simulation.all_states()[2])
