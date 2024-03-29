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
simulation = ComPASS.load_physics("water2ph")

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

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

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
print(f"Follows the cell states after initialization:")
for state in simulation.cell_states():
    print(state)

# You can also iterate jointly over states and positions
# but this will not be very efficient
rho = 1e3
g = simulation.get_gravity()
for state, position in zip(simulation.all_states(), simulation.all_positions()):
    state.p = rho * g * position[2]

# Or rather use the direct accesor over the whole array
simulation.all_states().p[:] = rho * g * simulation.all_positions()[:, 2]

print(f"Follows the node pressure states after modifications:")
print(simulation.node_states().p)

# You can acces a specific state and change its values
S = simulation.all_states()[2]
S.p = 3
S.S[:] = 0.5
print(f"Follows the specific state which has been modified:")
print(simulation.all_states()[2])

# You can specify a condition to restrain set() on a specific group
diphasic_state = simulation.build_state(simulation.Context.diphasic, p=pres, Sg=0.3)
vertices = simulation.vertices()
top_nodes = vertices[:, -1] > grid.origin[-1] + grid.extent[-1] / 2
simulation.node_states().set(top_nodes, diphasic_state)
print(f"Follows the node states after modification of top half (diphasic state):")
for z, ns in zip(vertices[:, -1], simulation.node_states()):
    print(f"with z coord {z}, state is {ns}")
