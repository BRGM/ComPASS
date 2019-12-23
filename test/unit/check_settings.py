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
import ComPASS.io.mesh as io


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
nx, ny, nz = 4, 4, 3

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

def select_dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
    return (x <= Ox) | (x >= Ox + Lx)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

# Create a specific state and set all degrees of freedom
initial_state = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(initial_state)

dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)
# Set varying dirichlet conditions
x = simulation.vertices()[:, 0]
y = simulation.vertices()[:, 1]
dirichlet.p[:] = ((pres + 1E5) * (x - Ox) + pres * (Ox + Lx - x)) / Lx
dirichlet.T[:] = ((Tres + 100) * (y - Oy) + Tres * (Oy + Ly - y)) / Ly

io.write_mesh(simulation, "mesh_alone")
io.write_mesh(
    simulation,
    "simulation_mesh",
    pointdata={
        "dirichlet": simulation.dirichlet_nodes(),
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": simulation.temperature_dirichlet_values(),
    },
    celldata={
        "initial pressure": simulation.cell_states().p,
        "zcell": simulation.compute_cell_centers()[:, 2],
    },
)
