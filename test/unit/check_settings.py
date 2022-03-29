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
from ComPASS.utils.various import tensor_coordinates

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

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)


def select_dirichlet_nodes():
    # x1 | y1 | z1
    # x2 | y2 | z2
    # .  | .  | .
    # .  | .  | .
    # xn | yn | zn
    xyz = simulation.global_vertices()  # an array of shape (n vertices, 3)
    x = xyz[:, 0]  # retrieve the first column, i.e. x
    return (x <= x.min()) | (x >= x.max())


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
dirichlet.p[:] = ((pres + 1e5) * (x - Ox) + pres * (Ox + Lx - x)) / Lx
dirichlet.T[:] = ((Tres + 100) * (y - Oy) + Tres * (Oy + Ly - y)) / Ly

# The following will export all meshes:
# - as a vtu file if the simulation is sequential (a single mesh)
# - as a pvtu file and a set of vtu files (in the vtu folder) otherwise
io.write_mesh(simulation, "mesh_alone")
petrophysics = simulation.petrophysics()
pointdata = {
    "dirichlet": simulation.dirichlet_nodes(),
    "dirichlet pressure": simulation.pressure_dirichlet_values(),
    "dirichlet temperature": simulation.temperature_dirichlet_values(),
}
celldata = {
    "initial pressure": simulation.cell_states().p,
    "zcell": simulation.compute_cell_centers()[:, 2],
    "phi": petrophysics.cell_porosity,
}
celldata.update(
    tensor_coordinates(petrophysics.cell_permeability, "k", diagonal_only=True)
)
celldata.update(
    tensor_coordinates(
        petrophysics.cell_permeability, "k"
    )  # write whole array as a tensor array
)
celldata.update(
    tensor_coordinates(petrophysics.cell_thermal_conductivity, "K", diagonal_only=True)
)
# Creates the file simulation_mesh(.vtu if the simulation is sequential,
# .pvtu otherwise) in the execution directory
io.write_mesh(simulation, "simulation_mesh", pointdata=pointdata, celldata=celldata)
