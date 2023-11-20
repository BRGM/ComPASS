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
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

p0 = 1.0 * bar
T0 = degC2K(20.0)
k_matrix = 1e-12  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K
gravity = 9.81

# define the geometry of the mesh
H = 300.0  # column height
nx, ny, nz = 5, 5, 30  # discretization, nb of cells in each direction
ds = nx * H / nz

# all outputs will be written in directory named output-"namefile"
ComPASS.set_output_directory_and_logfile(__file__)
# "diphasic" = possibility to have 2 phases (liquid and gas)
# and 2 components (air and water)
simulation = ComPASS.load_physics("diphasic")
simulation.set_gravity(gravity)  # by default gravity is 0

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(ds, ds, H),
    origin=(-0.5 * ds, -0.5 * ds, -H),
)

simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=simulation.top_boundary(grid),
)

# build state at equilibrium with only liquid phase, at given pressure and temperature
X0 = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)
# set T0 values and Dirichlet values at X0 state
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

# set a cell heat source close to source_pos
source_pos = np.array([0.0, 0.0, -H / 2.0])
distance_from_source = np.linalg.norm(
    simulation.compute_cell_centers() - source_pos, axis=1
)
index_closest_cell = np.argmin(distance_from_source)
thermal_sources = simulation.cellthermalsource()
cell_geom = np.array(grid.extent) / np.array(grid.shape)
cell_volume = cell_geom[0] * cell_geom[1] * cell_geom[2]
# init thermal source value
thermal_sources[index_closest_cell] = cell_volume * 0.08  # m3 * W/m3

# can impose a heat flux at nodes or whatever object (cell, node, fracture)
# using instead nodethermalsource, all_thermal_sources
# together with corresponding objects positions:
# simulation.vertices(), simulation.all_positions()

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 20, lsolver)

final_time = 100 * year
simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    initial_timestep=1 * day,
    output_period=5 * year,
)

# prepare for paraview visu, convert temperature in celsius
simulation.postprocess(convert_temperature=True)
