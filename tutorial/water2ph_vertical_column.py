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
from ComPASS.timeloops import TimeStepManager
from ComPASS.utils.grid import on_zmin, on_zmax

p0 = 1.0 * bar  # initial reservoir pressure
T0 = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1e-17  # column permeability in m^2
phi_matrix = 0.15  # column porosity
K_matrix = 1.0  # bulk thermal conductivity in W/m/K

H = 800.0  # column height
nx, ny, nz = 1, 1, 100  # discretization
Lx = nx * H / nz
Ly = ny * H / nz

simulation = ComPASS.load_physics("water2ph")  # water component in gas or liquid phase
simulation.set_gravity(0)  # no gravity

ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, H),
    origin=(-Lx / 2, -Ly / 2, -H),
)


def select_top_bottom_nodes():
    vertices = simulation.global_vertices()
    return on_zmin(grid)(vertices) | on_zmax(grid)(vertices)


# some initializations and partition of the mesh (if run in parallel)
simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=select_top_bottom_nodes,
)


def set_initial_states(states):
    gasPhase = simulation.phase_index(simulation.Phase.gas)
    liquidPhase = simulation.phase_index(simulation.Phase.liquid)
    states.context[:] = simulation.Context.liquid
    states.p[:] = p0
    states.T[:] = T0
    states.S[:, liquidPhase] = 1.0
    states.S[:, gasPhase] = 0.0
    states.C[:] = 1.0


for states in [
    simulation.dirichlet_node_states(),
    simulation.node_states(),
    simulation.cell_states(),
]:
    set_initial_states(states)

final_time = 1e4 * year
output_period = 1e3 * year
tsmger = TimeStepManager(
    initial_timestep=100 * day,
    maximum_timestep=output_period,
    minimum_timestep=1 * day,  # execution stops if smaller dt
    increase_factor=2.5,  # if success : next dt = increase_factor * previous dt
)

simulation.standard_loop(
    final_time=final_time,
    time_step_manager=tsmger,
    output_period=output_period,
)
