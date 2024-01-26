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
from ComPASS.timeloops import standard_loop, TimeStepManager


p0 = 1.0 * bar  # initial and top reservoir pressure
T0 = degC2K(
    20.0
)  # initial and top reservoir temperature - convert Celsius degrees to Kelvin
Tbot = T0 + 30  # bottom temperature
k_matrix = 1e-18  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K

H = 3000.0  # column height
nx, ny, nz = 1, 1, 300  # discretization

# simulation = ComPASS.load_physics("water2ph")
simulation = ComPASS.load_physics("diphasic")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(10.0, 10.0, H),
    origin=(-5, -5, -H),
)


def top_bott_nodes():
    return simulation.bottom_boundary(grid)() | simulation.top_boundary(grid)()


bottom_flag = 2
top_flag = 3


def set_physical_flags():
    nodeflags = simulation.global_nodeflags()

    nodeflags[simulation.top_boundary(grid)()] = top_flag
    nodeflags[simulation.bottom_boundary(grid)()] = bottom_flag


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_global_flags=set_physical_flags,
    set_pressure_dirichlet_nodes=simulation.top_boundary(grid),
    set_temperature_dirichlet_nodes=top_bott_nodes,
)

node_flags = simulation.nodeflags()

X0 = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)
# at the bottom, only the T is Dirichlet, the pressure value is not used,
# it is set on purpose with a value far from the result P (3e7) to test
Xbot = simulation.build_state(simulation.Context.liquid, p=1e6, T=Tbot)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(node_flags == top_flag, X0)
simulation.dirichlet_node_states().set(node_flags == bottom_flag, Xbot)


final_time = 1e3 * year
output_period = 1e2 * year
standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=TimeStepManager(30 * day, 100 * year),
    output_period=output_period,
)

assert np.all(simulation.node_states().T[node_flags == bottom_flag] == Tbot)
