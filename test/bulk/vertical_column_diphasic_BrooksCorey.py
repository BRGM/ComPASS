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
from ComPASS.timestep_management import TimeStepManager

p0 = 1.0 * bar
T0 = degC2K(20.0)
k_matrix = 1e-12  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K
gravity = 9.81

H = 500.0  # column height
nx, ny, nz = 1, 1, 100  # discretization
ds = H / nz

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("diphasicCO2")
simulation.set_gravity(gravity)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(ds, ds, H),
    origin=(-0.5 * ds, -0.5 * ds, -H),
)


def select_global_rocktype():
    cellrocktype = simulation.global_cell_rocktypes()
    cellrocktype[:] = 1


simulation.init(
    mesh=grid,
    set_global_rocktype=select_global_rocktype,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=simulation.top_boundary(grid),
)


# Set petrophyics functions after initialization
simulation.set_extendedBrooksCorey_pc_SPE11a()
from data.kr import kr_functions

simulation.set_kr_functions(kr_functions)

# set the correct capillary pressure **before** using build_state
# p is the **liquid pressure**, takes the pc into account
X0 = simulation.build_state(simulation.Context.liquid, p=p0, T=T0, rocktype=1)
# p is the **gas pressure**, takes the pc into account
Xgas = simulation.build_state(simulation.Context.gas, p=p0, T=T0)
# p is the **reference pressure**, ie the gas pressure
Xdiph = simulation.build_state(
    simulation.Context.diphasic, p=p0, T=T0, Sg=0.5, rocktype=1
)


simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 20, lsolver)

final_time = 1000 * year
current_time = simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=TimeStepManager(
        1 * day, increase_factor=2.0, decrease_factor=0.1
    ),
    # output_period = 100 * year,
)

run_loop = lambda initial_time: simulation.standard_loop(
    initial_time=initial_time,
    final_time=initial_time + 10 * year,
    newton=newton,
    time_step_manager=TimeStepManager(
        1 * day, increase_factor=2.0, decrease_factor=0.1
    ),
    output_period=year,
)


simulation.dirichlet_node_states().set(Xdiph)
current_time = run_loop(current_time)

simulation.dirichlet_node_states().set(Xgas)
current_time = run_loop(current_time)

simulation.postprocess()
