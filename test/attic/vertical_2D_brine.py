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
from ComPASS.timestep_management import TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("brine")
final_time = 200 * year
output_period = 10.0 * year

T = degC2K(25.0)
H = 0.1 * km
ptop = 1 * MPa
kres = 1e-12  # permeability reservoir
phires = 0.15  # porosity
Kres = 2  # bulk cell thermal conductivity W/m/K
gravity = 9.81
nz = 21

simulation.set_gravity(gravity)
grid = ComPASS.Grid(shape=(nz, 1, nz), extent=(H, H / nz, H), origin=(0.0, 0.0, -H),)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.top_boundary(grid),
    cell_permeability=kres,
    cell_porosity=phires,
    cell_thermal_conductivity=Kres,
)

X0 = simulation.build_state(p=ptop, T=T)
simulation.dirichlet_node_states().set(X0)
states = simulation.all_states()
states.set(X0)
x = simulation.all_positions()[:, 0]
salt_column = x > 0.5 * H
states.C[salt_column, 0, 0] = 0.4
states.C[salt_column, 0, 1] = 0.6

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

timestep = TimeStepManager(
    initial_timestep=10 * day, increase_factor=2.0, decrease_factor=0.2,
)

simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=timestep,
    output_period=output_period,
)
