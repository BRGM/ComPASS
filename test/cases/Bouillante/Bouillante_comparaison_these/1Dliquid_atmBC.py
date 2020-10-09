#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 4000m in depth
# Imposed Dirichlet BC at the bottom
# Homogeneous Neumann BC at both sides
# atm BC at the top
# gravity = 0

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.simulation_context import SimulationContext
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.grid import on_zmin, on_zmax

Lz = 1000.0
nz = 200
dz = Lz / nz
Lx = 3 * dz
Ly = 3 * dz
Ox, Oy, Oz = 0.0, 0.0, 0.0
nx = 4
ny = 4
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12 * np.eye(3)  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Pporous = 1.0e7  # porous Pressure (used also to init the freeflow nodes)
Tporous = 200.0  # porous Temperature (used also to init the freeflow nodes)
CpRoche = 2.0e6
gravity = 0.0

bot_flag = 4
freeflow_flag = 30  # do not modify this number

simulation = ComPASS.load_eos("diphasic_FreeFlowBC")
simulation.set_gravity(gravity)
simulation.set_atm_pressure(Pporous)
simulation.set_atm_temperature(300.0)
simulation.set_atm_rain_flux(0.0)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.bottom_boundary(grid),
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

fc = simulation.compute_face_centers()
simulation.set_freeflow_faces(on_zmax(grid)(fc))

state = dict(p=Pporous, T=Tporous, S=[0, 1], C=[[1.0, 0.0], [0.01, 0.99]])
X0 = simulation.build_state(simulation.Context.liquid, **state)
X_bottom = simulation.build_state(X0)
X_bottom.C[1][:] = [0.0, 1.0]
stateTop = dict(p=Pporous, T=Tporous, S=[0.5, 0.5], C=[[1.0, 0.0], [0.01, 0.99]])
X_top = simulation.build_state(simulation.Context.diphasic_FF_liq_outflow, **stateTop)
X_top.T = 250  # between Tporous and Tatm

simulation.all_states().set(X0)
simulation.node_states().set(on_zmax(grid)(simulation.vertices()), X_top)
simulation.dirichlet_node_states().set(X_bottom)

tsmger = TimeStepManager(
    initial_timestep=100.0,
    minimum_timestep=1e-3,
    maximum_timestep=10.0 * year,
    increase_factor=1.2,
    decrease_factor=0.2,
)

final_time = 40.0 * year
output_period = 0.1 * final_time

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

current_time = simulation.standard_loop(
    newton=newton,
    final_time=final_time,
    time_step_manager=tsmger,
    output_period=output_period,
)

print("time after the time loop", current_time / year, "years")
