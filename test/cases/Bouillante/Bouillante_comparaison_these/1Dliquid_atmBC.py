#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 1000m in depth
# Imposed Dirichlet BC at the bottom
# Homogeneous Neumann BC at both sides
# atm BC at the top
# gravity = 0

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.grid import on_zmax

Lz = 1000.0
nz = 200
nx = ny = 3
dz = Lz / nz
Lx = Ly = nx * dz
Ox, Oy, Oz = 0.0, 0.0, 0.0
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12 * np.eye(3)  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Pporous = 1.0e7  # porous Pressure (used also to init the freeflow nodes)
Tatm = 300.0
Tporous = 200.0  # porous Temperature (used also to init the freeflow nodes)
Ttop = 250  # init atm BC temperature between Tporous and Tatm
CpRoche = 2.0e6
gravity = 0.0

simulation = ComPASS.load_eos("diphasic")
simulation.set_gravity(gravity)
simulation.set_liquid_capillary_pressure("Beaude2018")
simulation.set_atm_pressure(Pporous)
simulation.set_atm_temperature(Tatm)
simulation.set_atm_rain_flux(0.0)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.bottom_boundary(grid),
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

fc = simulation.compute_face_centers()
simulation.set_freeflow_faces(on_zmax(grid)(fc))
is_ff = simulation.get_freeflow_nodes()  # array of bool of size n_nodes

X0 = simulation.build_state(simulation.Context.liquid, p=Pporous, T=Tporous, Cal=0.01)
X_bottom = simulation.build_state(
    simulation.Context.liquid, p=Pporous, T=Tporous, Cal=0.0
)
X_top = simulation.build_state(
    simulation.Context.diphasic_FF_liq_outflow, p=Pporous, T=Ttop, Cal=0.01
)

simulation.all_states().set(X0)
simulation.node_states().set(is_ff, X_top)
simulation.dirichlet_node_states().set(X_bottom)

tsmger = TimeStepManager(
    initial_timestep=0.2 * year,
    minimum_timestep=1,
    maximum_timestep=10.0 * year,
    increase_factor=1.3,
    decrease_factor=0.6,
)

final_time = 40.0 * year
output_period = 0.1 * final_time

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-7, 15, lsolver)

run_loop = lambda final_time, no_output=True: simulation.standard_loop(
    reset_iteration_counter=True,
    initial_time=0,
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    # output_period=output_period,
    no_output=no_output,
)

run_loop(final_time)
