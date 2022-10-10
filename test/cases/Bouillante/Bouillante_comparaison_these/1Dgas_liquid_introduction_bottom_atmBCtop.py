#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 200m in depth
# Initialized with gas and
# imposed liquid Dirichlet BC at the bottom
# Homogeneous Neumann BC at both sides
# atm BC at the top
# gravity = 0

import ComPASS
import numpy as np
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.grid import on_zmax

Lz = 200.0
nz = 50
dz = Lz / nz
Lx = 4 * dz
Ly = 4 * dz
Ox, Oy, Oz = 0.0, 0.0, -Lz
nx = 3
ny = 3
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Tporous = 300.0  # porous Temperature (used also to init the freeflow nodes)
Tatm = 330.0
CpRoche = 2.0e6

simulation = ComPASS.load_physics("diphasic")
simulation.set_liquid_capillary_pressure("Beaude2018")
gravity = simulation.get_gravity()
Patm = Ptop = simulation.get_atm_pressure()
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

# init gas
X0 = simulation.build_state(
    simulation.Context.gas,
    p=Ptop,
    T=Tporous,
)
# top (freeflow) values
Xtop = simulation.build_state(
    simulation.Context.gas_FF_no_liq_outflow,
    p=Ptop,
    T=Tporous,
)
# bottom values (liquid)
Pbot = Patm + gravity * 1000 * Lz
Xbot = simulation.build_state(
    simulation.Context.liquid,
    p=Pbot,
    T=Tporous,
)

simulation.all_states().set(X0)
simulation.node_states().set(is_ff, Xtop)
simulation.dirichlet_node_states().set(Xbot)


# context = SimulationContext()
# context.abort_on_ksp_failure = False
# context.dump_system_on_ksp_failure = False
# context.abort_on_newton_failure = False

timestep = TimeStepManager(
    initial_timestep=100,
    minimum_timestep=10.0,
    maximum_timestep=1.0 * year,
    increase_factor=1.1,
    decrease_factor=0.2,
)

final_time = 20.0 * year
output_period = 0.01 * final_time

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-6, 14, lsolver)

current_time = simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=timestep,
    output_period=output_period,
    # nitermax=1,
)
