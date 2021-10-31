#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box of 1000m in depth
# Initialized with gaz and
# imposed liquid Dirichlet BC at the bottom
# Homogeneous Neumann BC at both sides
# atm BC at the top
# gravity = 9.81

import ComPASS
import numpy as np
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.timestep_management import TimeStepManager
from ComPASS.mpi import master_print
from ComPASS.utils.grid import on_zmax

Lz = 1000.0
nz = 200
dz = Lz / nz
Lx = 4 * dz
Ly = 4 * dz
Ox, Oy, Oz = 0.0, 0.0, -1000.0
nx = 3
ny = 3
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12 * np.eye(3)  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Ptop = 1.0e5  # porous top Pressure
# Pbot = 1.E7                       # porous bottom Pressure
Patm = 1.0e5  # atm Pressure
Tporous = 300.0  # porous Temperature (used also to init the freeflow nodes)
CpRoche = 2.0e6

simulation = ComPASS.load_eos("diphasic_FreeFlowBC")
gravity = simulation.get_gravity()
simulation.set_atm_pressure(Patm)
simulation.set_atm_temperature(330)
simulation.set_atm_rain_flux(0.0)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

if not ComPASS.mpi.is_on_master_proc:
    grid = omega_reservoir = k_reservoir = cell_thermal_cond = None

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.bottom_boundary(grid),
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

fc = simulation.compute_face_centers()
simulation.set_freeflow_faces(on_zmax(grid)(fc))

# init gas
g_state = dict(p=Ptop, T=Tporous, S=[1, 0], C=[[1.0, 0.0], [0.0, 1.0]])
X0 = simulation.build_state(simulation.Context.gas, **g_state)
# top (freeflow) values
X_top = simulation.build_state(simulation.Context.gas_FF_no_liq_outflow, **g_state)
# bottom values (liquid)
Pbot = 5.0e5 + gravity * 1000 * Lz
l_state = dict(p=Pbot, T=Tporous, S=[0, 1], C=[[1.0, 0.0], [0.0, 1.0]])
X_bottom = simulation.build_state(simulation.Context.liquid, **l_state)

simulation.all_states().set(X0)
simulation.node_states().set(on_zmax(grid)(simulation.vertices()), X_top)
simulation.dirichlet_node_states().set(X_bottom)


master_print("set initial and BC")

# context = SimulationContext()
# context.abort_on_ksp_failure = False
# context.dump_system_on_ksp_failure = False
# context.abort_on_newton_failure = False

timestep = TimeStepManager(
    initial_timestep=1.0e-5 * year,
    minimum_timestep=1.0,
    maximum_timestep=30.0 * year,
    increase_factor=1.1,
    decrease_factor=0.6,
)

final_time = 20.0 * year
output_period = 0.01 * final_time

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-6, 8, lsolver)

current_time = simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=timestep,
    output_period=output_period,
    # nitermax=1,
)

print("time after the time loop", current_time / year, "years")
