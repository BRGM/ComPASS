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

import ComPASS
import numpy as np
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.simulation_context import SimulationContext
from ComPASS.newton import Newton, default_direct_solver
from ComPASS.timestep_management import TimeStepManager
from ComPASS.mpi import master_print

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

liquid_context = simulation.Context.liquid
diphasic_with_liq_outflow = simulation.Context.diphasic_FF_liq_outflow

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

    def Dirichlet_node():
        vertices = np.rec.array(simulation.global_vertices())
        return vertices[:, 2] <= Oz

    def set_global_flags():
        # nodes
        vertices = np.rec.array(simulation.global_vertices())
        nodeflags = simulation.global_nodeflags()
        nodeflags[vertices[:, 2] >= Topz] = freeflow_flag
        nodeflags[vertices[:, 2] <= Oz] = bot_flag
        # freeflow faces, necessary to flag them
        facecenters = simulation.compute_global_face_centers()
        faceflags = simulation.global_faceflags()
        faceflags[facecenters[:, 2] >= Topz] = freeflow_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = (
        Dirichlet_node
    ) = omega_reservoir = k_reservoir = cell_thermal_cond = set_global_flags = None

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=Dirichlet_node,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
    set_global_flags=set_global_flags,
)


def set_Dirichlet_state(state):
    node_flags = simulation.nodeflags()
    # bottom
    state.context[node_flags == bot_flag] = liquid_context
    state.p[node_flags == bot_flag] = Pporous
    state.T[node_flags == bot_flag] = Tporous
    state.S[node_flags == bot_flag] = [0, 1]
    state.C[node_flags == bot_flag] = [[1.0, 0.0], [0.0, 1.0]]


def set_FreeFlow_state(state):
    node_flags = simulation.nodeflags()
    # top
    state.context[node_flags == freeflow_flag] = diphasic_with_liq_outflow
    state.p[node_flags == freeflow_flag] = Pporous
    state.T[node_flags == freeflow_flag] = 250  # between Tporous and Tatm
    state.S[node_flags == freeflow_flag] = [0, 1]
    state.C[node_flags == freeflow_flag] = [[1.0, 0.0], [0.01, 0.99]]
    state.FreeFlow_phase_flowrate[node_flags == freeflow_flag] = 0.0


def set_states(state, depths):
    state.context[:] = liquid_context
    state.p[:] = Pporous  # no gravity
    state.T[:] = Tporous
    state.S[:] = [0, 1]
    state.C[:] = [[1.0, 0.0], [0.01, 0.99]]


def set_variable_initial_bc_values():
    set_states(simulation.node_states(), Topz - simulation.vertices()[:, 2])
    set_states(simulation.cell_states(), Topz - simulation.compute_cell_centers()[:, 2])
    set_FreeFlow_state(simulation.node_states())
    set_Dirichlet_state(simulation.dirichlet_node_states())


master_print("set initial and BC")
set_variable_initial_bc_values()

context = SimulationContext()
context.abort_on_ksp_failure = False
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False

timestep = TimeStepManager(
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
lsolver = default_direct_solver(simulation)
newton = Newton(simulation, 1e-5, 8, lsolver)

current_time = standard_loop(
    simulation,
    newton=newton,
    final_time=final_time,
    context=context,
    time_step_manager=timestep,
    output_period=output_period,
)

print("time after the time loop", current_time / year, "years")
