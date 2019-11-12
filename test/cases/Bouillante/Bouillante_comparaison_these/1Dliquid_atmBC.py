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
from ComPASS.newton import Newton, LinearSolver
from ComPASS.timestep_management import TimeStepManager
from ComPASS.mpi import master_print

Lz=4000.
nz=200
dz=Lz/nz
Lx=dz
Ly=dz
Ox, Oy, Oz = 0.,     0.,    -3000.
nx = 1
ny = 1
Topz = Oz+Lz

omega_reservoir = 0.35            # reservoir porosity
k_reservoir = 1E-12 * np.eye(3)   # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.            # reservoir thermal conductivity : no thermal diffusion
Pporous = 1.E7                    # porous Pressure (used also to init the freeflow nodes)
Tporous = 200.                    # porous Temperature (used also to init the freeflow nodes)
CpRoche = 2.E6
gravity = 0.0

bot_flag = 4
freeflow_flag = 30 # do not modify this number

ComPASS.load_eos('diphasic_FreeFlowBC')
ComPASS.set_gravity(gravity)
ComPASS.set_atm_pressure(Pporous)
ComPASS.set_atm_temperature(300.)
ComPASS.set_atm_rain_flux(0.)
ComPASS.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.activate_direct_solver = True

liquid_context = ComPASS.Context.liquid
diphasic_with_liq_outflow = ComPASS.Context.diphasic_FF_liq_outflow

if ComPASS.mpi.is_on_master_proc:
    
    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
        origin = (Ox, Oy, Oz),
    )
    
    def Dirichlet_node():
        vertices = np.rec.array(ComPASS.global_vertices())
        return (vertices[:,2] <= Oz)
    
    def set_global_flags():
        # nodes
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[vertices[:,2]>=Topz] = freeflow_flag
        nodeflags[vertices[:,2]<=Oz] = bot_flag
        # freeflow faces, necessary to flag them
        facecenters = ComPASS.compute_global_face_centers()
        faceflags = ComPASS.global_faceflags()
        faceflags[facecenters[:,2]>=Topz] = freeflow_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = Dirichlet_node =  omega_reservoir = k_reservoir = cell_thermal_cond = set_global_flags = None

ComPASS.init(
    mesh = grid,
    set_dirichlet_nodes = Dirichlet_node,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = cell_thermal_cond,
    set_global_flags = set_global_flags,
)

def set_Dirichlet_state(state):
    node_flags = ComPASS.nodeflags()
    # bottom
    state.context[node_flags==bot_flag] = liquid_context
    state.p[node_flags==bot_flag] = Pporous
    state.T[node_flags==bot_flag] = Tporous 
    state.S[node_flags==bot_flag] = [0, 1]
    state.C[node_flags==bot_flag] = [[ 1., 0.], [0.01, .99]]  

def set_FreeFlow_state(state):
    node_flags = ComPASS.nodeflags()
    # top
    state.context[node_flags==freeflow_flag] = diphasic_with_liq_outflow  
    state.p[node_flags==freeflow_flag] = Pporous
    state.T[node_flags==freeflow_flag] = 250 # between Tporous and Tatm
    state.S[node_flags==freeflow_flag] = [0, 1]
    state.C[node_flags==freeflow_flag] = [[ 1., 0.], [0.01, .99]] 
    state.FreeFlow_phase_flowrate[node_flags==freeflow_flag] = 0.

def set_states(state, depths):
    state.context[:] = liquid_context
    state.p[:] = Pporous # no gravity
    state.T[:] = Tporous
    state.S[:] = [0, 1]
    state.C[:] = [[ 1., 0.], [0.01, .99]] 

def set_variable_initial_bc_values():
    set_states(ComPASS.node_states(), Topz-ComPASS.vertices()[:,2])
    set_states(ComPASS.cell_states(), Topz-ComPASS.compute_cell_centers()[:,2])
    set_FreeFlow_state(ComPASS.node_states())
    set_Dirichlet_state(ComPASS.dirichlet_node_states())



master_print('set initial and BC')
set_variable_initial_bc_values()

# set linear solver properties
newton = Newton(1e-7, 15, LinearSolver(1e-6, 50))

context = SimulationContext()
context.abort_on_ksp_failure = False 
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False 

timestep = TimeStepManager(initial_timestep = 100.,
minimum_timestep = 1E-3, maximum_timestep = 10.*year,
increase_factor = 1.2, decrease_factor = 0.2,
)

final_time = 40. * year
output_period = 0.1 * final_time

current_time = standard_loop(final_time = final_time, 
context=context, newton=newton,
time_step_manager = timestep,
output_period = output_period, 
)

print('time after the time loop', current_time/year, 'years')