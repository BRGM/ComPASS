#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Detente en gaz 1D
# Avec une physique trop simplifiée (énergie interne = enthalpie)
# la temperature restait constante en temps.
# En ajoutant energie interne = enthalpie - Pression, on met un couplage entre la temperature et 
# la pression, ce qui fait varier la temperature dans ce cas test.


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
Lx=40.
Ly = 1.
Ox, Oy, Oz = 0.,     0.,    -3000.
nx = 2
ny = 1
Topz = Oz+Lz

omega_reservoir = 0.35            # reservoir porosity
k_reservoir = 1E-12 * np.eye(3)   # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 0.            # reservoir thermal conductivity
pall = 1.E6                       # initial Pressure
pbot_dir = 1.E5                   # bottom Pressure
Tall = 700.                       # initial Temperature
Tbot_dir = 700.                   # bottom Temperature
gravity = 0.
CpRoche = 0.


bot_flag = 4

ComPASS.load_eos('diphasic_FreeFlowBC')
ComPASS.set_gravity(gravity)
ComPASS.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

gas_context = ComPASS.Context.gas

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
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[vertices[:,2]<=Oz] = bot_flag

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

# master_print('Maillage distribue')

def set_Dirichlet_state(state):
    node_flags = ComPASS.nodeflags()
    # bottom
    state.context[node_flags==bot_flag] = gas_context
    state.p[node_flags==bot_flag] = pbot_dir
    state.T[node_flags==bot_flag] = Tbot_dir
    state.S[node_flags==bot_flag] = [1, 0]
    state.C[node_flags==bot_flag] = [[ 0.8, 0.2], [0, 1.]]

def set_states(state):
    state.context[:] = gas_context
    state.p[:] = pall
    state.T[:] = Tall
    state.S[:] = [1, 0]
    state.C[:] = [[ 0.8, 0.2], [0, 1.]]

def set_initial_bc_values():
    set_states(ComPASS.node_states())
    set_states(ComPASS.cell_states())
    set_Dirichlet_state(ComPASS.dirichlet_node_states())

# master_print('set initial and BC')
set_initial_bc_values()

# set linear solver properties
newton = Newton(1e-7, 15, LinearSolver(1e-6, 50))

context = SimulationContext()
context.abort_on_ksp_failure = False 
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False 

timestep = TimeStepManager(initial_timestep = 100.,
minimum_timestep = 1E-1, maximum_timestep = 50.*year,
increase_factor = 1.2, decrease_factor = 0.2,
)

final_time = 1000. * year
output_period = 0.01 * final_time


current_time = standard_loop(final_time = final_time,
context=context, newton=newton,
time_step_manager = timestep,
output_period = output_period, 
)

print('time after the time loop', current_time/year, 'years')
