#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box similar in depth with Bouillante
# Imposed molar and heat flux at the bottom boundary

import ComPASS
import numpy as np
import MeshTools as MT
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
Lx=Ly=2*dz
Ox, Oy, Oz = 0.,     0.,    -3000.
nx=ny=2
Topz = Oz+Lz

omega_reservoir = 0.35            # reservoir porosity
k_reservoir = 1E-12               # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.            # reservoir thermal conductivity
ptop = 0.1 * MPa                  # top Pressure
Ttop = 300.                       # top Temperature
pbot = 2.5E7
Tbot = 400.                       # bottom Temperature
Tbot_input = 550.                       # bottom Temperature input
qbin = 2.E-2 #329./Lx/Ly
CpRoche = 2.E6
gravity = 10.
liq_molar_fraction_qbin_air = 1.E-4   # fixme: should be 0
liq_molar_fraction_qbin = np.array([liq_molar_fraction_qbin_air,1. - liq_molar_fraction_qbin_air])

geotherm = (Tbot-Ttop)/Lz

bot_flag = 4
top_flag = 5

simulation = ComPASS.load_eos('diphasic')
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

simulation.info.activate_direct_solver = True

gas_context = simulation.Context.gas
liquid_context = simulation.Context.liquid
diphasic_context = simulation.Context.diphasic

#assert simulation.Phase.gas==1
#assert simulation.Phase.liquid==2
#assert simulation.Component.air==1
#assert simulation.Component.water==2

if ComPASS.mpi.is_on_master_proc:
    
    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
        origin = (Ox, Oy, Oz),
    )
    
    def Dirichlet_node():
        vertices = np.rec.array(simulation.global_vertices())
        return (vertices[:,2] >= Topz)
    
    def set_global_flags():
        vertices = np.rec.array(simulation.global_vertices())
        nodeflags = simulation.global_nodeflags()
        nodeflags[vertices[:,2]>=Topz] = top_flag
        nodeflags[vertices[:,2]<=Oz] = bot_flag
        
        face_centers = np.rec.array(simulation.compute_global_face_centers())
        faceflags = simulation.global_faceflags()
        faceflags[:] = 0
        faceflags[face_centers[:,2]<=Oz] = bot_flag

if not ComPASS.mpi.is_on_master_proc:
    grid = Dirichlet_node =  omega_reservoir = k_reservoir = cell_thermal_cond = set_global_flags = None

simulation.init(
    mesh = grid,
    set_dirichlet_nodes = Dirichlet_node,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = cell_thermal_cond,
    set_global_flags = set_global_flags,
)

#master_print('Maillage distribue')
sys.stdout.write('Maillage distribue' + '\n')
sys.stdout.flush()


def to_array(xyz):
    array = np.zeros((xyz.shape[0],3))
    for i, elt in enumerate(xyz):
        array[i] = np.fromiter(elt,dtype=np.float)
    return array

def lininterp(depths, top, gradient):
    assert(np.all(depths)>=0),"depths is not positif"
    return top + (gradient)*(depths)

def set_molar_fraction(Pg, T):
    Hur = .5
    Cwg = Hur*simulation.Psat(T)/Pg
    Cag = 1. - Cwg

    T1 = 293.
    T2 = 353.
    H1 = 6.e+9
    H2 = 10.e+9
    Ha = H1 + (H2-H1)*(T-T1)/(T2-T1)

    Cal = Cag*Pg/Ha
    Cwl = 1. - Cal
    return [[Cag, Cwg], [Cal, Cwl]]

def set_Dirichlet_state(state):
    node_flags = simulation.nodeflags()
    # top
    state.context[node_flags==top_flag] = gas_context
    state.p[node_flags==top_flag] = ptop
    state.T[node_flags==top_flag] = Ttop
    state.S[node_flags==top_flag] = [1, 0]
    state.C[node_flags==top_flag] = set_molar_fraction(ptop, Ttop)

def set_states(state, depths):
    state.context[:] = liquid_context
    state.p[:] = lininterp(depths, ptop, gravity*1000.)
    state.T[:] = lininterp(depths, Ttop, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = [[ 1., 0.], [1.E-5, 1.-1.E-5]]   # FIXME air=0

def set_variable_initial_bc_values():
    set_states(simulation.node_states(), Topz-simulation.vertices()[:,2])
    set_states(simulation.cell_states(), Topz-simulation.compute_cell_centers()[:,2])
    set_Dirichlet_state(simulation.dirichlet_node_states())

# in some part of the bottom, molar  and heat flux.
def set_variable_boundary_heat_flux():
    face_flags = simulation.faceflags()

    neumann_heat_faces = np.zeros(len(face_flags),dtype = bool)
    neumann_heat_faces[face_flags==bot_flag] = True

    # heat source
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = liq_molar_fraction_qbin[:]*qbin
    Neumann.heat_flux = qbin * simulation.liquid_molar_enthalpy(pbot,Tbot_input,liq_molar_fraction_qbin) 
    simulation.set_Neumann_faces(neumann_heat_faces, Neumann)


#master_print('set initial and BC')
sys.stdout.write('set initial and BC' + '\n')
set_variable_initial_bc_values()
#master_print('set Neumann BC')
sys.stdout.write('set Neumann BC' + '\n')
set_variable_boundary_heat_flux()
sys.stdout.flush()

# set linear solver properties
newton = Newton(simulation, 1e-5, 10, LinearSolver(1e-6, 50))

context = SimulationContext()
context.abort_on_ksp_failure = False 
context.dump_system_on_ksp_failure = False
context.abort_on_newton_failure = False 

timestep = TimeStepManager(initial_timestep = 1000.,
minimum_timestep = 1E-8, maximum_timestep = 10.*year,
increase_factor = 1.2, decrease_factor = 0.2,
)

final_time = 300. * year
output_period = 0.01 * final_time


current_time = standard_loop(simulation, final_time = final_time, 
#output_every=500,
context=context, newton=newton,
time_step_manager = timestep,
#nitermax=1,
#iteration_callbacks=[ma_fonction],
output_period = output_period, 
#specific_outputs=[1. * day], output_every=20,
)

print('time after the time loop', current_time/year, 'years')
