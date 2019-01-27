#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid
# top Dirichlet monophasic gas Hur .5

import ComPASS
import numpy as np
import MeshTools as MT
import MeshTools.CGALWrappers as CGAL
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.mpi


omega_reservoir = 0.35            # reservoir porosity  .1
k_reservoir = 1E-12               # reservoir permeability in m^2
cell_thermal_cond = 2.            # reservoir thermal conductivity
p0 = 0.1 * MPa                    # top Pressure
T0 = degC2K(20)                   # top Temperature
Tbot = degC2K(40)
CpRoche = 2.E6
gravity = 10.
geotherm = 0.02


Lx, Ly, Lz = 100., 100.,1.E3
Ox, Oy, Oz = 0., 0., 0.
nx, ny, nz = 3, 3, 20
Topz = Oz+Lz


gas_context = 1
liquid_context = 2
diphasic_context = 3


top_flag = 2
bottom_flag = 3


ComPASS.load_eos('diphasic')
ComPASS.set_gravity(gravity)
ComPASS.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

if ComPASS.mpi.is_on_master_proc:
    
    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
        origin = (Ox, Oy, Oz),
    )
    
    def Dirichlet_node():
        vertices = np.rec.array(ComPASS.global_vertices())
        result = np.zeros(len(vertices[:,0]),dtype = bool)
        result[vertices[:,2] >= Oz+Lz] = True
        result[vertices[:,2] <= Oz] = True
        return result
    
    def set_global_flags():
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[:] = 0
        nodeflags[(vertices[:,2] >= Oz+Lz)] = top_flag
        nodeflags[(vertices[:,2] <= Oz)] = bottom_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = Dirichlet_node = set_global_flags = None

ComPASS.init(
    mesh = grid,
    set_dirichlet_nodes = Dirichlet_node,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = cell_thermal_cond,
    set_global_flags = set_global_flags,
)

sys.stdout.write('Maillage distribue' + '\n')
sys.stdout.flush()


def to_array(xyz):
    array = np.zeros((xyz.shape[0],3))
    for i, elt in enumerate(xyz):
        array[i] = np.fromiter(elt,dtype=np.float)
    return array

def lininterp(depths, top, gradient):
    return top + (gradient)*(depths)

def molar_fraction_balance(Pg, T):
    Hur = .5
    Cwg = Hur*ComPASS.Psat(T)/Pg
    Cag = 1. - Cwg
    Ha = 1.E8
    Cal = Cag*Pg/Ha
    Cwl = 1. - Cal
    return [[Cag, Cwg], [Cal, Cwl]]

def set_Dirichlet_state(state): #top nodes
    # context 1:Ig; 2:Il; 3:Ig+Il
    node_flags = ComPASS.nodeflags()
    # top monophasic gas with water vapour (call fonction to have hur ?)
    state.context[node_flags==top_flag] = gas_context
    state.p[node_flags==top_flag] = p0
    state.T[node_flags==top_flag] = T0
    state.S[node_flags==top_flag] = [1, 0]
    state.C[node_flags==top_flag] = molar_fraction_balance(p0, T0)
    # bot monophasic liquid
    state.context[node_flags==bottom_flag] = liquid_context
    state.p[node_flags==bottom_flag] = p0 + 800.*gravity*(Lz-Oz)
    state.T[node_flags==bottom_flag] = Tbot
    state.S[node_flags==bottom_flag] = [0, 1] # gas, liquid
    state.C[node_flags==bottom_flag] = [[ 1., 0.], [1.E-5, 1.-1.E-5]] # air, water


def set_states(state, depths):
    # context 1:Ig; 2:Il; 3:Ig+Il
    state.context[:] = liquid_context
    state.p[:] = lininterp(depths, p0, 800*gravity)
    state.T[:] = lininterp(depths, T0, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = [[ 1., 0.], [1.E-5, 1.-1.E-5]]
    
def set_variable_initial_bc_values():
    set_states(ComPASS.node_states(), abs(ComPASS.vertices()[:,2]-Topz))
    set_states(ComPASS.cell_states(), abs(ComPASS.compute_cell_centers()[:,2]-Topz))
    set_Dirichlet_state(ComPASS.dirichlet_node_states())


sys.stdout.write('set initial and BC' + '\n')
set_variable_initial_bc_values()
sys.stdout.flush()


init_dt = .01 * hour
final_time = 1000. * year
output_period = 0.01 * final_time
ComPASS.set_maximum_timestep(.7*year)


def ma_fonction(it_timestep,time):
    if(it_timestep>300):
        ComPASS.mpi.abort()


standard_loop(initial_timestep = init_dt, final_time = final_time, output_every=20,
iteration_callbacks=[ma_fonction],
#output_period = output_period, specific_outputs=[1. * day], output_every=20,
)

