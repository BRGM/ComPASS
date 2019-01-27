#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, initialized with gaz at the bottom and liquid at the top
# needs lot of Newton iterations to start the convergence

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
CpRoche = 2.E6
gravity = 10.
bottom_heat_flux = 0.1
geotherm = 0.02 #bottom_heat_flux / cell_thermal_cond # gradient Temperature/m


Lx, Ly, Lz = 100., 100., 1.E3
Ox, Oy, Oz = 0., 0., 0.
nx, ny, nz = 8, 8, 20

gas_flag = 1
liquid_flag = 2

gas_context = 1
liquid_context = 2
diphasic_context = 3


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
    
    def set_global_flags():
        boundary = Lz/3.
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[:] = 0
        nodeflags[vertices[:,2]-Oz >= boundary] = liquid_flag
        nodeflags[vertices[:,2]-Oz < boundary] = gas_flag
        
        cell_centers = np.rec.array(ComPASS.compute_global_cell_centers())
        cellflags = ComPASS.global_cellflags()
        cellflags[:] = 0
        cellflags[cell_centers[:,2]-Oz >= boundary] = liquid_flag
        cellflags[cell_centers[:,2]-Oz < boundary] = gas_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = set_global_flags = None

ComPASS.init(
    mesh = grid,
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


def set_states(state, flag, depths):
    state.context[flag==gas_flag] = gas_context
    state.p[flag==gas_flag] = lininterp(depths[flag==gas_flag], p0, 800.*gravity)
    state.T[flag==gas_flag] = lininterp(depths[flag==gas_flag], T0, geotherm)
    state.S[flag==gas_flag] = [1, 0]
    state.C[flag==gas_flag] = [[ .9, .1], [0., 1.]]
    
    state.context[flag==liquid_flag] = liquid_context
    state.p[flag==liquid_flag] = lininterp(depths[flag==liquid_flag], p0, 800.*gravity)
    state.T[flag==liquid_flag] = lininterp(depths[flag==liquid_flag], T0, geotherm)
    state.S[flag==liquid_flag] = [0, 1]
    state.C[flag==liquid_flag] = [[ 1., 0.], [1.E-5, 1.-1.E-5]]

def set_variable_initial_bc_values():
    set_states(ComPASS.node_states(), ComPASS.nodeflags(), Oz+Lz-ComPASS.vertices()[:,2])
    set_states(ComPASS.cell_states(), ComPASS.cellflags(), Oz+Lz-ComPASS.compute_cell_centers()[:,2])


sys.stdout.write('set initial and BC' + '\n')
set_variable_initial_bc_values()
sys.stdout.flush()


init_dt = .01 # * hour
final_time = 100. * year
output_period = 0.01 * final_time
ComPASS.set_maximum_timestep(.7*year)


def ma_fonction(it_timestep,time):
    if(it_timestep>1):
        ComPASS.mpi.abort()


standard_loop(initial_timestep = init_dt, final_time = final_time, output_every=20,
#iteration_callbacks=[ma_fonction],
#output_period = output_period, specific_outputs=[1. * day], output_every=20,
)

