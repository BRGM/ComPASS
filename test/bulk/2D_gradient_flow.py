#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

import numpy as np

pleft, pright = 30 * MPa, 10 * MPa
Tleft, Tright = degC2K(60), degC2K(100)
omega_reservoir = 0.2 # reservoir porosity
k_reservoir = 1E-12 # reservoir permeability in m^2
K_reservoir = 2                   # bulk thermal conductivity in W/m/K
a=10

Lx = 100.
Ly=50
Lz=1
nx = 100
ny=50
nz=1

onecomp = True
exact_sol = True

mu = 3E-4 # dynamic viscosity of pur water around 100°C (will change with temperature)
U = ((k_reservoir / mu) * (pleft - pright) / Lx)
print('Average Darcy velocity:', U * year, 'm/year')
print('                  i.e.: %.2f%%' % (100 * U * year/ Lx), 'of the simulation domain in one year.')
final_time =2* Lx /(U/omega_reservoir)
print('Final time is set to: %.2f years' % (final_time/year))
## Velocity (with omega is 2.4 m / hour, good final time is 8 to 10 hours  ! Initial time step is 7e-3 seconds

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
)

on_the_left = lambda x: x <= grid.origin[0]
on_the_left_in_the_middle = lambda x,y: (x <= grid.origin[0] ) & (grid.extent[1]/2-a <y) & (y< grid.extent[1]/2+a)
on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]

def select_dirichlet_nodes():
    x = ComPASS.global_vertices()[:,0]
    return on_the_left(x) | on_the_right(x)

def set_boundary_conditions():
    def set_states(states, x,y):
        left = on_the_left(x)
        states.p[left] = pleft
        leftmid = on_the_left_in_the_middle(x, y)
        states.T[left] = Tright #Tleft
        states.T[leftmid] = Tleft
        right = on_the_right(x)
        states.p[right] = pright
        states.T[right] = Tright
        both = (left | right)
        states.context[both] = 1
        states.S[both] = 1
        if onecomp:
            states.C[both] = 1.
        else:
            states.C[left] = (1, 0) #(0., 1)
            states.C[leftmid] = (0., 1.)
            states.C[right] = (1, 0)
    verts = ComPASS.vertices()
    set_states(ComPASS.dirichlet_node_states(), verts[:,0], verts[:,1])

def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pright # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.T[:] = Tright
        states.S[:] = 1
        if onecomp:
            states.C[:] = 1.
        else:
            states.C[:] = (1, 0)
    verts = ComPASS.vertices()
    cellcenters = ComPASS.compute_cell_centers()
    set_states(ComPASS.node_states(),  verts[:,0])
    set_states(ComPASS.cell_states(), cellcenters[:,0])

# %%% Simulation %%%

ComPASS.set_output_directory_and_logfile(__file__)

if onecomp:
    if exact_sol:
        ComPASS.load_eos('linear_liquid_water')
    else:
        ComPASS.load_eos('liquid_water')
else:
    ComPASS.load_eos('water_with_tracer')

ComPASS.init(
    grid = grid,
    set_dirichlet_nodes = select_dirichlet_nodes,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = K_reservoir,
)

set_initial_values()
set_boundary_conditions()

cell_temperatures = []
def collect_node_temperature(iteration, t):
    if ComPASS.mpi.communicator().size>1:
        if ComPASS.mpi.is_on_master_proc:
            print('WARNING - No output in parallel mode')
        return
    print('Collecting temperature at iteration', iteration)
    print('                           and time', t/year, 'years')
    states = ComPASS.cell_states()
    cell_temperatures.append((t, np.copy(states.T)))

standard_loop(final_time = final_time, output_period = final_time/50, initial_timestep = final_time/1e8,
              output_callbacks=(collect_node_temperature,))

if ComPASS.mpi.communicator().size==1:
    assert ComPASS.mpi.is_on_master_proc
    xy = ComPASS.compute_cell_centers()[:,0:2]
    XX = xy[:, 0].reshape(ny, nx) 
    YY = xy[:, 1].reshape(ny, nx) 
    # with open(ComPASS.to_output_directory('cell_temperatures.csv'), 'w') as f:
    #     s = ';'.join(['%f' %(xi) for xi in x])
    #     print('"time (years)\\x";' + s, file=f)
    #     for tT in cell_temperatures:
    #         t, T = tT
    #         T = K2degC(T)
    #         print('%f;' %(t/year) + ';'.join(['%f' %(Ti) for Ti in T]), file=f)
    try:
        import matplotlib

        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('WARNING - matplotlib was not found - no graphics will be generated')
        plt = None
    if plt is not None:
        plt.clf()
        for tT in cell_temperatures:
            t, T = tT
            T = K2degC(T)
        plt.contourf(XX,YY, np.reshape(cell_temperatures[18][1], [ny,nx]))
        plt.xlabel('x in meters')
        plt.ylabel('temperature in Celsius degrees')
        plt.savefig(ComPASS.to_output_directory('cell_temperatures'))

#etude de convergence par calcul d'erreur
