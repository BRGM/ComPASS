#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
#import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
import numpy as np

rhow =1# 1E3
b = 1#4.2e+3          # specific heat in J/kg/K
rhocp =1# 2000*800   # volumetric heat capacity
muf = 1#3E-4
pleft, pright = 300,200#30 * MPa, 10 * MPa
Tleft, Tright = degC2K(60), degC2K(100)
omega_reservoir = 0.2 # reservoir porosity
k_reservoir =1 # 1E-12 # reservoir permeability in m^2
K_reservoir = 10 #2     # bulk thermal conductivity in W/m/K
a=8

Lx = 100.
Ly = 80
Lz = 1
nx = 100
ny = 80
nz = 1
onecomp = True
exact_sol = True

if onecomp:
    if exact_sol:
        ComPASS.load_eos('linear_water')
    else:
        ComPASS.load_eos('liquid_water')
else:
    ComPASS.load_eos('water_with_tracer')
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhow
fluid_properties.volumetric_heat_capacity = b
fluid_properties.dynamic_viscosity = muf
ComPASS.set_rock_volumetric_heat_capacity(rhocp)

#mu = 3E-4 # dynamic viscosity of pur water around 100°C (will change with temperature)
U = ((k_reservoir / muf) * (pleft - pright) / Lx)

print('Average Darcy velocity:', U * year, 'm/year')
print('                  i.e.: %.2f%%' % (100 * U * year/ Lx), 'of the simulation domain in one year.')
final_time =5* Lx /(U/omega_reservoir)
print('Final time is set to: %.2f years' % (final_time/year))
## Velocity (with omega is 2.4 m / hour, good final time is 8 to 10 hours  ! Initial time step is 7e-3 seconds

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (0, 0 , 0),)

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

ComPASS.init(
    mesh = grid,
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

standard_loop(
    final_time = final_time, output_period = final_time/50,
    time_step_manager = TimeStepManager(final_time/1e3, 1.),
    output_callbacks=(collect_node_temperature,),
)

if ComPASS.mpi.communicator().size==1:
    assert ComPASS.mpi.is_on_master_proc
    xy = ComPASS.compute_cell_centers()[:,0:2]
    XX = xy[:, 0].reshape(ny, nx)
    YY = xy[:, 1].reshape(ny, nx)
    import ComPASS.utils.mpl_backends as mpl_backends
    plt = mpl_backends.import_pyplot(False)
    if plt:
        plt.clf()
        for tT in cell_temperatures:
            t, T = tT
            #T = K2degC(T)
            fig = plt.figure(1)
            #plt.subplot(211)
            cs = plt.contourf(XX,YY,np.reshape(K2degC(T),[ny,nx]))
            fig.colorbar(cs)
            plt.axis('scaled')

            #plt.subplot(212)
            #plt.plot(x,c1)
            plt.title('t='+str(t))
            plt.xlabel('x in meters')
            plt.ylabel('y in meters')
            plt.draw()
            plt.pause(0.1)
            plt.savefig(ComPASS.to_output_directory('cell_temperatures_'+str(t)),format='png')
            plt.clf()

        #plt.contourf(XX,YY, np.reshape(cell_temperatures[18][1], [ny,nx]))
            #plt.xlabel('x in meters')
            #plt.ylabel('temperature in Celsius degrees')
