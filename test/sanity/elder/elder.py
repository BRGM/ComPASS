#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.load_eos('water2ph')
ComPASS.set_gravity(9.81)
ComPASS.set_rock_volumetric_heat_capacity(0)

p0 = 1. * bar # initial reservoir pressure
T0 = degC2K( 12. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
T1 = degC2K( 20. )
k_matrix = 1.21E-10             # reservoir permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.1                # reservoir porosity
K_matrix = 2                    # bulk thermal conductivity in W/m/K
                                
Lx, Ly, Lz = 600., 10., 150.
#Ox, Oy, Oz = -1500., -1000., -1600.
nx, ny, nz = 100, 1, 100

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    #origin = (Ox, Oy, Oz),
)

def bottom_dirichlet_nodes(vertices):
    on_bottom = vertices[:,2] <= 0.
    on_bottom&= (vertices[:,0]>=0.25*Lx) & (vertices[:,0]<=0.75*Lx)
    return on_bottom

def pressure_base_node():
    vertices = ComPASS.global_vertices()
    zmax = vertices[:,2].max()
    xmin = vertices[:,0].min()
    xmax = vertices[:,0].max()
    on_corner = (vertices[:,2] == zmax) & ((vertices[:,0] == xmin) | (vertices[:,0] == xmax))
    return on_corner

def dirichlet_nodes():
    vertices = np.rec.array(ComPASS.global_vertices())
    zmax = vertices[:,2].max()
    on_top = vertices[:,2] >= zmax
    return on_top | bottom_dirichlet_nodes(vertices)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    cell_thermal_conductivity = K_matrix,
    set_temperature_dirichlet_nodes = dirichlet_nodes,
    set_pressure_dirichlet_nodes = pressure_base_node,
)

def set_initial_states(states, z):
    g = ComPASS.get_gravity()
    rho = ComPASS.liquid_molar_density(p0, T0)
    states.context[:] = 2
    states.p[:] = p0 + rho * g * (Lz - z)
    states.T[:] = T0
    states.S[:] = [0, 1]
    states.C[:] = 1.
set_initial_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:,2])
set_initial_states(ComPASS.node_states(), ComPASS.vertices()[:,2])
set_initial_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,2])

states = ComPASS.dirichlet_node_states()
where = bottom_dirichlet_nodes(ComPASS.vertices())
states.T[where] = T1

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as tri

def my_graph(n, t):
    vertices = ComPASS.vertices()
    states = ComPASS.node_states()
    where = vertices[:,1]==0
    x = vertices[:,0][where]
    z = vertices[:,2][where]
    T = K2degC(states.T[where])
    plt.clf()
    plt.tricontourf(x, z, T, np.linspace(K2degC(T0), K2degC(T1), 10))
    plt.title('mon graph au temps %.2f ans' % (t / year))
    plt.colorbar()
    plt.savefig('graph_%05d.png' %n)

final_time = 5E1 * year
output_period = 1E0 * year
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 30 * day, final_time = final_time, output_period = output_period,
                                 output_callbacks=[my_graph,])
