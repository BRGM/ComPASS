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

p_right = 1. * bar              # initial reservoir pressure
T_right = degC2K( 20. )         # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1E-12           # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15          # column porosity
mass_flux = 1E-1

Lx, Ly, Lz = 100., 10., 10.                  # column dimensions
nx, ny, nz = 100, 1, 1  # discretization

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)
ComPASS.set_gravity(0) # no gravity

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (0, -0.5 * Ly, -0.5 * Lz),
)

def left_nodes():
    vertices = np.rec.array(ComPASS.global_vertices())
    return vertices.x <= 0

def right_nodes():
    vertices = np.rec.array(ComPASS.global_vertices())
    return vertices.x >= Lx

ComPASS.init(
    grid = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    #set_pressure_dirichlet_nodes = right_nodes,
    #set_temperature_dirichlet_nodes = lambda: left_nodes() | right_nodes(),
    set_dirichlet_nodes = right_nodes,
)

def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p_right
    states.T[:] = T_right
    states.S[:] = [0, 1]
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
               ComPASS.node_states(),
               ComPASS.cell_states()]:
    set_initial_states(states)

def set_boundary_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = mass_flux # one component
    Neumann.heat_flux = mass_flux * ComPASS.liquid_molar_enthalpy(p_right, T_right)
    face_centers = np.rec.array(ComPASS.face_centers())   
    ComPASS.set_Neumann_faces(face_centers.x <= 0, Neumann) 
set_boundary_flux()

final_time = 1E4 * year
output_period = 0.1 * final_time
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 1 * year, final_time = final_time, output_period = output_period)

if ComPASS.mpi.communicator().size==1:
    assert ComPASS.mpi.is_on_master_proc
    states = ComPASS.cell_states()
    print(np.min(states.p) / bar, "bar <= pressure <=", np.max(states.p) / bar, "bar")
    print(K2degC(np.min(states.T)), "deg C <= temperature <=", K2degC(np.max(states.T)), "deg C")
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('WARNING - matplotlib was not found - no graphics will be generated')
        plt = None
    else:
        cell_centers = np.rec.array(ComPASS.cell_centers())
        x = cell_centers.x
        mu = ComPASS.liquid_dynamic_viscosity(states.p, states.T)
        rho = ComPASS.liquid_molar_density(states.p, states.T)
        plt.clf()
        plt.subplot(121)
        plt.plot(x, mu)
        plt.xlabel('x in meters')
        plt.ylabel('dynamic viscosity')
        plt.subplot(122)
        plt.plot(x, rho)
        plt.xlabel('x in meters')
        plt.ylabel('specific mass')
        plt.savefig(ComPASS.to_output_directory('cell_properties'))
        plt.clf()
        plt.subplot(111)
        amean = lambda a: 0.5 * (a[:-1] + a[1:])
        gradx_cells = lambda a: (a[1:] - a[:-1])/(x[1:] - x[:-1])
        plt.plot((x[0], x[-1]), (mass_flux, mass_flux), 'r')
        plt.plot(amean(x), - k_matrix * amean(rho/mu) * gradx_cells(states.p), 'xk')
        plt.ylim(0.9 * mass_flux, 1.1 * mass_flux)
        plt.xlabel('x in meters')
        plt.ylabel('mass flux')
        plt.savefig(ComPASS.to_output_directory('mass_flux'))
