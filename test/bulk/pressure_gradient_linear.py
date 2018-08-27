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

rhof = 1E3               # specific mass in kg/m^3
cpf = 4200               # specific heat in J/kg/K
hf = rhof * cpf          # specific enthalpy
muf = 1E-3               # viscosity Pa.s
p_right = 1. * bar       # initial reservoir pressure
T_right = degC2K( 20. )  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1E-12         # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15        # column porosity
K_matrix = 2             # bulk cell thermal conductivity W/m/K
mass_flux = 1E-1

Lx, Ly, Lz = 100., 10., 10.   # column dimensions
nx, ny, nz = 100, 1, 1        # discretization

ComPASS.load_eos('linear_water')
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.specific_enthalpy = hf
fluid_properties.viscosity = muf

ComPASS.set_output_directory_and_logfile(__file__)
ComPASS.set_gravity(0) # no gravity

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (0, -0.5 * Ly, -0.5 * Lz),
)

def left_nodes():
    return ComPASS.global_vertices()[:, 0] <= 0

def right_nodes():
    return ComPASS.global_vertices()[:, 0] >= Lx

ComPASS.init(
    mesh = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    cell_thermal_conductivity = K_matrix,
    set_dirichlet_nodes = right_nodes,
)

def set_initial_states(states):
    states.context[:] = 1
    states.p[:] = p_right
    states.T[:] = T_right
    states.S[:] = 1.
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
               ComPASS.node_states(),
               ComPASS.cell_states()]:
    set_initial_states(states)

def set_boundary_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = mass_flux # one component
    Neumann.heat_flux = mass_flux * hf
    face_centers = np.rec.array(ComPASS.face_centers())   
    ComPASS.set_Neumann_faces(face_centers[:, 0] <= 0, Neumann) 
set_boundary_flux()

final_time = 1E4 * year
output_period = 0.1 * final_time
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 1E-5, final_time = final_time, output_period = output_period)

x = ComPASS.cell_centers()[:, 0]
amean = lambda a: 0.5 * (a[:-1] + a[1:])
gradx_cells = lambda a: (a[1:] - a[:-1])/(x[1:] - x[:-1])
assert np.abs(- k_matrix * (rhof/muf) * gradx_cells(states.p) - mass_flux).max() < 1E-12

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
        plt.clf()
        plt.plot((x[0], x[-1]), (mass_flux, mass_flux), 'r')
        plt.plot(amean(x), - k_matrix * (rhof/muf) * gradx_cells(states.p), 'xk')
        plt.ylim(0.9 * mass_flux, 1.1 * mass_flux)
        plt.xlabel('x in meters')
        plt.ylabel('mass flux')
        plt.savefig(ComPASS.to_output_directory('mass_flux'))

