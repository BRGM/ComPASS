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

p0 = 1. * bar              # initial reservoir pressure
T0 = degC2K( 20. )         # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
bottom_heat_flux = 0.08    # W/m2                                  
k_matrix = 1E-18           # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15          # column porosity
K_matrix = 2.              # bulk thermal conductivity in W/m/K

H = 3000.                  # column height
nx, ny, nz = 1, 1, 300     # discretization

rhof = 1E3                 # specific mass in kg/m^3
cpf = 4200                 # specific heat in J/kg/K
rhofcpf = rhof * cpf       # volumetric heat capacity
muf = 1E-3                 # fluid dynamic viscosity

ComPASS.load_eos('linear_water')
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.volumetric_heat_capacity = rhofcpf
fluid_properties.dynamic_viscosity = muf


ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (10., 10., H),
    origin = (-5, -5, -H),
)

def top_nodes():
    return ComPASS.global_vertices()[:, 2] >= 0

ComPASS.init(
    mesh = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    cell_thermal_conductivity = K_matrix,
    set_dirichlet_nodes = top_nodes,
)

def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = 1. 
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
               ComPASS.node_states(),
               ComPASS.cell_states()]:
    set_initial_states(states)

def set_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux
    face_centers = ComPASS.face_centers()   
    ComPASS.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann) 
set_boundary_heat_flux()

final_time = 1E4 * year
output_period = 1E3 * year
ComPASS.set_maximum_timestep(10 * year)
standard_loop(initial_timestep= 30 * day, final_time = final_time, output_period = output_period)
