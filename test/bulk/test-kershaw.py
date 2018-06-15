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

from ComPASS.FVCAFile import FVCAFile
from ComPASS.RawMesh import RawMesh

# relies on Kershaw mesh from https://www.latp.univ-mrs.fr/latp_numerique/?q=node/11
fcva = FVCAFile('dkershaw08.msh')
mesh = RawMesh.convert(fcva)
epsilon = 1e-6 # to locate faces

p0 = 1. * bar              # initial reservoir pressure
T0 = degC2K( 20. )         # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
bottom_heat_flux = 0.08    # W/m2                                  
k_matrix = 1E-18           # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15          # column porosity
K_matrix = 2               # bulk thermal conductivity in W/m/K

ComPASS.load_eos('water2ph')
ComPASS.set_gravity(0)
ComPASS.set_output_directory_and_logfile(__file__)

def top_nodes():
    return ComPASS.global_vertices()[:,2] + epsilon >= 1

ComPASS.init(
    mesh = mesh,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    set_dirichlet_nodes = top_nodes,
    cell_thermal_conductivity = K_matrix,
)

def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = [0, 1]
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
               ComPASS.node_states(),
               ComPASS.cell_states()]:
    set_initial_states(states)

def set_boundary_heat_flux():
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = bottom_heat_flux
    ComPASS.set_Neumann_faces(ComPASS.face_centers()[:,2] - epsilon<= 0, Neumann) 
set_boundary_heat_flux()

final_time = 1 * hour
output_period = 0.1 * final_time
ComPASS.set_maximum_timestep(output_period)
standard_loop(
    initial_timestep= 0.1*output_period,
    final_time = final_time,
    output_period = output_period
)
