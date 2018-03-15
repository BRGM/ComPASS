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

p0 = 1. * bar # initial reservoir pressure
T0 = degC2K( 20. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
k_matrix = 1E-18               # reservoir permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15            # reservoir porosity
                                  
H = 3000.
nx, ny, nz = 1, 1, 300

def top_nodes():
    vertices = np.rec.array(ComPASS.global_vertices())
    on_top = vertices.z >= H
    return on_top

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (1., 1., H),
    origin = (-0.5, -0.5, -H),
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    set_dirichlet_nodes = top_nodes,
)

def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = [0, 1]
    states.C[:] = 1.
set_initial_states(ComPASS.dirichlet_node_states())
set_initial_states(ComPASS.node_states())
set_initial_states(ComPASS.cell_states())

def set_boundary_heat_flux():
    face_centers = np.rec.array(ComPASS.face_centers())
    bottom_faces = np.nonzero(face_centers.z <= -H)[0]
    bottom_faces = np.asarray(bottom_faces + 1, dtype=np.int32) # Fortran indexing starts at 1   
    Neumann = ComPASS.NeumannBC()
    Neumann.heat_flux = 0.08 # W/m2
    ComPASS.set_Neumann_faces(bottom_faces, Neumann) 
set_boundary_heat_flux()

final_time = 1E6 * year
output_period = 1E4 * year
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 30 * day, final_time = final_time, output_period = output_period)
