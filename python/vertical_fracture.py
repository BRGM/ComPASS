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
qmass = 1E-1               # 
Tbottom = degC2K( 250. )   # temperature influx
k_matrix = 1E-15           # domain permeability in m^2
phi_matrix = 0.15          # domain porosity
k_fracture = 1E-12         # fracture permeability in m^2
phi_fracture = 0.3         # fracture porosity

H = 1000.                  # domain height
nH = 50                    # discretization
nx, ny, nz = 2*nH, 1, nH
Lx, Ly, Lz = 2*H, 0.1*H, H

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

# thermodynamic functions are only available once the eos is loaded
pbottom = ComPASS.get_gravity() * H * 900.
hbottom = ComPASS.liquid_molar_enthalpy(pbottom, Tbottom)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (-0.5*Lx, -0.5*Ly, -H)
)

def top_nodes():
    vertices = np.rec.array(ComPASS.global_vertices())
    return vertices.z >= 0

def select_fractures():
    centers = ComPASS.compute_global_face_centers()
    xc = centers[:, 0]
    zc = centers[:, 2]
    return (xc == 0) & (zc < -0.5 * H)

ComPASS.init(
    grid = grid,
    cell_permeability = k_matrix,
    cell_porosity = phi_matrix,
    fracture_faces = select_fractures,
    fracture_permeability = k_fracture,
    fracture_porosity = phi_fracture,
    set_dirichlet_nodes = top_nodes,
)

def set_initial_states(states):
    states.context[:] = 2
    states.p[:] = p0
    states.T[:] = T0
    states.S[:] = [0, 1]
    states.C[:] = 1.
for states in [ComPASS.dirichlet_node_states(),
               ComPASS.node_states(),
               ComPASS.cell_states(),
               ComPASS.fracture_states()]:
    set_initial_states(states)

def set_boundary_fluxes():
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = qmass
    Neumann.heat_flux = qmass * hbottom
    face_centers = np.rec.array(ComPASS.face_centers())   
    bottom_fracture_edges = ComPASS.find_fracture_edges(face_centers.z <= -H)
    ComPASS.set_Neumann_fracture_edges(bottom_fracture_edges, Neumann) 
set_boundary_fluxes()

final_time = 50 * year
output_period = 0.05 * final_time
ComPASS.set_maximum_timestep(output_period)
standard_loop(initial_timestep= 1 * hour, final_time = final_time, output_period = output_period)
