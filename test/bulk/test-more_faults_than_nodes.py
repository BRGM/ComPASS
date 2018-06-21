#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import MeshTools as MT
import MeshTools.GridTools as GT
#import vtkwriters as vtkw

L = 1E3
n = 1
gridextent = (L,) * 3
gridshape = (n,) * 3
vertices, tets = GT.grid2tets(gridshape, gridextent)
tets = np.asarray(tets, dtype=MT.idtype())

mesh = MT.TetMesh.make(vertices, tets)

# The plane x-y=0
def on_plane_xy(pts, epsilon):
    return np.abs(pts[:, 0] - pts[:,1])<epsilon
# The plane x-z=0
def on_plane_xz(pts, epsilon):
    return np.abs(pts[:, 0] - pts[:,2])<epsilon
# The plane y-z=0
def on_plane_yz(pts, epsilon):
    return np.abs(pts[:, 1] - pts[:,2])<epsilon
def on_plane(pts, epsilon=1E-10*L):
    return on_plane_xy(pts, epsilon) | on_plane_xz(pts, epsilon) | on_plane_yz(pts, epsilon)

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)
print('Gravity:', ComPASS.get_gravity())

def toponodes():
    vertices = ComPASS.global_vertices()
    z = vertices[:, -1]
    return z==z.max()

def select_fractures():
    where = np.zeros(mesh.nb_faces, dtype=np.bool)
    centers = ComPASS.compute_global_face_centers()
    where = on_plane(centers)
    print('Nb of nodes', mesh.nb_vertices, 'vs', np.count_nonzero(where), 'fracture faces')
    return where

def select_dirichlet_nodes():
    vertices = ComPASS.global_vertices()
    z = vertices[:, -1]
    on_top = toponodes()
    print('Selecting', np.sum(on_top), 'top nodes.')
    return on_top

def set_global_flags():
    nodeflags = ComPASS.global_nodeflags()
    nodeflags[:] = 0
    nodeflags[toponodes()] = 1

def set_boundary_conditions():
    dirichlet = ComPASS.dirichlet_node_states()
    topnodes = ComPASS.nodeflags() == 1 
    dirichlet.p[topnodes] = 1E5
    dirichlet.T[topnodes] = degC2K(30)
    dirichlet.context[:] = 2
    dirichlet.S[:] = [0, 1]
    dirichlet.C[:] = 1.

def set_initial_values():
    for state in [ComPASS.node_states(), ComPASS.fracture_states(), ComPASS.cell_states()]:
        state.context[:] = 2
        state.p[:] = 1E5
        state.T[:] = degC2K(30)
        state.S[:] = [0, 1]
        state.C[:] = 1.

ComPASS.init(
    mesh = mesh,
    fracture_faces = select_fractures,
    set_dirichlet_nodes = select_dirichlet_nodes,
    cell_permeability = 1e-15,
    cell_porosity = 0.1,
    cell_thermal_conductivity = 2,
    fracture_permeability = 1e-15,
    fracture_porosity = 0.5,
    fracture_thermal_conductivity = 2,
    set_global_flags = set_global_flags,
)
set_boundary_conditions()
set_initial_values()


standard_loop(final_time = 1E3 * year, output_period = 1E2 * year)
