import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.GridTools as GT
MT = ComPASS.ComPASS.MeshTools
#import vtkwriters as vtkw

L = 1E3
n = 1
gridextent = (L,) * 3
gridshape = (n,) * 3
vertices, tets = GT.grid2tets(gridshape, gridextent)
tets = np.asarray(tets, dtype=MT.idtype())

mesh = MT.tetmesh(vertices, tets)

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

zmax = mesh.vertices[:, -1].max()
topnodes = np.nonzero(mesh.vertices[:, -1]==zmax)[0]

ComPASS.set_output_directory_and_logfile(__file__)

def select_fractures():
    where = np.zeros(mesh.nb_faces(), dtype=np.bool)
    centers = ComPASS.compute_global_face_centers()
    where = on_plane(centers)
    print('Nb of nodes', mesh.nb_vertices(), 'vs', np.count_nonzero(where), 'fracture faces')
    return where

def select_dirichlet_nodes():
    print('Selecting', topnodes.shape[0], 'top nodes.')
    on_top = np.zeros(mesh.nb_vertices(), dtype=np.bool)
    on_top[topnodes] = True
    return on_top

def set_boundary_conditions():
    dirichlet = ComPASS.dirichlet_node_states()
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

print('Gravity:', ComPASS.gravity())

ComPASS.init(
    mesh = mesh,
    fracture_faces = select_fractures,    
    set_dirichlet_nodes = select_dirichlet_nodes
)
set_boundary_conditions()
set_initial_values()


standard_loop(final_time = 1E3 * year, output_frequency = 1E2 * year)
