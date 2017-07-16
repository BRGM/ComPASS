import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.GridTools as GT
MT = ComPASS.ComPASS.MeshTools
#import vtkwriters as vtkw

gridshape = (1, 1, 1)
gridextent = (1E3, 1E3, 1E3)
vertices, tets = GT.grid2tets(gridshape, gridextent)
tets = np.asarray(tets, dtype=MT.idtype())

mesh = MT.tetmesh(vertices, tets)

zmax = mesh.vertices[:, -1].max()
topnodes = np.nonzero(mesh.vertices[:, -1]==zmax)[0]

ComPASS.set_output_directory_and_logfile(__file__)

def select_dirichlet_nodes():
    print('Selecting', topnodes.shape[0], 'top nodes.')
    on_top = np.zeros(mesh.nb_vertices(), dtype=np.bool)
    on_top[topnodes] = True
    return on_top

def set_boundary_conditions():
    dirichlet.p[topo_nodes] = 1E5
    dirichlet.T[topo_nodes] = 30 + 273.15 
    dirichlet.context[:] = 2
    dirichlet.S[:] = [0, 1]
    dirichlet.C[:] = 1.

def set_initial_values():
    for state in [ComPASS.node_states(), ComPASS.fracture_states(), ComPASS.cell_states()]:
        state.context[:] = 2
        state.p[:] = 1E5
        state.T[:] = 30 + 273.15 
        state.S[:] = [0, 1]
        state.C[:] = 1.

print('Gravity:', ComPASS.gravity())

ComPASS.init(
    mesh = mesh,
    set_dirichlet_nodes = select_dirichlet_nodes
)
set_boundary_conditions()
set_initial_values()


#standard_loop(final_time = 30 * year, output_frequency = year)
