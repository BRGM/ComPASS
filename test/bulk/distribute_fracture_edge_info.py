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
from ComPASS.timeloops import standard_loop, TimeStepManager
import ComPASS.mpi as mpi

k = 1E-15           # dummy value
phi = 0.15          # dummy value
thermal_cond = 2.   # dummy value

nx, ny, nz = (4,) * 3

ComPASS.load_eos('linear_water')
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    origin = (-0.5, -0.5, -0.5),
)

# This could be set from an external property
# We flag all nodes that belong to fractures
# We could flag only a set of nodes out of a list of edges, for example
# We use bitwise or to code different edge family 
nb_edge_families = 3
def set_node_flags():
    xyz = ComPASS.global_vertices()
    x, y, z = [xyz[:, j] for j in range(3)]
    edge_families = [
        (x == 0) & ((np.abs(y) == 0.5) | (np.abs(z) == 0.5)),
        (y == 0) & ((np.abs(x) == 0.5) | (np.abs(z) == 0.5)),
        (z == 0) & ((np.abs(x) == 0.5) | (np.abs(y) == 0.5)),
    ]
    flags = ComPASS.global_nodeflags()
    flags[:] = 0
    for i, edge_nodes in enumerate(edge_families):
        print('Family', i, 'has', np.sum(edge_nodes), 'nodes')
        flags[edge_nodes] = np.bitwise_or(flags[edge_nodes], 2**i)

def select_fractures():
    centers = ComPASS.compute_global_face_centers()
    xc, yc, zc = [centers[:, j] for j in range(3)]
    where = (xc == 0) | (yc == 0) | (zc == 0)
    print(np.sum(where), 'fracture faces')
    return where

ComPASS.init(
    mesh = grid,
    cell_permeability = k,                        # dummy value
    cell_porosity = phi,                          # dummy value
    cell_thermal_conductivity = thermal_cond,     # dummy value
    fracture_permeability = k,                    # dummy value
    fracture_porosity = phi,                      # dummy value
    fracture_thermal_conductivity = thermal_cond, # dummy value
    fracture_faces = select_fractures,
    set_global_flags = set_node_flags,
)

assert ComPASS.all_fracture_edges_tagged()
# Reselect fracture edges on local procs
fracture_edges = ComPASS.find_fracture_edges(ComPASS.frac_face_id())
print(fracture_edges.shape[0], 'fracture edges on proc', mpi.proc_rank)
flags = ComPASS.nodeflags()[fracture_edges]
for fi in range(nb_edge_families):
    mask = np.bitwise_and(flags, 2**fi)
    family_edges = fracture_edges[np.nonzero(mask[:,0] & mask[:, 1])]
    print('Family', fi, 'has', family_edges.shape[0], 'edges',
          'with', np.unique(family_edges).shape[0], 'nodes',
          'on proc', mpi.proc_rank, 
          # ':\n', family_edges,
    )
