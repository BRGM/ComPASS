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
import ComPASS.GridTools as GT
MT = ComPASS.ComPASS.MeshTools

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
    set_dirichlet_nodes = select_dirichlet_nodes
)
set_boundary_conditions()
set_initial_values()


standard_loop(final_time = 1E3 * year, output_period = 1E2 * year)
