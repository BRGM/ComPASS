# -*- coding: utf-8 -*-
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
import thirdparties
import MeshTools as MT
import vtkwriters as vtkw
import gmsh_reader

ComPASS.load_eos('water2ph')

filename = 'case-tests/andra/andra_gallery.msh'

nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)

# filter out 2D elements
elements = [(elt, tags) for elt, tags in elements
            if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)]
tags = [tags for elt, tags in elements]
physical = np.array([tag[0] for tag in tags])
elements = [elt for elt, tags in elements]

mesh = MT.HybridMesh.create(nodes, elements)

def set_physical_flags():
    flags = ComPASS.global_cellflags()
    flags[:] = physical

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = mesh,
    set_global_flags = set_physical_flags
)
#set_initial_values()
#set_boundary_conditions()
#
#standard_loop(final_time = 30 * year, output_frequency = year)

# dump mesh info
connectivity = ComPASS.get_connectivity()
fracture_faces = ComPASS.frac_face_id()
fracture_nodes = [np.array(connectivity.NodebyFace[fk]) - 1 for fk in fracture_faces]  # switch first node indexing from 1 to 0 
fracturenodes_offsets = np.cumsum([len(a) for a in fracture_nodes])
print('offset', fracturenodes_offsets)
fracturenodes_values = np.hstack(fracture_nodes) if fracturenodes_offsets else np.array([])
fracture_types = ComPASS.facetypes()[fracture_faces]
np.savez(ComPASS.to_output_directory('mesh_proc_%04d' % ComPASS.mpi.proc_rank),
    vertices =  ComPASS.vertices().view(dtype=np.double).reshape((-1, 3)),
    cellnodes_offsets = connectivity.NodebyCell.offsets()[1:], # VTK does not use the first 0 offset
    cellnodes_values = connectivity.NodebyCell.contiguous_content() - 1, # switch first node indexing from 1 to 0 
    celltypes = ComPASS.celltypes(),
    fracturenodes_offsets = fracturenodes_offsets,
    fracturenodes_values = fracturenodes_values,
    fracture_types = fracture_types,
)

# dump states
node_states = ComPASS.node_states()
cell_states = ComPASS.cell_states()
fracture_states = ComPASS.fracture_states()
np.savez(ComPASS.to_output_directory('state_proc_%04d' % ComPASS.mpi.proc_rank),
    node_pressure = node_states.p,
    node_temperature = node_states.T,
    cell_pressure = cell_states.p,
    cell_temperature = cell_states.T,
    fracture_pressure = fracture_states.p,
    fracture_temperature = fracture_states.T,
)
