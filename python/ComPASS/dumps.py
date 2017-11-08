import os
import numpy as np

import ComPASS

@ComPASS.mpi.on_master_proc
def dump_own_element_numbers():
    filename = ComPASS.to_output_directory('own_elements')
    communicator = ComPASS.mpi.communicator()
    with open(filename, 'w') as f:
        print('# Number of procs', file=f)
        print('nb_procs =', communicator.size, file=f)
        print('nb_own_cells =', ComPASS.nb_cells_own(), file=f)
        print('nb_own_nodes =', ComPASS.nb_nodes_own(), file=f)
        print('nb_own_faces =', ComPASS.nb_faces_own(), file=f)

def dump_mesh(subdirectory=''):
    connectivity = ComPASS.get_connectivity()
    fracture_faces = ComPASS.frac_face_id()
    fracture_nodes = [np.array(connectivity.NodebyFace[fk]) - 1 for fk in fracture_faces]  # switch first node indexing from 1 to 0 
    fracturenodes_offsets = np.cumsum([len(a) for a in fracture_nodes])
    print('offset', fracturenodes_offsets)
    fracturenodes_values = np.hstack(fracture_nodes) if fracturenodes_offsets else np.array([])
    fracture_types = ComPASS.facetypes()[fracture_faces]
    directory = ComPASS.to_output_directory(subdirectory)
    assert os.path.isdir(directory)
    np.savez(os.path.join(directory, 'mesh_proc_%04d' % ComPASS.mpi.proc_rank),
        vertices =  ComPASS.vertices().view(dtype=np.double).reshape((-1, 3)),
        cellnodes_offsets = connectivity.NodebyCell.offsets()[1:], # VTK does not use the first 0 offset
        cellnodes_values = connectivity.NodebyCell.contiguous_content() - 1, # switch first node indexing from 1 to 0 
        celltypes = ComPASS.celltypes(),
        fracturenodes_offsets = fracturenodes_offsets,
        fracturenodes_values = fracturenodes_values,
        fracture_types = fracture_types,
    )

def dump_states(tag='', subdirectory=''):
    if tag:
        tag+= '_'
    directory = ComPASS.to_output_directory(subdirectory)
    assert os.path.isdir(directory)
    node_states = ComPASS.node_states()
    cell_states = ComPASS.cell_states()
    fracture_states = ComPASS.fracture_states()
    np.savez(os.path.join(directory, 'state_%sproc_%04d' % (tag, ComPASS.mpi.proc_rank)),
        node_pressure = node_states.p,
        node_temperature = node_states.T,
        cell_pressure = cell_states.p,
        cell_temperature = cell_states.T,
        fracture_pressure = fracture_states.p,
        fracture_temperature = fracture_states.T,
    )
