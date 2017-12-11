import os
import numpy as np

import ComPASS
mpi = ComPASS.mpi
create_directories = ComPASS.utils.create_directories

class Dumper:
    
    def __init__(self, output_directory=None):
        if output_directory is None:
            output_directory = ComPASS.to_output_directory('')
        self.output_directory = os.path.abspath(output_directory)
        self.mesh_directory = os.path.join(self.output_directory, 'mesh')
        self.states_directory = os.path.join(self.output_directory, 'states')
        self.simulation_running = False
        self.proc_label_pattern = 'proc_{:06d}' #IMPROVE: hard coded maximum nb of procs

    def to_output_directory(self, filename=''):
        return os.path.join(self.output_directory, filename)
        
    def to_mesh_directory(self, filename=''):
        return os.path.join(self.mesh_directory, filename)
        
    def to_states_directory(self, filename=''):
        return os.path.join(self.states_directory, filename)

    def proc_label(self, proc):
        return self.proc_label_pattern.format(proc)

    def mesh_filename(self, proc):
        return self.to_mesh_directory('mesh_%s.npz' % self.proc_label(proc))

    def states_filename(self, proc, tag=''):
        if tag:
            tag+= '_'
        return self.to_states_directory('state_%s%s.npz' % (tag, self.proc_label(proc)))

    def start_simulation(self):
        assert not self.simulation_running
        create_directories(self.to_output_directory())
        create_directories(self.to_mesh_directory())
        create_directories(self.to_states_directory())
        self.simulation_running = True
        self.dump_own_element_numbers()
        self.dump_mesh()

    @mpi.on_master_proc
    def dump_own_element_numbers(self):
        filename = self.to_output_directory('own_elements')
        communicator = mpi.communicator()
        with open(filename, 'w') as f:
            print('# Number of procs', file=f)
            print('nb_procs =', communicator.size, file=f)
            print('nb_own_cells =', ComPASS.nb_cells_own(), file=f)
            print('nb_own_nodes =', ComPASS.nb_nodes_own(), file=f)
            print('nb_own_faces =', ComPASS.nb_faces_own(), file=f)

    def dump_mesh(self):
        connectivity = ComPASS.get_connectivity()
        fracture_faces = ComPASS.frac_face_id()
        fracture_nodes = [np.array(connectivity.NodebyFace[fk]) - 1 for fk in fracture_faces]  # switch first node indexing from 1 to 0 
        fracturenodes_offsets = np.cumsum([len(a) for a in fracture_nodes])
        fracturenodes_values = np.hstack(fracture_nodes) if fracturenodes_offsets else np.array([])
        fracture_types = ComPASS.facetypes()[fracture_faces]
        np.savez(self.mesh_filename(mpi.proc_rank),
            vertices =  ComPASS.vertices().view(dtype=np.double).reshape((-1, 3)),
            cellnodes_offsets = connectivity.NodebyCell.offsets()[1:], # VTK does not use the first 0 offset
            cellnodes_values = connectivity.NodebyCell.contiguous_content() - 1, # switch first node indexing from 1 to 0 
            celltypes = ComPASS.celltypes(),
            fracturenodes_offsets = fracturenodes_offsets,
            fracturenodes_values = fracturenodes_values,
            fracture_types = fracture_types,
        )

    def dump_states(self, tag=''):
        assert self.simulation_running
        node_states = ComPASS.node_states()
        cell_states = ComPASS.cell_states()
        fracture_states = ComPASS.fracture_states()
        np.savez(self.states_filename(mpi.proc_rank, tag),
            node_pressure = node_states.p,
            node_temperature = node_states.T,
            cell_pressure = cell_states.p,
            cell_temperature = cell_states.T,
            fracture_pressure = fracture_states.p,
            fracture_temperature = fracture_states.T,
        )

