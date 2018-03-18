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
        create_directories(self.to_output_directory())
        create_directories(self.to_mesh_directory())
        create_directories(self.to_states_directory())

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
        self.simulation_running = True
        self.dump_own_element_numbers()
        self.dump_mesh()

    @mpi.on_master_proc
    def dump_own_element_numbers(self):
        filename = self.to_output_directory('own_elements')
        communicator = mpi.communicator()
        np_linewidth_backup = np.get_printoptions()['linewidth']
        np.set_printoptions(linewidth=np.inf) # we want all outputs below on a single line
        with open(filename, 'w') as f:
            print('# Number of procs', file=f)
            print('nb_procs =', communicator.size, file=f)
            print('nb_own_cells =', ComPASS.nb_cells_own(), file=f)
            print('nb_own_nodes =', ComPASS.nb_nodes_own(), file=f)
            print('nb_own_faces =', ComPASS.nb_faces_own(), file=f)
        np.set_printoptions(linewidth=np_linewidth_backup)

    def dump_mesh(self):
        connectivity = ComPASS.get_connectivity()
        fracture_nodes_coc = ComPASS.get_nodes_by_fractures()
        fracture_nodes = [np.array(nodes) - 1 for nodes in fracture_nodes_coc]  # switch first node indexing from 1 to 0 
        fracturenodes_offsets = np.cumsum([len(a) for a in fracture_nodes])
        fracturenodes_values = np.hstack(fracture_nodes) if len(fracture_nodes)>0 else np.array([])
        fracture_faces = ComPASS.frac_face_id()
        fracture_types = ComPASS.facetypes()[fracture_faces - 1] # switch first node indexing from 1 to 0 
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
        dumped_states = {
            'node_pressure': node_states.p,
            'node_temperature': node_states.T,
            'cell_pressure': cell_states.p,
            'cell_temperature': cell_states.T,
            'fracture_pressure': fracture_states.p,
            'fracture_temperature': fracture_states.T,        
        }
        for phase in range(ComPASS.number_of_phases()):
            dumped_states['cell_saturation_phase_%d'%(phase+1)] = cell_states.S[:, phase] 
            dumped_states['node_saturation_phase_%d'%(phase+1)] = node_states.S[:, phase] 
            dumped_states['fracture_saturation_phase_%d'%(phase+1)] = fracture_states.S[:, phase] 
        for phase in range(ComPASS.number_of_phases()):
            for comp in range(ComPASS.number_of_components()):
                dumped_states['cell_comp%d_in_phase_%d'%(comp+1, phase+1)] = cell_states.C[:, phase, comp] 
                dumped_states['node_comp%d_in_phase_%d'%(comp+1, phase+1)] = node_states.C[:, phase, comp] 
                dumped_states['fracture_comp%d_in_phase_%d'%(comp+1, phase+1)] = fracture_states.C[:, phase, comp] 
        np.savez(self.states_filename(mpi.proc_rank, tag), **dumped_states)

def dump_mesh():
    dumper = Dumper()
    dumper.dump_own_element_numbers()
    dumper.dump_mesh()

def dump_states(tag=''):
    Dumper().dump_states(tag)

