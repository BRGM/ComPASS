import sys
import numpy as np
from mpi4py import MPI

from .ComPASS import *
import ComPASS.utils.filenames
import ComPASS.runtime as runtime

def set_output_directory_and_logfile(case_name):
    runtime.output_directory, runtime.logfile = ComPASS.utils.filenames.output_directory_and_logfile(case_name)

class Grid:
    def __init__(self, shape, extent=None, origin=None):
        shape = tuple(shape)
        dim = len(shape)
        assert dim<=3
        if origin is None:
            origin = (0.,) * dim
        origin = tuple(origin)
        assert len(origin)==dim
        if extent is None:
            extent = (1.,) * dim
        extent = tuple(extent)
        assert len(extent)==dim
        if dim<3:
            shape+= (1,) * (3 - dim) # grid is always in 3D
            origin+= (0.,) * (3 - dim) # grid is always in 3D
            shape+= (1.,) * (3 - dim) # grid is always in 3D
        self.shape = shape
        self.extent = extent
        self.origin = origin

def on_master_proc(f):
    def call(*args, **kwargs):
        comm = MPI.COMM_WORLD
        if comm.rank==0:
            f(*args, **kwargs)
    return call

def synchronize():
    comm.Barrier() # wait for every process to synchronize

def init(
    meshfile=None, grid=None,
    wells = lambda: [],
    fracture_faces = lambda: None,
    cells_porosity = lambda: None,
    fracture_porosity = lambda: None,
    cells_permeability = lambda: None,
    fracture_permeability = lambda: None,
):
    assert meshfile is None or grid is None
    assert not(meshfile is None and grid is None)
    if meshfile:
        print('Loading mesh from file is not implemented here')
        # FIXME: This should be something like MPI.Abort()
        sys.exit(-1)
    elif grid:
        comm = MPI.COMM_WORLD
        ComPASS.init_warmup(runtime.logfile)
        if comm.rank==0:
            ComPASS.build_grid(shape = grid.shape, origin = grid.origin, extent = grid.extent)
    else:
        print('No mesh!')
        # FIXME: This should be something like MPI.Abort()
        sys.exit(-1)
    if comm.rank==0:
        well_list = list(wells())
        ComPASS.set_well_geometries(well_list)
        ComPASS.global_mesh_mesh_bounding_box()
        ComPASS.global_mesh_compute_all_connectivies()
        ComPASS.global_mesh_set_frac()
        ComPASS.global_mesh_node_of_frac()
        ComPASS.global_mesh_set_dir_BC()
        ComPASS.global_mesh_frac_by_node()
        ComPASS.global_mesh_make_post_read_set_poroperm()
        ComPASS.global_mesh_make_post_read_well_connectivity_and_ip()
        ComPASS.set_well_data(well_list)
        ComPASS.compute_well_indices()
    #comm.Barrier() # wait for every process to synchronize
    ComPASS.init_phase2(runtime.output_directory)
    comm.Barrier() # wait for every process to synchronize


def get_vertices():
   return np.array(ComPASS.get_vertices_buffer(), copy = False)

def get_id_faces():
   return np.array(ComPASS.get_id_faces_buffer(), copy = False)

def get_cell_permeability():
   return np.array(ComPASS.get_cell_permeability_buffer(), copy = False)

def get_fracture_permeability():
   return np.array(ComPASS.get_fracture_permeability_buffer(), copy = False)

def get_cell_porosity():
   return np.array(ComPASS.get_cell_porosity_buffer(), copy = False)

def get_fracture_porosity():
   return np.array(ComPASS.get_fracture_porosity_buffer(), copy = False)

def compute_face_centers():
    vertices = get_vertices()
    connectivity = get_connectivity()
    centers = np.array([
        vertices[
          np.array(face_nodes, copy=False) - 1 # fortran indexes start at 1
        ].mean(axis=0)
        for face_nodes in connectivity.NodebyFace
      ])
    return centers

def set_fractures(faces):
    idfaces = get_id_faces()
    idfaces[faces] = -2
    global_mesh_set_frac() # this will collect faces with flag -2 as fracture faces

@on_master_proc
def timestep_summary():
    ComPASS.summarize_timestep()

def output_visualization_files(iteration):
    ComPASS.output_visu(iteration, runtime.output_directory)

