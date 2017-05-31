import sys
import numpy as np

# We must load mpi4py first so that MPI is  initialized before calling PETSC_Initialize
from mpi4py import MPI

from .ComPASS import *
import ComPASS.utils.filenames
import ComPASS.runtime as runtime

proc_rank = MPI.COMM_WORLD.rank
is_on_master_proc = proc_rank==0

def on_master_proc(f):
    def call(*args, **kwargs):
        if is_on_master_proc:
            f(*args, **kwargs)
    return call

def synchronize():
    MPI.COMM_WORLD.Barrier() # wait for every process to synchronize

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

def init(
    meshfile=None, grid=None,
    wells = lambda: [],
    fracture_faces = lambda: None,
    cells_porosity = lambda: None,
    faces_porosity = lambda: None,
    cells_permeability = lambda: None,
    faces_permeability = lambda: None,
    fractures_permeability = lambda: None,
):
    assert meshfile is None or grid is None
    if meshfile:
        print('Loading mesh from file is not implemented here')
        # FIXME: This should be something like MPI.Abort()
        sys.exit(-1)
    else: 
        ComPASS.init_warmup(runtime.logfile)
        if is_on_master_proc:
            assert grid is not None
            ComPASS.build_grid(shape = grid.shape, origin = grid.origin, extent = grid.extent)
    if is_on_master_proc:
        well_list = list(wells())
        ComPASS.set_well_geometries(well_list)
        ComPASS.global_mesh_mesh_bounding_box()
        ComPASS.global_mesh_compute_all_connectivies()
        fractures = fracture_faces()
        if fractures is not None:
            set_fractures(fractures)
        else:
            ComPASS.global_mesh_set_frac()
        ComPASS.global_mesh_node_of_frac()
        ComPASS.global_mesh_set_dir_BC()
        ComPASS.global_mesh_frac_by_node()
        # The following line is necessary to allocate arrays in the fortran code
        ComPASS.global_mesh_make_post_read_set_poroperm()
        cellperm = cells_permeability()
        if cellperm is not None:
            ComPASS.get_cell_permeability()[:] = np.ascontiguousarray( cellperm )
        faceperm = faces_permeability()
        fracperm = fractures_permeability()
        if fractures is not None:
            if faceperm is not None:
                assert fracperm is None
                # the following assert is annoying when we just want to broadcast a values (typically a scalar value)
                # anyway assignement through the numpy.ndarray interface will fail
                # assert ComPASS.get_face_permeability().shape==faceperm.shape
                ComPASS.get_face_permeability()[:] = np.ascontiguousarray( faceperm )
            elif fracperm is not None:
                assert faceperm is None
                #the following assert is annoying when we just want to broadcast a values (typically a scalar value) 
                # anyway assignement through the numpy.ndarray interface will fail
                #assert fracperm.shape==tuple(np.count(fractures))
                ComPASS.get_face_permeability()[fractures] = np.ascontiguousarray( fracperm )
        ComPASS.global_mesh_make_post_read_well_connectivity_and_ip()
        ComPASS.set_well_data(well_list)
        ComPASS.compute_well_indices()
    ComPASS.init_phase2(runtime.output_directory)
    synchronize() # wait for every process to synchronize


def get_vertices():
   return np.array(ComPASS.get_vertices_buffer(), copy = False)

def get_id_faces():
   return np.array(ComPASS.get_id_faces_buffer(), copy = False)

def get_cell_permeability():
   return np.array(ComPASS.get_cell_permeability_buffer(), copy = False)

def get_face_permeability():
   return np.array(ComPASS.get_face_permeability_buffer(), copy = False)

def get_cell_porosity():
   return np.array(ComPASS.get_cell_porosity_buffer(), copy = False)

def get_face_porosity():
   return np.array(ComPASS.get_face_porosity_buffer(), copy = False)

def compute_cell_centers():
    vertices = get_vertices()
    connectivity = get_connectivity()
    centers = np.array([
        vertices[
          np.array(cell_nodes, copy=False) - 1 # fortran indexes start at 1
        ].mean(axis=0)
        for cell_nodes in connectivity.NodebyCell
      ])
    return centers

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

def compute_face_normals():
    vertices = get_vertices()
    connectivity = get_connectivity()
    # fortran indexes start at 1
    face_nodes = [np.array(nodes, copy=False) - 1 for nodes in connectivity.NodebyFace]
    normals = np.array([
        np.cross(
            vertices[nodes[1]]-vertices[nodes[0]],
            vertices[nodes[2]]-vertices[nodes[0]])
        for nodes in face_nodes
    ])
    # normalize
    norms = np.linalg.norm(normals, axis=1)
    norms.shape = (-1, 1)
    normals /= norms
    return normals

def set_fractures(faces):
    idfaces = get_id_faces()
    idfaces[faces] = -2
    global_mesh_set_frac() # this will collect faces with flag -2 as fracture faces

@on_master_proc
def timestep_summary():
    ComPASS.summarize_timestep()

def output_visualization_files(iteration):
    ComPASS.output_visu(iteration, runtime.output_directory)

