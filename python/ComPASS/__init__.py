import sys
import atexit
import importlib

# We must load mpi wrapper first so that MPI is  initialized before calling PETSC_Initialize
import ComPASS.mpi as mpi
import ComPASS.runtime as runtime
import ComPASS.utils.filenames

import numpy as np

# FIXME: This is to import thirdparties
import os
sys.path.append(os.path.abspath('..'))
import thirdparties
import MeshTools as MT

initialized = False

# CHECKME: There might be a more elegant way to do this
kernel = None

def load_eos(eosname):
    global kernel
    assert kernel is None
    kernel = importlib.import_module('ComPASS.eos.%s' % eosname)
    # CHECKME: we replace the behavior: from import ComPASS.eos.eosname import *
    #          there might be a more elegant way to do this
    gdict = globals()
    kdict = vars(kernel)
    for key in kdict:
        if key not in gdict:
            gdict[key] = kdict[key]

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

# FIXME: grid is kept for backward compatibility, should be deprecated
#        then mesh should not default to None and be explicitely provided
def init(
    mesh=None, grid=None,
    wells = lambda: [],
    fracture_faces = lambda: None,
    cells_porosity = lambda: None,
    faces_porosity = lambda: None,
    cells_permeability = lambda: None,
    faces_permeability = lambda: None,
    fractures_permeability = lambda: None,
    set_dirichlet_nodes = lambda: None,
    set_global_flags = None
):
    # FUTURE: This could be managed through a context manager ? 
    global initialized
    assert not initialized
    # FIXME: grid is kept for backward compatibility, should be deprecated
    assert not (grid is None and mesh is None)
    if grid is not None:
        assert mesh is None
        mesh = grid
    if type(mesh) is str:
        if not os.path.exists(mesh):
            print('Mesh file (%s) not found!' % mesh)
        print('Loading mesh from file is not implemented here')
        # FIXME: This should be something like MPI.Abort()
        sys.exit(-1) 
    elif type(mesh) is Grid:
        kernel.init_warmup(runtime.logfile)
        kernel.global_mesh_set_cartesian_mesh()
        if mpi.is_on_master_proc:
            kernel.build_grid(shape = grid.shape, origin = grid.origin, extent = grid.extent)
    else:
#    elif type(mesh) in [MeshTools.TetMesh, MeshTools.HexMesh]:
        kernel.init_warmup(runtime.logfile)
        if mpi.is_on_master_proc:
            print("!!!")
            print("!!! VTK output will not work with this mesh type!")
            print("!!!")
            vertices = MT.as_coordinate_array(mesh.vertices)
            cells_nodes, cells_faces, faces_nodes = mesh.COC_data()
            kernel.create_mesh(vertices,
                               *cells_nodes,
                               *cells_faces,
                               *faces_nodes
                               )
            # distribute cell types to be able to reconstruct local meshes
            celltypes = ComPASS.global_celltypes()
            celltypes[:] = mesh.vtk_ids().astype(np.int8, copy=False)
#    else:
#        print('Mesh type not understood!')
#        # FIXME: This should be something like MPI.Abort()
#        sys.exit(-1)
    if mpi.is_on_master_proc and set_global_flags is not None:
        assert callable(set_global_flags)
        set_global_flags()
    if mpi.is_on_master_proc:
        well_list = list(wells())
        kernel.set_well_geometries(well_list)
        kernel.global_mesh_mesh_bounding_box()
        kernel.global_mesh_compute_all_connectivies()
        fractures = fracture_faces()
        if fractures is not None:
            set_fractures(fractures)
        kernel.global_mesh_node_of_frac()
        #kernel.global_mesh_set_dir_BC()
        kernel.global_mesh_allocate_id_nodes()
        # Node information is reset first
        info = np.rec.array(global_node_info(), copy=False)
        for a in [info.pressure, info.temperature]:
            a[:] = ord('i')
        dirichlet = set_dirichlet_nodes()
        if dirichlet is not None:
            for a in [info.pressure, info.temperature]:
                a[dirichlet] = ord('d')
        kernel.global_mesh_count_dirichlet_nodes()
        kernel.global_mesh_frac_by_node()
        # The following line is necessary to allocate arrays in the fortran code
        kernel.global_mesh_make_post_read_set_poroperm()
        cellperm = cells_permeability()
        if cellperm is not None:
            kernel.get_cell_permeability()[:] = np.ascontiguousarray( cellperm )
        faceperm = faces_permeability()
        fracperm = fractures_permeability()
        if fractures is not None:
            if faceperm is not None:
                assert fracperm is None
                # the following assert is annoying when we just want to broadcast a values (typically a scalar value)
                # anyway assignement through the numpy.ndarray interface will fail
                # assert kernel.get_face_permeability().shape==faceperm.shape
                kernel.get_face_permeability()[:] = np.ascontiguousarray( faceperm )
            elif fracperm is not None:
                assert faceperm is None
                #the following assert is annoying when we just want to broadcast a values (typically a scalar value) 
                # anyway assignement through the numpy.ndarray interface will fail
                #assert fracperm.shape==tuple(np.count(fractures))
                kernel.get_face_permeability()[fractures] = np.ascontiguousarray( fracperm )
        kernel.global_mesh_make_post_read_well_connectivity_and_ip()
        kernel.set_well_data(well_list)
        kernel.compute_well_indices()
    kernel.init_phase2(runtime.output_directory)
    mpi.synchronize() # wait for every process to synchronize
    # FUTURE: This could be managed through a context manager ? 
    initialized = True
    atexit.register(finalize)

def get_id_faces():
   return np.array(kernel.get_id_faces_buffer(), copy = False)

def get_cell_permeability():
   return np.array(kernel.get_cell_permeability_buffer(), copy = False)

def get_face_permeability():
   return np.array(kernel.get_face_permeability_buffer(), copy = False)

def get_cell_porosity():
   return np.array(kernel.get_cell_porosity_buffer(), copy = False)

def get_face_porosity():
   return np.array(kernel.get_face_porosity_buffer(), copy = False)

def _compute_centers(points, elements):
    return np.array([
        points[
          np.array(element, copy=False) - 1 # fortran indexes start at 1
        ].mean(axis=0)
        for element in elements
      ])

def compute_global_cell_centers():
    return _compute_centers(
        global_vertices().view(dtype=np.double).reshape((-1, 3)),
        get_global_connectivity().NodebyCell
    )

def compute_global_face_centers():
    return _compute_centers(
        global_vertices().view(dtype=np.double).reshape((-1, 3)),
        get_global_connectivity().NodebyFace
    )

def compute_cell_centers():
    return _compute_centers(
        vertices().view(dtype=np.double).reshape((-1, 3)),
        get_connectivity().NodebyCell
    )

def compute_face_centers():
    return _compute_centers(
        vertices().view(dtype=np.double).reshape((-1, 3)),
        get_connectivity().NodebyFace
    )

#def compute_face_normals():
#    vertices = global_vertices().view(dtype=np.double).reshape((-1, 3))
#    connectivity = get_connectivity()
#    face_nodes = [np.array(nodes, copy=False) - 1 for nodes in connectivity.NodebyFace] # fortran indexes start at 1
#    normals = np.array([
#        np.cross(
#            vertices[nodes[1]]-vertices[nodes[0]],
#            vertices[nodes[2]]-vertices[nodes[0]])
#        for nodes in face_nodes
#    ])
#    # normalize
#    norms = np.linalg.norm(normals, axis=1)
#    norms.shape = (-1, 1)
#    normals /= norms
#    return normals

def get_boundary_faces():
    connectivity = get_connectivity()
    return np.array([
        np.array(face_cells, copy=False).shape[0]==1
        for face_cells in connectivity.CellbyFace
      ])

def get_boundary_vertices():
    connectivity = get_connectivity()
    boundary_faces = get_boundary_faces()
    vertices_id = np.unique(
            np.hstack([
                np.array(face_nodes, copy=False)
                for boundary, face_nodes in zip(boundary_faces, connectivity.NodebyFace) if boundary       
            ])
        )
    vertices_id -= 1 # Fortran index...
    nbnodes = len(connectivity.CellbyNode)
    res = np.zeros(nbnodes, dtype=np.bool)
    res[vertices_id] = True
    return res

def set_fractures(faces):
    idfaces = get_id_faces()
    idfaces[faces] = -2
    global_mesh_set_frac() # this will collect faces with flag -2 as fracture faces

@mpi.on_master_proc
def timestep_summary():
    kernel.summarize_timestep()

def output_visualization_files(iteration):
    kernel.output_visu(iteration, runtime.output_directory)

def to_output_directory(filename):
    assert os.path.isdir(runtime.output_directory)
    return os.path.abspath(os.path.join(runtime.output_directory, filename))


