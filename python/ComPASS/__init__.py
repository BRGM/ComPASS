#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import atexit
import importlib
import copy

# We must load mpi wrapper first so that MPI is  initialized before calling PETSC_Initialize
import ComPASS.mpi as mpi
import ComPASS.runtime as runtime
import ComPASS.utils.filenames
import ComPASS.dumps

import numpy as np

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

# This is temporary but will be generalized in the future
# here Properties will just be used as a namespace
class Properties:
    pass

# FIXME: grid is kept for backward compatibility, should be deprecated
#        then mesh should not default to None and be explicitely provided
def init(
    mesh=None, grid=None,
    wells = lambda: [],
    fracture_faces = lambda: None,
    set_dirichlet_nodes = lambda: None,
    set_pressure_dirichlet_nodes = lambda: None,
    set_temperature_dirichlet_nodes = lambda: None,
    set_global_flags = None,
    set_global_rocktype = None,
    **kwargs
):
    # here properties will just be used as a namespace
    properties = {}
    def make_property_accessor(value):
        return lambda: value
    for location in ['cell', 'fracture']:
        for property in ['porosity', 'permeability']:
            name = '%s_%s' % (location, property)
            try:
                f = kwargs[name]
                if callable(f):
                    properties[name] = f
                else:
                    properties[name] = make_property_accessor(f)
            except KeyError:
                properties[name] = lambda: None
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
        print('Loading mesh from file is desactivated here.')
        # FIXME: This should be something like MPI.Abort()
        sys.exit(-1)
    elif type(mesh) is Grid:
        kernel.init_warmup(runtime.logfile)
        kernel.global_mesh_set_cartesian_mesh()
        if mpi.is_on_master_proc:
            kernel.build_grid(shape = grid.shape, origin = grid.origin, extent = grid.extent)
            celltypes = ComPASS.global_celltypes()
            celltypes[:] = 11 # VTK_VOXEL
            facetypes = ComPASS.global_facetypes()
            #facetypes[:] = 8 # VTK_PIXEL
            facetypes[:] = 9 # VTK_QUAD
    else:
#    elif type(mesh) in [MeshTools.TetMesh, MeshTools.HexMesh]:
        kernel.init_warmup(runtime.logfile)
        if mpi.is_on_master_proc:
            vertices = MT.as_coordinate_array(mesh.vertices)
            cells_nodes, cells_faces, faces_nodes = mesh.COC_data()
            kernel.create_mesh(vertices,
                               *cells_nodes,
                               *cells_faces,
                               *faces_nodes
                               )
            # distribute cell types to be able to reconstruct local meshes
            celltypes = ComPASS.global_celltypes()
            celltypes[:] = mesh.cells_vtk_ids().astype(np.int8, copy=False)
            facetypes = ComPASS.global_facetypes()
            facetypes[:] = mesh.faces_vtk_ids().astype(np.int8, copy=False)
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
        else:
            dirichlet = set_pressure_dirichlet_nodes()
            if dirichlet is not None:
                info.pressure[dirichlet] = ord('d')
            dirichlet = set_temperature_dirichlet_nodes()
            if dirichlet is not None:
                info.temperature[dirichlet] = ord('d')
        kernel.global_mesh_count_dirichlet_nodes()
        kernel.global_mesh_frac_by_node()
        # The following line is necessary to allocate arrays in the fortran code
        kernel.global_mesh_allocate_rocktype()
        if set_global_rocktype is not None:
            assert callable(set_global_rocktype)
            set_global_rocktype()
        kernel.global_mesh_make_post_read_set_poroperm()
        for location in ['cell', 'fracture']:
            for property in ['porosity', 'permeability']:
                value = properties[location + '_' + property]()
                if value is not None:
                    dim = 3
                    if location=='fracture':
                        assert fractures is not None
                        dim = 2
                    buffer = np.array(getattr(kernel, 'get_%s_%s_buffer' % (location, property))(), copy = False)
                    value = np.ascontiguousarray( value )
                    # CHECKME: fracture permeability is scalar
                    if property=='permeability' and location=='cell':
                        n = buffer.shape[0]
                        if value.shape==(1,): # scalar value
                            value = np.tile(value[0] * np.eye(dim), (n, 1 ,1)) 
                        if value.shape==(dim, dim): # tensor value
                            value = np.tile(value, (n, 1, 1))
                        assert value.shape==(n, dim, dim)
                        assert buffer.shape==value.shape
                    buffer[:] = value
        kernel.global_mesh_make_post_read_well_connectivity_and_ip()
        kernel.set_well_data(well_list)
        kernel.compute_well_indices()
    kernel.init_phase2(runtime.output_directory)
    mpi.synchronize() # wait for every process to synchronize
    # FUTURE: This could be managed through a context manager ?
    initialized = True
    atexit.register(finalize)

def get_global_id_faces():
   return np.array(kernel.get_global_id_faces_buffer(), copy = False)

def get_cell_permeability():
   return np.array(kernel.get_cell_permeability_buffer(), copy = False)

def get_fracture_permeability():
   return np.array(kernel.get_fracture_permeability_buffer(), copy = False)

def get_cell_porosity():
   return np.array(kernel.get_cell_porosity_buffer(), copy = False)

def get_fracture_porosity():
   return np.array(kernel.get_fracture_porosity_buffer(), copy = False)

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
    centers = ComPASS.face_centers()
    dtype = set(centers.dtype[k] for k in range(len(centers.dtype)))
    centers = centers.view(dtype=dtype.pop()).reshape((-1, 3))
    assert not dtype # dtype should be empty here
    return centers

def compute_fracture_centers():
    return compute_face_centers()[ComPASS.frac_face_id() -1] # Fortran indexes start at 1

def old_compute_face_centers():
    return _compute_centers(
        vertices().view(dtype=np.double).reshape((-1, 3)),
        get_connectivity().NodebyFace
    )

def old_compute_fracture_centers():
    return _compute_centers(
        vertices().view(dtype=np.double).reshape((-1, 3)),
        ComPASS.get_nodes_by_fractures()
    )

def compute_global_face_normals():
    vertices = global_vertices().view(dtype=np.double).reshape((-1, 3))
    connectivity = get_global_connectivity()
    face_nodes = [np.array(nodes, copy=False)[:3] - 1 # fortran indexes start at 1
                  for nodes in connectivity.NodebyFace] 
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
    idfaces = get_global_id_faces()
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

def add_output_subdirectory(path):
    target = ComPASS.to_output_directory(path)
    if mpi.is_on_master_proc:
       if not os.path.exists(target):
            os.makedirs(target)
    mpi.synchronize()
    assert os.path.isdir(target)

def find_fracture_edges(faces):
    faces = np.asarray(faces)
    if faces.dtype==np.bool:
        faces = np.nonzero(faces)[0]
    face_nodes = get_connectivity().NodebyFace
    # we do not want to store twice the same edge
    fracture_edges = set()
    info = np.rec.array(node_info(), copy=False)
    for face in faces:
        nodes = np.array(face_nodes[face], copy=False) - 1 # Fortran indexing...
        in_fracture = info.frac[nodes] == ord('y')
        for i in range(nodes.shape[0]):
            if in_fracture[i-1] and in_fracture[i]:
                if nodes[i-1]<nodes[i]:
                    fracture_edges.add((nodes[i-1], nodes[i]))
                else:
                    fracture_edges.add((nodes[i], nodes[i-1]))
    return  np.array(list(fracture_edges))

def set_Neumann_faces(faces, Neumann):
    faces = np.asarray(faces)
    if faces.dtype==np.bool:
        faces = np.nonzero(faces)[0]
    faces+= 1 # Fortran indexing starts at 1   
    kernel.set_Neumann_faces(faces, Neumann)

def set_Neumann_fracture_edges(edges, Neumann):
    edges = np.asarray(edges)
    edges+= 1 # Fortran indexing starts at 1   
    kernel.set_Neumann_fracture_edges(edges, Neumann)