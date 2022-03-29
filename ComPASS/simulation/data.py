""" simulation data access
"""

import warnings
import numpy as np

from .utils import (
    compute_centers,
    get_boundary_faces as _get_boundary_faces,
    get_boundary_vertices as _get_boundary_vertices,
)
from .._kernel import simulation_wrapper as _sw

# CHECKME: does (re)loading module change _simulation_object.self values?
from . import _simulation_object


def collect_all_edges():
    assert not _simulation_object.self.mesh_is_local, "mesh is assumed to be global"
    face_nodes = _sw.get_global_connectivity().NodebyFace
    edges = []
    for face in face_nodes:
        nodes = np.array(face, copy=False)
        for i in range(nodes.shape[0]):
            edges.append((nodes[i - 1], nodes[i]))
    edges = np.array(edges)
    assert edges.ndim == 2 and edges.shape[1] == 2
    edges = np.sort(edges, axis=1)
    return np.unique(edges, axis=0)


def get_global_id_faces():
    return np.array(_sw.get_global_id_faces_buffer(), copy=False)


def get_cell_permeability():
    warnings.warn("Use simulation.petropyhics()", DeprecationWarning)
    return _sw.petrophysics().cell_permeability


def get_fracture_permeability():
    warnings.warn("Use simulation.petropyhics()", DeprecationWarning)
    return _sw.petrophysics().fracture_permeability


def get_cell_porosity():
    warnings.warn("Use simulation.petropyhics()", DeprecationWarning)
    return _sw.petrophysics().cell_porosity


def get_fracture_porosity():
    warnings.warn("Use simulation.petropyhics()", DeprecationWarning)
    return _sw.petrophysics().fracture_porosity


def get_cell_heat_source():
    return np.array(_sw.get_cell_heat_source_buffer(), copy=False)


def compute_global_cell_centers():
    return compute_centers(
        _sw.global_vertices().view(dtype=np.double).reshape((-1, 3)),
        _sw.get_global_connectivity().NodebyCell,
    )


def compute_global_face_centers():
    return compute_centers(
        _sw.global_vertices().view(dtype=np.double).reshape((-1, 3)),
        _sw.get_global_connectivity().NodebyFace,
    )


def compute_cell_centers():
    return compute_centers(
        _sw.vertices().view(dtype=np.double).reshape((-1, 3)),
        _sw.get_connectivity().NodebyCell,
    )


def compute_face_centers():
    # centers = _sw.face_centers()
    # dtype = set(centers.dtype[k] for k in range(len(centers.dtype)))
    # centers = centers.view(dtype=dtype.pop()).reshape((-1, 3))
    # assert not dtype # dtype should be empty here
    return _sw.face_centers()


def compute_fracture_centers():
    return compute_face_centers()[_sw.frac_face_id() - 1]  # Fortran indexes start at 1


def compute_dof_locations():
    import warnings

    warnings.warn(
        DeprecationWarning(
            '''"compute_dof_locations()" is deprecated and will be removed soon, use "all_positions() instead."'''
        )
    )
    return all_positions()


def all_positions():
    """Returns all position of degrees of freedom stacked in the same order
    as states.
    """
    return np.vstack(
        [_sw.vertices(), compute_fracture_centers(), compute_cell_centers()]
    )


def old_compute_face_centers():
    return compute_centers(
        _sw.vertices().view(dtype=np.double).reshape((-1, 3)),
        _sw.get_connectivity().NodebyFace,
    )


def old_compute_fracture_centers():
    return compute_centers(
        _sw.vertices().view(dtype=np.double).reshape((-1, 3)),
        _sw.get_nodes_by_fractures(),
    )


def compute_global_face_normals():
    vertices = _sw.global_vertices().view(dtype=np.double).reshape((-1, 3))
    connectivity = _sw.get_global_connectivity()
    face_nodes = [
        np.array(nodes, copy=False)[:3] - 1  # fortran indexes start at 1
        for nodes in connectivity.NodebyFace
    ]
    normals = np.array(
        [
            np.cross(
                vertices[nodes[1]] - vertices[nodes[0]],
                vertices[nodes[2]] - vertices[nodes[0]],
            )
            for nodes in face_nodes
        ]
    )
    # normalize
    norms = np.linalg.norm(normals, axis=1)
    norms.shape = (-1, 1)
    normals /= norms
    return normals


def get_global_boundary_faces():
    return _get_boundary_faces(_sw.get_global_connectivity())


def get_boundary_faces():
    return _get_boundary_faces(_sw.get_connectivity())


def get_global_boundary_vertices():
    return _get_boundary_vertices(_sw.get_global_connectivity())


def get_boundary_vertices():
    return _get_boundary_vertices(_sw.get_connectivity())


def set_fractures(faces):
    idfaces = get_global_id_faces()
    assert faces.shape == idfaces.shape or (
        faces.max() < idfaces.shape[0] and faces.min() >= 0
    )
    idfaces[faces] = -2
    _sw.global_mesh_set_frac()  # this will collect faces with flag -2 as fracture faces


def set_Neumann_faces(faces, Neumann):
    faces = np.asarray(faces)
    if faces.dtype == np.bool:
        faces = np.nonzero(faces)[0]
    faces += 1  # Fortran indexing starts at 1
    _sw.set_Neumann_faces(faces, Neumann)


def set_Neumann_fracture_edges(edges, Neumann):
    edges = np.asarray(edges)
    edges += 1  # Fortran indexing starts at 1
    _sw.set_Neumann_fracture_edges(edges, Neumann)
