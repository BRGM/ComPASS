""" simulation utils : "pure" functions
"""

import re
import numpy as np
from .. import mpi
from ..postprocess import postprocess as postprocess_command
from ..exceptions import CompassException


def reshape_as_scalar_array(value, n):
    value = np.ascontiguousarray(value)
    if value.shape == (1,):  # scalar value
        value = np.tile(value[0], n)
    assert value.shape == (n,)
    return value


def reshape_as_tensor_array(value, n, dim):
    value = np.ascontiguousarray(value)
    if value.shape == (1,):  # scalar value
        value = np.tile(value[0] * np.eye(dim), (n, 1, 1))
    if value.shape == (dim, dim):  # tensor value
        value = np.tile(value, (n, 1, 1))
    if value.shape == (n,):  # scalar values
        value = value[:, None, None] * np.eye(dim)
    assert value.shape == (n, dim, dim)
    return value


def compute_centers(points, elements):
    return np.array(
        [
            points[
                np.array(element, copy=False) - 1  # fortran indexes start at 1
            ].mean(axis=0)
            for element in elements
        ]
    )


def get_boundary_faces(connectivity):
    return np.array(
        [
            np.array(face_cells, copy=False).shape[0] == 1
            for face_cells in connectivity.CellbyFace
        ]
    )


def get_boundary_vertices(connectivity):
    boundary_faces = get_boundary_faces(connectivity)
    vertices_id = np.unique(
        np.hstack(
            [
                np.array(face_nodes, copy=False)
                for boundary, face_nodes in zip(boundary_faces, connectivity.NodebyFace)
                if boundary
            ]
        )
    )
    vertices_id -= 1  # Fortran index...
    nbnodes = len(connectivity.CellbyNode)
    res = np.zeros(nbnodes, dtype=np.bool)
    res[vertices_id] = True
    return res


def coordinates(a):
    assert len(a.shape) == 2 and a.shape[1] == 3
    return (a[:, j] for j in range(a.shape[1]))


def get_faces_nodes(connectivity, faces_selection):
    selection = np.asarray(faces_selection)
    assert selection.ndim == 1
    if selection.dtype == np.bool:
        selection = np.nonzero(selection)[0]
    NodebyFace = connectivity.NodebyFace
    assert len(NodebyFace) > 0
    # Fortran -> C indexing
    faces = [np.array(NodebyFace[f], copy=False) - 1 for f in selection]
    return faces


def facenodes(simulation, faces_selection):
    assert simulation.is_local
    faces = get_faces_nodes(simulation.get_connectivity(), faces_selection)
    return np.unique(faces)


def postprocess(simulation, **kwargs):
    mpi.synchronize()
    if mpi.is_on_master_proc:
        postprocess_command(simulation.runtime.output_directory, **kwargs)


def physics_name(simulation):
    kernel = simulation.get_kernel()
    m = re.match("ComPASS\.physics\.(\w*)", kernel.__name__)
    if not m:
        raise CompassException("Not a valid ComPASS physics!")
    assert len(m.groups()) == 1
    return m.groups()[0]


def collect_all_edges(simulation):
    assert not simulation.mesh_is_local, "mesh is assumed to be global"
    face_nodes = simulation.get_global_connectivity().NodebyFace
    edges = []
    for face in face_nodes:
        nodes = np.array(face, copy=False)
        for i in range(nodes.shape[0]):
            edges.append((nodes[i - 1], nodes[i]))
    edges = np.array(edges)
    assert edges.ndim == 2 and edges.shape[1] == 2
    edges = np.sort(edges, axis=1)
    return np.unique(edges, axis=0)
