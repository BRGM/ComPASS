""" simulation utils : "pure" functions
"""


import numpy as np


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
