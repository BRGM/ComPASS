# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 09:15:07 2018

@author: lopez
"""

import sys
import pickle
import numpy as np
import MeshTools as MT
from MeshTools.io.medit import reader


def process(filename):

    L = 3e3  # target cube extent

    mesh_info = reader(filename)

    vertices = mesh_info.vertices

    vmin, vmax = vertices.min(), vertices.max()
    vertices = (L / (vmax - vmin)) * (vertices - vmin)
    # set top nodes at 0 elevation
    vertices[:, 2] -= vertices[:, 2].max()

    mesh = MT.TetMesh.make(vertices, MT.idarray(mesh_info.tets))

    special_faces = np.array(
        [
            mesh.connectivity.faces.id(MT.Triangle(MT.idarray(face)))
            for face in mesh_info.triangles
        ]
    )

    zmin = vertices[:, 2].min()
    fracture_nodes = np.unique(mesh_info.triangles)
    bottom_fracture_nodes = fracture_nodes[vertices[fracture_nodes, 2] <= zmin]

    with open(filename + ".pkl", "wb") as f:
        pickle.dump(
            {
                "mesh": mesh,
                "special_faces": special_faces,
                "special_faces_flags": mesh_info.triangles_index,
                "cells_flags": mesh_info.tets_index,
                "bottom_fracture_nodes": bottom_fracture_nodes,
            },
            f,
        )


for filename in sys.argv[1:]:
    process(filename)
