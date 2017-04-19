import numpy as np

from .ComPASS import *

def get_vertices():
   return np.array(ComPASS.get_vertices_buffer(), copy = False)

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



