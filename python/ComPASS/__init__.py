import numpy as np

from .ComPASS import *

def get_vertices():
   return np.array(ComPASS.get_vertices_buffer(), copy = False)

def get_id_faces():
   return np.array(ComPASS.get_id_faces_buffer(), copy = False)

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


