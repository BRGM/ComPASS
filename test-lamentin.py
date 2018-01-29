import os
import numpy as np
import MeshTools as MT

import vtkwriters as vtkw

dirname = 'Lamentin'
def retrieve(filename, dtype):
    return np.loadtxt(os.path.join(dirname, filename), dtype=dtype)

vertices = retrieve('vertices', np.double)
cells = retrieve('cells', MT.idtype())
faces = retrieve('faces', MT.idtype())
fault_faces = retrieve('faces_fault', MT.idtype())
fault_faces_nodes = faces[fault_faces]

mesh = MT.tetmesh(vertices, cells)
fault_faces_id = mesh.faces_ids(fault_faces_nodes)

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, cells),
    os.path.join(dirname, 'mesh.vtu')
)

ffcenters = mesh.face_centers(fault_faces_id)

vtkw.write_vtu(
    vtkw.vtu_doc(ffcenters, np.reshape(np.arange(ffcenters.shape[0]), (-1,1))),
    os.path.join(dirname, 'fault_faces_centers.vtu')
)

fcenters = mesh.face_centers(mesh.faces_ids(faces))
patches = retrieve('patches', np.int)
assert fcenters.shape[0]==patches.shape[0]

vtkw.write_vtu(
    vtkw.vtu_doc(fcenters, np.reshape(np.arange(fcenters.shape[0]), (-1,1)), pointdata={'faceid': patches}),
    os.path.join(dirname, 'facets_centers.vtu')
)





