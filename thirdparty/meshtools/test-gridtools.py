import numpy as np
import GridTools as GT
import MeshTools as MT

import vtkwriters as vtkw

shape = (1, 1, 1)

vertices, tets = GT.grid2tets(shape)

mesh = MT.TetMesh.make(vertices, MT.idarray(tets))

print('Nb vertices:', mesh.nb_vertices)
print('Nb cells:', mesh.nb_cells)
print('Nb faces:', mesh.nb_faces)

def print_collection(iterable):
    for item in iterable:
        print(item)

print('Vertices')
print_collection(mesh.vertices)
print('Cell face')
print_collection(mesh.connectivity.cells.faces)
print('Cell nodes')
print_collection(mesh.connectivity.cells.nodes)
print('Face nodes')
print_collection(mesh.connectivity.faces.nodes)

print('All cell centers:')
print_collection(mesh.cell_centers())

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, tets, celldata={'id': np.arange(1, tets.shape[0]+1)}),
    'cubes.vtu'
)



