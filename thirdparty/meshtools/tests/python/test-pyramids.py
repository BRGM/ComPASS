import numpy as np
import MeshTools as MT
import vtkwriters as vtkw

vertices = np.array([
    (0, 0, 0),
    (2, 0, 0),
    (2, 2, 0),
    (0, 2, 0),
    (1, 1, 1),
    (1, 1, -1),
], dtype='d')

mesh = MT.HybridMesh.Mesh()

mesh.set_vertices(vertices)
cellnodes = mesh.connectivity.cells.nodes
Pyramid = MT.Pyramid
cellnodes.append( Pyramid((0, 1, 2, 3, 4)) )
cellnodes.append( Pyramid((0, 1, 2, 3, 5)) )
mesh.connectivity.update_from_cellnodes()

offsets, cellsnodes = mesh.cells_nodes_as_COC()
vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        mesh.vertices_array(), 
        np.array(offsets[1:], copy=False), # no first zero offset for vtk 
        np.array(cellsnodes, copy=False),
        mesh.cells_vtk_ids(),
    ),
    'pyramids.vtu'
)
