import numpy as np
import MeshTools as MT
import MeshTools.GridTools as GT
import MeshTools.vtkwriters as vtkw

def test_hexmesh():

    vertices, cells = GT.grid2hexs(shape=(3, 2, 4))

    cells = np.ascontiguousarray(cells, dtype=MT.idtype())

    mesh = MT.HexMesh.make(vertices, cells)

    vertices = MT.as_coordinate_array(mesh.vertices)
    cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

    vtkw.write_vtu(
        vtkw.vtu_doc(vertices, cellnodes),
        'hexs.vtu'
    )

if __name__=='__main__':
    test_hexmesh()
