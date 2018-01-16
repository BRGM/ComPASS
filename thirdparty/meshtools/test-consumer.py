import numpy as np
import MeshTools as MT
import GridTools as GT
import vtkwriters as vtkw
import consumer

vertices, cells = GT.grid2hexs((3, 2, 4))

cells = np.ascontiguousarray(cells, dtype=MT.idtype())

mesh = MT.HexMesh.make(vertices, cells)

vertices, cells = mesh.cells_nodes_as_COC()

consumer.consume_vector(cells)