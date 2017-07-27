import numpy as np
import ComPASS.ComPASS.MeshTools as MT
import ComPASS.GridTools as GT
import vtkwriters as vtkw

vertices, cells = GT.grid2hexs((3, 2, 4))

cells = np.ascontiguousarray(cells, dtype=MT.idtype())

mesh = MT.hexmesh(vertices, cells)

vtkw.write_vtu(
    vtkw.vtu_doc(mesh.vertices, mesh.cellnodes),
    'hexs.vtu'
)
