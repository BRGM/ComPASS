#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import MeshTools as MT
import GridTools as GT
import vtkwriters as vtkw

vertices, cells = GT.grid2hexs((3, 2, 4))

cells = np.ascontiguousarray(cells, dtype=MT.idtype())

mesh = MT.HexMesh.make(vertices, cells)

vertices = MT.as_coordinate_array(mesh.vertices)
cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, cellnodes),
    'hexs.vtu'
)
