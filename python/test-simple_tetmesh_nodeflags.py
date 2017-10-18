#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.GridTools as GT
MT = ComPASS.ComPASS.MeshTools
#import vtkwriters as vtkw

gridshape = (1, 1, 1)
gridextent = (1E3, 1E3, 1E3)
vertices, tets = GT.grid2tets(gridshape, gridextent)
tets = np.asarray(tets, dtype=MT.idtype())

mesh = MT.tetmesh(vertices, tets)

def set_global_nodeflags():
    flags = ComPASS.global_nodeflags()
    n = flags.shape[0]
    assert vertices.shape[0]==n
    flags[:] = np.arange(n)

def check_nodeflags():
    local_vertices = ComPASS.vertices().view(np.double).reshape((-1, 3))
    flags = ComPASS.nodeflags()
    result = np.all(vertices[flags]==local_vertices)
    assert result
    return result

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = mesh,
    set_global_flags = set_global_nodeflags
)

if not check_nodeflags():
    print('Bad distribution of node flags!')
