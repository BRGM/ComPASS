#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np
import MeshTools as MT
import MeshTools.GridTools as GT
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop


def test_nodeflags_on_tetmesh():

    simulation = ComPASS.load_eos('water2ph')

    gridshape = (1, 1, 1)
    gridextent = (1E3, 1E3, 1E3)
    vertices, tets = GT.grid2tets(gridshape, gridextent)
    mesh = MT.TetMesh.make(vertices, MT.idarray(tets))

    def set_global_nodeflags():
        flags = simulation.global_nodeflags()
        n = flags.shape[0]
        assert vertices.shape[0]==n
        flags[:] = np.arange(n)

    def check_nodeflags():
        local_vertices = simulation.vertices()
        flags = simulation.nodeflags()
        return np.all(vertices[flags]==local_vertices)

    ComPASS.set_output_directory_and_logfile(__file__)

    simulation.init(
        mesh = mesh,
        set_global_flags = set_global_nodeflags,
        cell_porosity = 0.5,               # dummy value, no simulation
        cell_permeability = 1E-12,         # dummy value, no simulation
        cell_thermal_conductivity = 2,     # dummy value, no simulation
    )

    assert check_nodeflags()

if __name__ == '__main__':
    test_nodeflags_on_tetmesh()
