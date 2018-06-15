#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS

ComPASS.set_output_directory_and_logfile(__file__)
ComPASS.load_eos('water2ph')

nx, ny, nz = 3, 2, 1

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
)

def dirichlet_nodes():
    vertices = ComPASS.global_vertices()
    newz = 10.
    z = vertices[:, 2]
    z[z!=0] = newz
    assert np.all(np.logical_or(ComPASS.global_vertices()[:, 2] == 0, ComPASS.global_vertices()[:, 2] == newz))

ComPASS.init(
    mesh = grid,
    set_dirichlet_nodes = dirichlet_nodes,
    cell_porosity = 0.1,
    cell_permeability = 1.E-12,
    cell_thermal_conductivity = 2.,
)
