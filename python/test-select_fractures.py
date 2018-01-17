#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop


ComPASS.load_eos('water2ph')

def fractures_factory(grid):
    def select_fractures():
        face_centers = ComPASS.compute_global_face_centers()
        dz = grid.extent[2] / grid.shape[2]
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = grid.origin[2] + 0.5 * grid.extent[2]
        return np.abs(face_centers[:, 2] - zfrac) < 0.5 * dz
    return select_fractures

grid = ComPASS.Grid(
    origin = (-1500., -1000., -1600.),
    extent = (3000., 2000., 100.),
    shape = (2, 2, 2)
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    fracture_faces = fractures_factory(grid)
)

ComPASS.dumps.dump_mesh()
