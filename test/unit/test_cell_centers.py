#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS


def test():
    nx, ny, nz = 1, 1, 1

    ComPASS.load_physics("water2ph")
    ComPASS.set_output_directory_and_logfile(__file__)
    grid = ComPASS.Grid(shape=(nx, ny, nz))
    ComPASS.init_and_load_mesh(grid)

    assert np.allclose(ComPASS.compute_global_cell_centers(), [[0.5, 0.5, 0.5]])
    assert np.allclose(
        ComPASS.compute_global_face_centers(),
        [
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 1.0],
            [1.0, 0.5, 0.5],
            [0.5, 1.0, 0.5],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ],
    )
