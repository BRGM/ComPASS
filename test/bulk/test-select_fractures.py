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


simulation = ComPASS.load_eos("water2ph")


def fractures_factory(grid):
    def select_fractures():
        face_centers = simulation.compute_global_face_centers()
        dz = grid.extent[2] / grid.shape[2]
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = grid.origin[2] + 0.5 * grid.extent[2]
        return np.abs(face_centers[:, 2] - zfrac) < 0.25 * dz

    return select_fractures


grid = ComPASS.Grid(
    origin=(-1500.0, -1000.0, -1600.0), extent=(3000.0, 2000.0, 100.0), shape=(2, 2, 2)
)

ComPASS.set_output_directory_and_logfile(__file__)

simulation.init(
    mesh=grid,
    fracture_faces=fractures_factory(grid),
    cell_porosity=0.5,  # dummy value, no simulation
    cell_permeability=1e-12,  # dummy value, no simulation
    cell_thermal_conductivity=2,  # dummy value, no simulation
    fracture_porosity=0.5,  # dummy value, no simulation
    fracture_permeability=1e-12,  # dummy value, no simulation
    fracture_thermal_conductivity=2,  # dummy value, no simulation
)

ComPASS.dumps.dump_mesh(simulation)
