#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.set_output_directory_and_logfile(__file__)

def cell_permeability_factory(grid):
    Ox, Oy, Oz = grid.origin
    Lx, Ly, Lz = grid.extent
    def cell_permeability():
        cell_centers = ComPASS.compute_cell_centers()
        zc = cell_centers[:, 2]
        nbcells = cell_centers.shape[0]
        # tensor array
        cellperm = np.empty((nbcells, 3, 3), dtype=np.double)
        # matrix permeability
        cellperm[:] = 1E-16 * np.eye(3)
        # anisotropic reservoir permeability
        zc-= Oz
        reservoir = (zc > Lz/3.) & (zc < 2*Lz/3.)
        kres = 1E-12 * np.array([
            [np.sqrt(2), np.sqrt(2), 0],
            [-np.sqrt(2), np.sqrt(2), 0],
            [0, 0, 0.1]], dtype=np.double)
        cellperm[reservoir] = kres
        return cellperm
    return cell_permeability

grid = ComPASS.Grid(
    shape = (61, 41, 12),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
    cells_permeability = cell_permeability_factory(grid),
)

standard_loop(final_time = 30 * year, output_frequency = year)
