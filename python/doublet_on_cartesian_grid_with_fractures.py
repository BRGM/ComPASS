import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

import numpy as np

def fractures_factory(grid):
    def select_fractures():
        fractures_centers = ComPASS.compute_face_centers()
        dz = grid.extent[2] / grid.shape[2]
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = grid.origin[2] + 0.5 * grid.extent[2]
        return np.abs(fractures_centers[:, 2] - zfrac) < 0.5 * dz
    return select_fractures

grid = ComPASS.Grid(
    origin = (-1500., -1000., -1600.),
    extent = (3000., 2000., 100.),
    shape = (31, 21, 11)
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
    fracture_faces = fractures_factory(grid)
)

standard_loop(final_time = 30 * year, output_frequency = year)
