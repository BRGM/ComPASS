import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
)

standard_loop(final_time = 30 * year, output_frequency = year)
