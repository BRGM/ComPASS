import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.load_eos('water2ph')

pres = 10. * MPa
Tres = 70. # degC

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    set_dirichlet_nodes = doublet_utils.select_boundary_factory(grid),
    wells = doublet_utils.make_wells_factory(grid),
)
doublet_utils.init_states(pres, Tres)

standard_loop(final_time = 30 * year, output_frequency = year)
