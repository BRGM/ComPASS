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
    grid = grid
)

@ComPASS.mpi.on_master_proc
def display_own_elements_number():
    print("cells own:", ComPASS.kernel.nb_cells_own())
    print("faces own:", ComPASS.kernel.nb_faces_own())
    print("nodes own:", ComPASS.kernel.nb_nodes_own())
    print("fractures own:", ComPASS.kernel.nb_fractures_own())

display_own_elements_number()
