import ComPASS
import doublet_utils
from ComPASS.utils.units import *

ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (31, 21, 1),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
)

def output_states(X):
    print(X.context[:2])
    print(X.p[:2])
    print(X.T[:2])
    print(X.C[:2])
    print(X.S[:2])
    print(X.accumulation[:2])

@ComPASS.on_master_proc
def output_all_states():
    print('*'*5, 'boundary nodes')
    output_states(ComPASS.boundary_node_states())
    print('*'*5, 'nodes')
    output_states(ComPASS.node_states())
    print('*'*5, 'fractures')
    output_states(ComPASS.fracture_states())
    print('*'*5, 'cells')
    output_states(ComPASS.cell_states())

output_all_states()

ComPASS.finalize()
