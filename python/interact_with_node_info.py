import numpy as np

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

def array_description(a):
    return {val: np.count_nonzero(a==val) for val in np.unique(a)} 

def output_node_info(info):
    # This a way to access field as attibutes
    info = np.rec.array(info)
    print('Nb of nodes:', info.proc.shape[0])
    print(array_description(info.proc.view('c')))
    print(array_description(info.frac.view('c')))
    print(array_description(info.pressure.view('c')))
    print(array_description(info.temperature.view('c')))

@ComPASS.on_master_proc
def output_global_node_info():
    print('Proc info (own/ghost) is not relevant at global scale.') 
    info = output_node_info(ComPASS.global_node_info())

def output_local_node_info():
    print('Node info for proc:', ComPASS.proc_rank) 
    info = output_node_info(ComPASS.node_info())

# The following line will crash as the mesh is distributed at the end of ComPASS.init
# consequently the global node info array is deallocated
# output_global_node_info()

output_local_node_info()

ComPASS.finalize()
