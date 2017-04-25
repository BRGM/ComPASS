import ComPASS
import doublet_utils
from ComPASS.utils.units import *

ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
)

@ComPASS.on_master_proc
def print_iteration_info():
    print()
    print('Time Step (iteration):', n)
    print('Current time: %.1f years' % (ComPASS.get_current_time() / year), ' -> final time:', final_time / year)
    print('Timestep: %.3f days' % (ComPASS.get_timestep() / day))
    
final_time = 30 * year
n = 0
output_frequency = 1 * year
t_output = 0
while ComPASS.get_current_time() <= final_time:
    n+= 1
    print_iteration_info()
    ComPASS.make_timestep()
    t = ComPASS.get_current_time()
    if t > t_output:
        ComPASS.output_visualization_files(n)
    # WARNING / CHECKME we may loose some outputs
    while (t_output < t):
        t_output = t_output + output_frequency
    ComPASS.timestep_summary()

ComPASS.finalize()
