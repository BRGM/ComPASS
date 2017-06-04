import numpy as np
import ComPASS
import doublet_utils
from ComPASS.utils.units import *

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
