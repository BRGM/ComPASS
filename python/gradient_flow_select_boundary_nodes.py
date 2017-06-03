import ComPASS
import doublet_utils
from ComPASS.utils.units import *
np = ComPASS.np

pleft, pright = 30 * MPa, 10 * MPa
p0 = pright
Tleft, Tright = degC2K(60), degC2K(100)
T0 = Tright

grid = ComPASS.Grid(
    #shape = (121, 11, 1),
    shape = (31, 11, 1),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

on_the_left = lambda pts: pts.x <= grid.origin[0]
on_the_right = lambda pts: pts.x >= grid.origin[0] + grid.extent[0]

def select_dirichlet_nodes():
    vertices = np.rec.array(ComPASS.global_vertices(), copy=False)
    return on_the_left(vertices) | on_the_right(vertices)

def set_boundary_conditions():
    vertices = np.rec.array(ComPASS.vertices(), copy=False)
    #print('>'*10, 'on proc', ComPASS.proc_rank, ':\n', vertices.x.min(), '< X <', vertices.x.max(), '\n', vertices.y.min(), '< Y <', vertices.y.max())
    dirichlet = ComPASS.dirichlet_node_states()
    left = on_the_left(vertices)
    dirichlet.p[left] = pleft 
    dirichlet.T[left] = Tleft 
    right = on_the_right(vertices)
    dirichlet.p[right] = pright
    dirichlet.T[right] = Tright
    dirichlet.context[:] = 2
    dirichlet.S[:] = [0, 1]
    dirichlet.C[:] = 1.

def set_initial_values():
    for state in [ComPASS.node_states(), ComPASS.fracture_states(), ComPASS.cell_states()]:
        state.context[:] = 2
        state.p[:] = p0
        state.T[:] = T0
        state.S[:] = [0, 1]
        state.C[:] = 1.

@ComPASS.on_master_proc
def print_iteration_info():
    print()
    print('Time Step (iteration):', n)
    print('Current time: %.1f years' % (ComPASS.get_current_time() / year), ' -> final time:', final_time / year)
    print('Timestep: %.3f days' % (ComPASS.get_timestep() / day))

# %%% Simulation %%% 

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    set_dirichlet_nodes = select_dirichlet_nodes
)

set_boundary_conditions()
set_initial_values()

final_time = 30 * year
output_frequency = 1 * year
nitermax = 1E10
t = 0
n = 0
t_output = 0
while t <= final_time and n < nitermax:
    if t >= t_output:
        ComPASS.output_visualization_files(n)
        # WARNING / CHECKME we may loose some outputs
        while (t_output < t):
            t_output = t_output + output_frequency
    n+= 1
    print_iteration_info()
    ComPASS.make_timestep()
    t = ComPASS.get_current_time()
    ComPASS.timestep_summary()
# Output final time
ComPASS.output_visualization_files(n)

ComPASS.finalize()
