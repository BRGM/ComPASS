from doublet_on_cartesian_grid_setup import *

from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver

"""
This script runs the default time loop with with a Legacy iterative solver
"""

lsolver = linear_solver(simulation, legacy=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    initial_timestep=30 * day,
    final_time=100 * day,
    output_period=year,
    newton=newton,
    context=context,
)

print(lsolver)
