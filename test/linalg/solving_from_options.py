from doublet_on_cartesian_grid_setup import *

from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver

"""
This script runs the default time loop with a linear solver set from command line options
"""

command_line_lsolver = linear_solver(simulation, from_options=True)
newton = Newton(simulation, 1e-5, 8, command_line_lsolver)

simulation.standard_loop(
    initial_timestep=30 * day, final_time=100 * day, output_period=year, newton=newton,
)

# You can try running this script with different options, for example this will setup a Petsc direct solver and display a short view :

# python3 solving_from_options.py --legacy_linear_solver False --direct_linear_solver True --linear_solver_view True
