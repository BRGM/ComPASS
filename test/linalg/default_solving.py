from doublet_on_cartesian_grid_setup import *

"""
This script runs the default time loop with no specific solver settings
"""


simulation.standard_loop(
    initial_timestep=30 * day,
    final_time=100 * day,
    output_period=year,
    context=context,
)
