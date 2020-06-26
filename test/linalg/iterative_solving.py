from doublet_on_cartesian_grid_setup import *

"""
This script runs the default time loop with with a Petsc iterative solver
"""

lsolver = PetscIterativeSolver(
    PetscLinearSystem(simulation), IterativeSolverSettings(1.0e-6, 150, 30)
)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    newton=newton,
    initial_timestep=30 * day,
    final_time=100 * day,
    output_period=year,
    context=context,
)
