# Documentation : https://charms.gitlabpages.inria.fr/ComPASS/
import numpy as np
import ComPASS
from ComPASS.utils.units import *  # contains MPa, degC2K, km, year...
from ComPASS.utils.grid import on_zmin, on_zmax

#  |------------------------------------|
#  |            Overburden              |
#  |------------------------------------|
#  |                                    |
#  |                                    |
#  |             Reservoir              |
#  |                                    |
#  |                                    |
#  |------------------------------------|
#  |            Underburden             |
#  |------------------------------------|


# Geometry
total_depth = 3 * km
Hburden = 0.1 * total_depth


# Reservoir petrophysics
k_reservoir = 1e-12  # reservoir permeability in m^2
omega = 0.15  # reservoir and burden porosity
K_reservoir = 2  # reservoir and burden bulk thermal conductivity in W/m/K


# Burden petrophysics
k_burden = 1e-18  # burden permeability in m^2


# -------------------------------------------------------------------
# Set output informations
ComPASS.set_output_directory_and_logfile(__file__)
# Load the water2ph physics : it contains the water component
# which can be in liquid and/or gas phase
simulation = ComPASS.load_physics("water2ph")


# -------------------------------------------------------------------
# Create a Cartesian grid with cubic cells
nx = nz = 51
ny = 1
d = total_depth / nz
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(nx * d, ny * d, nz * d),  # creates cubic cells
    origin=(0.0, 0.0, -total_depth),
)


# -------------------------------------------------------------------
# Identify the Dirichlet nodes
def dirichlet_boundaries_function():
    pts = simulation.global_vertices()
    return on_zmax(grid)(pts) | on_zmin(grid)(pts)


# -------------------------------------------------------------------
# Initialize the regionalized values and distribute the domain
def cell_permeability_function():
    cell_centers = simulation.compute_global_cell_centers()
    zc = cell_centers[:, 2]
    nbcells = cell_centers.shape[0]
    permeability = np.empty(nbcells, dtype=np.double)
    permeability[:] = k_burden
    reservoir = (zc < -Hburden) & (zc > -total_depth + Hburden)
    permeability[reservoir] = k_reservoir
    return permeability


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=dirichlet_boundaries_function,
    cell_porosity=omega,
    cell_thermal_conductivity=K_reservoir,
    cell_permeability=cell_permeability_function,
)


# -------------------------------------------------------------------
# Initialize the domain with liquid phase
# at bottom_pressure pressure and bottom_temperature temperature.
top_temperature = degC2K(25.0)
top_pressure = 1 * bar
bottom_temperature = degC2K(180.0)
bottom_pressure = top_pressure + simulation.get_gravity() * 900.0 * total_depth
# First construct the state (at thermodynamic equilibrium)
X0 = simulation.build_state(
    simulation.Context.liquid, p=bottom_pressure, T=bottom_temperature
)
# then apply this state everywhere to init the domain (constant values)
simulation.all_states().set(X0)


# Modify the initial values to apply a linear gradient
# in pressure and temperature
def set_linear_gradients_state(states, depth):
    # depth must be positive
    def linear_gradients(bottom_value, top_value, domain_depth, depth):
        return top_value + (bottom_value - top_value) / domain_depth * depth

    states.p[:] = linear_gradients(bottom_pressure, top_pressure, total_depth, depth)
    states.T[:] = linear_gradients(
        bottom_temperature, top_temperature, total_depth, depth
    )


z_top = grid.origin[2] + grid.extent[2]
set_linear_gradients_state(
    simulation.all_states(), z_top - simulation.all_positions()[:, 2]
)

# -------------------------------------------------------------------
# Set the Dirichlet boundary condition values
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(X0)
set_linear_gradients_state(dirichlet, z_top - simulation.vertices()[:, 2])

# -------------------------------------------------------------------
# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

# -------------------------------------------------------------------
# Execute the time loop with solver parameters
simulation.standard_loop(
    initial_timestep=1 * day,
    final_time=200 * year,
    output_period=5 * year,
    newton=newton,
)


# -------------------------------------------------------------------
# Some postprocesses, it allows to visualize with Paraview
simulation.postprocess()
