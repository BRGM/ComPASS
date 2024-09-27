# Documentation : https://compass.gitlab.io/v4/doc/

#%% import the ComPASS library and some utilities
import ComPASS
from ComPASS.utils.units import *  # contains MPa, degC2K, km, year...

#  |------------------------------------|
#  |                                    |
#  |                                    |
#  |             Reservoir              |
#  |                                    |
#  |                                    |
#  |------------------------------------|


#%% Load the water2ph physics : water component in liquid and/or gas phase
simulation = ComPASS.load_physics("water2ph")


#%% Create a Cartesian grid with cubic cells
# Geometry
total_depth = 240.0  # m
nx = nz = 31
ny = 1
d = total_depth / nz
grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(nx * d, ny * d, nz * d),  # creates cubic cells
    origin=(0.0, 0.0, -total_depth),
)


#%% Initialize the regionalized values and distribute the domain
# Reservoir petrophysics
k_reservoir = 1e-12  # reservoir permeability in m^2
omega = 0.15  # reservoir porosity
K_reservoir = 2  # reservoir bulk thermal conductivity in W/m/K

simulation.init(
    mesh=grid,
    cell_porosity=omega,
    cell_thermal_conductivity=K_reservoir,
    cell_permeability=k_reservoir,
)


#%% Initialize the domain
# Initialize the domain with liquid phase
# at 10 bar and 5 Celsius degrees.
temperature = degC2K(5.0)
pressure = 10 * bar
# First construct the state (at thermodynamic equilibrium)
X0 = simulation.build_state(
    simulation.Context.liquid,
    p=pressure,
    T=temperature,
)
# then apply this state everywhere to init the domain (constant values)
simulation.all_states().set(X0)


#%% Execute the time loop
final_time = 0.01 * year
simulation.standard_loop(
    initial_timestep=2 * day,
    final_time=final_time,
    output_period=0.5 * final_time,
)


#%% Some postprocesses, it allows to visualize with Paraview
simulation.postprocess(time_unit="day")
