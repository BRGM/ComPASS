# Documentation : https://compass.gitlab.io/v4/doc/

import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timeloops import TimeStepManager
from ComPASS.utils.salome import SalomeWrapper


# Reservoir petrophysics
k_matrix = 1e-13  # reservoir permeability in m^2
omega_matrix = 0.15  # reservoir porosity
K = 2.0  # bulk thermal conductivity in W/m/K (same in fracture)
# Fracture
k_fracture = 1e-10  # fracture permeability in m^2
omega_fracture = 0.5  # fracture porosity


# -------------------------------------------------------------------
# Load the water2ph physics : it contains the water component
# which can be in liquid and/or gas phase
simulation = ComPASS.load_physics("diphasic")


# -------------------------------------------------------------------
# Import a very coarse mesh created with Salome
sw = SalomeWrapper(
    simulation,
    nodes_file="NODES.txt",
    tets_file="TETRAS.txt",
    groups_module="GROUPS",
)


# -------------------------------------------------------------------
# You can write the mesh importation (matrix and faces blocks)
# to visualize it (only master proc know the mesh)
# if mpi.is_on_master_proc:
#     sw.info.to_vtu_block("salome-block")  # creates salome-block.vtu
#     sw.info.faces_to_multiblock("salome-faults")  # creates salome-faults.vtm


# -------------------------------------------------------------------
# Fracture factory
# set the fracture thickness (global variable)
fracture_thickness = 1  # m
simulation.set_fracture_thickness(fracture_thickness)


# -------------------------------------------------------------------
# Initialize the regionalized values and distribute the domain
simulation.init(
    salome_wrapper=sw,
    cell_permeability=k_matrix,
    cell_porosity=omega_matrix,
    cell_thermal_conductivity=K,
    fracture_faces=lambda: sw.info.fault.faces,
    fracture_permeability=k_fracture,
    fracture_porosity=omega_fracture,
    fracture_thermal_conductivity=K,
)


# -------------------------------------------------------------------
# Initialize the domain with liquid phase
top_pressure = 1.0 * bar  # top reservoir pressure
top_temperature = degC2K(20.0)  # top reservoir temperature

# First construct and apply everywhere the state (at thermodynamic equilibrium)
Xtop = simulation.build_state(
    simulation.Context.liquid,
    p=top_pressure,
    T=top_temperature,
)
simulation.all_states().set(Xtop)

# -------------------------------------------------------------------
# Apply Dirichlet BC at the top nodes
# the Dirichlet node states are copied from the actual states values
simulation.reset_dirichlet_nodes(sw.info.top.nodes)


# -------------------------------------------------------------------
# Define time step parameters
tsm = TimeStepManager(
    initial_timestep=30 * day,
    increase_factor=1.5,
    decrease_factor=0.5,
)

# define the time loop solver options
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 10, lsolver)


# -------------------------------------------------------------------
# Execute the time loop
run_loop = lambda initial_time=0, final_time=1000 * year, nb_output=10, nitermax=None: simulation.standard_loop(
    initial_time=initial_time,
    final_time=final_time,
    time_step_manager=tsm,
    nb_output=nb_output,
    nitermax=nitermax,
    newton=newton,
)

# execute first time loop without Neumann flux (reservoir init)
current_time = run_loop(
    final_time=1 * year,
    nb_output=4,
)
