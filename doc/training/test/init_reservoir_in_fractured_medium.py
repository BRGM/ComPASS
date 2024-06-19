# Documentation : https://charms.gitlabpages.inria.fr/ComPASS/
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import TimeStepManager
from ComPASS.utils.salome import SalomeWrapper
import ComPASS.mpi as mpi


# Reservoir petrophysics
k_matrix = 1e-13  # reservoir permeability in m^2
omega_matrix = 0.15  # reservoir porosity
K = 2.0  # bulk thermal conductivity in W/m/K (same in fracture)
# Fracture
k_fracture = 1e-10  # fracture permeability in m^2
omega_fracture = 0.3  # fracture porosity


# -------------------------------------------------------------------
# Set output informations
ComPASS.set_output_directory_and_logfile(__file__)
# Load the water2ph physics : it contains the water component
# which can be in liquid and/or gas phase
simulation = ComPASS.load_physics("water2ph")


# -------------------------------------------------------------------
# Import a very coarse mesh created with Salome


# -------------------------------------------------------------------
# If necessary you can write the mesh importation (matrix and faces blocks)
# to visualize it
# sw.info.to_vtu_block("salome-block")
# sw.info.faces_to_multiblock("salome-faults")

# -------------------------------------------------------------------
# Fracture factory
# set the fracture thickness (global variable)
fracture_thickness = 0.3  # m
simulation.set_fracture_thickness(fracture_thickness)
# the list of fracture faces is in sw.info.fault.faces

# -------------------------------------------------------------------
# Initialize the regionalized values and distribute the domain
simulation.init(
    mesh=sw.mesh,
    cell_permeability=k_matrix,
    cell_porosity=omega_matrix,
    cell_thermal_conductivity=K,
    fracture_permeability=k_fracture,
    fracture_porosity=omega_fracture,
    fracture_thermal_conductivity=K,
    set_global_flags=sw.flags_setter,
)


# -------------------------------------------------------------------
# Initialize the domain with liquid phase
top_pressure = 10 * bar  # top reservoir pressure
top_temperature = degC2K(70.0)  # top reservoir temperature
bottom_temperature = degC2K(150.0)  # bottom reservoir temperature

# First construct and apply everywhere the state (at thermodynamic equilibrium)
X0 = simulation.build_state(
    simulation.Context.liquid,
    p=top_pressure,
    T=top_temperature,
)
simulation.all_states().set(X0)

# Modify the initial values to apply a linear gradient
# in temperature
def linear_gradients(bottom_value, top_value, domain_depth, depth):
    return top_value + ((bottom_value - top_value) / domain_depth) * depth


z = simulation.all_positions()[:, 2]
z_top = np.max(z)
total_depth = z_top - np.min(z)
assert total_depth > 0.0
simulation.all_states().T[:] = linear_gradients(
    bottom_temperature, top_temperature, total_depth, z_top - z
)


# -------------------------------------------------------------------
# Identify the Dirichlet nodes (at the top)
sw.rebuild_locally()
dirichlet_nodes = np.zeros(sw.mesh.nb_vertices, dtype=bool)
dirichlet_nodes[sw.info.top.nodes] = True

# the Dirichlet node states are copied from the actual states values
simulation.reset_dirichlet_nodes(dirichlet_nodes)


# -------------------------------------------------------------------
# Define time step parameters
tsm = TimeStepManager(
    initial_timestep=2 * hour,
    increase_factor=1.2,
    decrease_factor=0.5,
)

simulation.standard_loop(
    final_time=10 * year,
    output_period=1 * year,
    time_step_manager=tsm,
)


# -------------------------------------------------------------------
# Some postprocesses, it allows to visualize with Paraview
simulation.postprocess()
