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
sw = SalomeWrapper(
    simulation,
    nodes_file="NODES.txt",
    tets_file="TETRAS.txt",
    groups_module="GROUPS",
)


# -------------------------------------------------------------------
# If necessary you can write the mesh importation (matrix and faces blocks)
# to visualize it
# sw.info.to_vtu_block("salome-block")
# sw.info.faces_to_multiblock("salome-faults")


# -------------------------------------------------------------------
# Create the two wells using the list of there ordered nodes
injector_id = 0  # you chose the id you want
producer_id = 1  # you chose the id you want
Qm = 200.0 * ton / hour  # flowrate
well_radius = 0.115  # well radius
injection_temperature = degC2K(40.0)


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
    fracture_faces=lambda: sw.info.fault.faces,
    fracture_permeability=k_fracture,
    fracture_porosity=omega_fracture,
    fracture_thermal_conductivity=K,
    set_global_flags=sw.flags_setter,
)


# -------------------------------------------------------------------
# Initialize the domain with the last states saved in snapshot_directory

# -------------------------------------------------------------------
# Identify the Dirichlet nodes (all vertical walls)
sw.rebuild_locally()
dirichlet_nodes = np.zeros(sw.mesh.nb_vertices, dtype=bool)
dirichlet_nodes[sw.info.east.nodes] = True
dirichlet_nodes[sw.info.north.nodes] = True
dirichlet_nodes[sw.info.west.nodes] = True
dirichlet_nodes[sw.info.south.nodes] = True

# the Dirichlet node states are copied from the actual states values
simulation.reset_dirichlet_nodes(dirichlet_nodes)


# -------------------------------------------------------------------
# Define well connections between source well (producer) and target well (injector)


# Update imposed_flowrate and injection_temperature
# of the target well (injector one) knowing the information
# of the connected source well (producer one).
# def chain_wells(tick):
#     for source, target in doublet:
#         injector_data = simulation.get_well_data(target)
#         producer_wellhead = simulation.well_connections[source]
#         molar_flowrate = producer_wellhead.molar_flowrate[
#             0
#         ]  # "0" because there is a single component
#         assert (
#             molar_flowrate >= 0
#         ), f"source well {source} with flowrate {molar_flowrate} should be a producer"
#         injector_data.imposed_flowrate = ???
#         injector_data.injection_temperature = ???


# -------------------------------------------------------------------
# Define time step parameters
tsm = TimeStepManager(
    initial_timestep=2 * hour,
    increase_factor=1.2,
    decrease_factor=0.5,
)

# -------------------------------------------------------------------
# Execute the time loop
final_time = 10 * year
output_period = 0.5 * year
simulation.standard_loop(
    final_time=final_time,
    time_step_manager=tsm,
    output_period=output_period,
)


# -------------------------------------------------------------------
# Some postprocesses, it allows to visualize with Paraview
simulation.postprocess()
