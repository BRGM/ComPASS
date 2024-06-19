import numpy as np
import os
import pandas as pd

import ComPASS
import ComPASS.mpi as mpi
import ComPASS.io.mesh as io
import vtkwriters as vtkw
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.utils.various import tensor_coordinates
from ComPASS.newton import Newton
from ComPASS.linalg.petsc_linear_solver import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.utils.salome import SalomeWrapper


ComPASS.set_output_directory_and_logfile(__file__)

Serie = "Orison avec faille"


""" Constant Parameters """
p0: float = 10 * bar
T0: float = degC2K(20.0)
phi_fracture: float = 0.1  # fracture porosity
phi_matrix: float = 0.2  # domain porosity
lambda_th_water: float = 0.65  # Thermal conductivity of liquid phase [kg.m.s^-3.K^-1]
lambda_th_rock: float = 2.0  # Thermal conductivity of solid phase [kg.m.s^-3.K^-1]
heat_capacity_rock: float = 1000.0  # spe_heat_capacity_rock [J.kg-1.K-1]
rho_rock_density: float = 2600.0  # [kg.m-3]
gravity: float = 9.80665
zmax = 0.0
zmin = -1000.0
Tgrad = 0.1  # 100°C/km Geothermal gradient


if mpi.is_on_master_proc:

    # Print params of simulation
    print()
    print("PARAMS of SIMULATION:")
    print("SERIE :", Serie)
    print("NUM :", Num)
    print("KSOCLE [m2]:", K_socle)
    print("KMESO [m2]:", K_meso)
    print("KTERQUAT [m2]:", K_tert_quat)
    print("K_FRACTURE_SOCLE [m2]", K_Fault_socle)
    print("K_FRACTURE_MESO [m2]", K_Fault_meso)
    print("K_FRACTURE_TERQUAT [m2]", K_Fault_tert_quat)
    print("LAMBDASOCLE [kg.m/(K.s^3)]", lambda_th_rock)
    print("LAMBDAMESO [kg.m/(K.s^3)]:", kth_meso)
    print("LAMBDATERQUAT [kg.m/(K.s^3)]:", kth_tert_quat)
    print("FRACTURE_THICKNESS [m]", f_thickness)
    print()

""" ComPASS simulation """
simulation = ComPASS.load_physics("water2ph")

sw = SalomeWrapper(simulation)

""" global parameters """
simulation.set_gravity(gravity)
simulation.set_fracture_thickness(f_thickness)
simulation.set_rock_volumetric_heat_capacity(rho_rock_density * heat_capacity_rock)
bulk_thermal_conductivity = (
    lambda_th_rock * (1 - phi_matrix) + lambda_th_water * phi_matrix
)

if mpi.is_on_master_proc:
    vertices = sw.mesh.vertices_array()
    zmax, zmin = vertices[:, -1].max(), vertices[:, -1].min()

    print("zmax=", zmax, "zmin=", zmin)

""" Functions to define permeabilities: Units and Faults """


def subdomains_permeability(display=False):
    k_cells = np.zeros(sw.info.mesh.nb_cells)
    k_cells[sw.info.unit_socle.cells] = K_socle
    k_cells[sw.info.unit_middle.cells] = K_meso
    k_cells[sw.info.unit_top.cells] = K_tert_quat
    if display:
        MT.to_vtu(sw.mesh, "k_cells.vtu", celldata={"perm": k_cells})
    return k_cells


def subdomains_thconductivity(display=False):
    th_cells = np.full(sw.info.mesh.nb_cells, lambda_th_rock)
    th_cells[sw.info.unit_middle.cells] = kth_meso
    th_cells[sw.info.unit_top.cells] = kth_tert_quat
    if display:
        MT.to_vtu(sw.mesh, "th_cells.vtu", celldata={"th_cond": th_cells})
    return th_cells


def fracture_sets_permeability(display=False):
    k_frac = np.zeros(sw.mesh.nb_faces)
    k_frac[sw.info.Fault_socle.faces] = K_Fault_socle
    k_frac[sw.info.Fault_middle.faces] = K_Fault_meso
    k_frac[sw.info.Fault_top.faces] = K_Fault_tert_quat
    return k_frac[sw.info.Fault.faces]


""" Functions to define Dirichlet boundary condition """

""" Dirichlet Temperature at the top and Bottom if no neumman at the bottom boundary """


def topbottom_dirichlet_nodes():
    where = np.zeros(sw.mesh.nb_vertices, dtype=bool)
    where[sw.info.Top.nodes] = True
    where[sw.info.Bottom.nodes] = True
    return where


""" Dirichlet Pressure at the top """


def top_dirichlet_nodes():
    where = np.zeros(sw.mesh.nb_vertices, dtype=bool)
    where[sw.info.Top.nodes] = True
    return where


simulation.init(
    mesh=sw.mesh,
    cell_permeability=subdomains_permeability,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=subdomains_thconductivity,
    fracture_thermal_conductivity=lambda_th_water,
    fracture_faces=lambda: sw.info.Fault.faces,
    fracture_permeability=fracture_sets_permeability,
    fracture_porosity=phi_fracture,
    set_temperature_dirichlet_nodes=topbottom_dirichlet_nodes,
    set_pressure_dirichlet_nodes=top_dirichlet_nodes,
    # set_dirichlet_nodes=set_dirichlet_nodes,
    set_global_flags=sw.flags_setter,
)

sw.rebuild_locally()


""" Set global initial conditions: P0, TO = f(z)"""
""" Function: linear gradient useful for Pre and Temp"""


def lininterp(z, top, gradient):
    return top + (gradient) * (z)


def set_initial_states():
    gravity = simulation.get_gravity()

    def set_states(states, z):
        states.context[:] = simulation.Context.liquid
        states.p[:] = lininterp(
            zmax - z,
            p0,
            gravity * simulation.liquid_molar_density(p0, T0),
        )
        states.T[:] = lininterp(zmax - z, T0, Tgrad)
        states.S[:] = [0, 1]
        states.C[:] = 1

    z = simulation.vertices()[:, 2]
    set_states(simulation.node_states(), z)
    set_states(simulation.dirichlet_node_states(), z)
    z = simulation.compute_cell_centers()[:, 2]
    set_states(simulation.cell_states(), z)
    z = simulation.compute_fracture_centers()[:, 2]
    set_states(simulation.fracture_states(), z)


set_initial_states()


""" Function to export all initial values as vtufile"""


def export_initial_states():
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    rho = simulation.liquid_molar_density
    petrophysics = simulation.petrophysics()

    pointdata = {
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
        "initial pressure": node_states.p,
        "initial temperature": K2degC(node_states.T),
        # "initial context": node_states.C,
        "initial saturation": node_states.S[:, 0],
        "liquid density": rho(node_states.p, node_states.T),
        # "Psat for T reservoir": simulation.Psat(node_states.T),
        # "Tsat for p reservoir": K2degC(simulation.Tsat(node_states.p)),
    }
    celldata = {
        "initial pressure": cell_states.p,
        "initial saturation": cell_states.S,
        "initial temperature": K2degC(cell_states.T),
        "phi": petrophysics.cell_porosity,
    }
    celldata.update(
        tensor_coordinates(petrophysics.cell_permeability, "k", diagonal_only=True)
    )
    # celldata.update(
    # tensor_coordinates(petrophysics.cell_thermal_conductivity, "K", diagonal_only=True)
    # )
    io.write_mesh(simulation, "initial_states", pointdata=pointdata, celldata=celldata)


export_initial_states()


""" Solver parameters """
lsolver = linear_solver(simulation, direct=False)
newton = Newton(simulation, 1e-5, 8, lsolver)
tsmger = TimeStepManager(
    initial_timestep=1 * day,  # 10
    increase_factor=2,
    decrease_factor=1 / 1.5,
    karma_threshold=10,
    karma_reset=False,
)


# loading wells data___________________________________________________________________________________________________
wells_file_path = "/home/haouchine/orison_avec_faille/Observables_sampling_noise.txt"
wells_data = np.loadtxt(
    wells_file_path, skiprows=1, delimiter=" ", dtype=str
)  # X, Y, Z, T


""" Simulation """
final_time = 1e3 * year
simulation.standard_loop(
    final_time=final_time,
    time_step_manager=tsmger,
    output_period=0.1 * final_time,
    # newton = newton,
    output_after_loop=True,
    output_before_start=True,
)

# simulation.postprocess()    # =========================================================================================


""" Function to find the index of the node for which the coordinates X, Y, Z are the closest to the coordinates of the observed data """


def closest_index(xobs, yobs, zobs, xnode, ynode, znode):
    distances = np.sqrt((xnode - xobs) ** 2 + (ynode - yobs) ** 2 + (znode - zobs) ** 2)
    closest_idx = np.argmin(distances)
    return closest_idx


# Access nodes coordinates
node_coordinates = simulation.vertices()

# Access nodes and cells state
node_states = simulation.node_states()

# Access predicted temperatures at nodes and cells
predicted_node_temperatures = node_states.T

# Decompose nodes coordinates
node_x_coordinates = node_coordinates[:, 0]
node_y_coordinates = node_coordinates[:, 1]
node_z_coordinates = node_coordinates[:, 2]

# Access well coordinate and temperature
x_obs = np.array(wells_data[:, 0], dtype=float)
y_obs = np.array(wells_data[:, 1], dtype=float)
z_obs = np.array(wells_data[:, 2], dtype=float)
T_obs = np.array(wells_data[:, 3], dtype=float)

# Initialize arrays
T_pred = []
predicted_coordinates = []
T_misfits = []

# Find the predicted temperature corresponding to the wells coordinates
for x, y, z in zip(x_obs, y_obs, z_obs):
    # Finding the index of the closest node
    closest_idx = closest_index(
        x, y, z, node_x_coordinates, node_y_coordinates, node_z_coordinates
    )

    # Extracting the predicted temperature and pressure of the node
    predicted_temperature = predicted_node_temperatures[closest_idx]
    predicted_coordinates.append([x, y, z])
    T_pred.append(predicted_temperature)

# Calculate the misfit for each pair of observed and predicted temperatures
for T_observed, T_predicted in zip(T_obs, T_pred):

    T_misfit = T_observed - T_predicted
    T_misfit = np.absolute(T_misfit)
    T_misfits.append(T_misfit)

T_misfits_array = np.array(T_misfits)

# Calculating Mean Squared Error (MSE)
Tmse = np.mean((T_misfits_array) ** 2)

# Calculating Mean Absolute Error (MAE)
Tmae = np.mean(np.abs(T_misfits_array))

# Calculating Correlation Coefficient
Tcorrelation_coef = np.corrcoef(T_obs, T_pred)[0, 1]

# Save data to a csv
Dataframe = pd.DataFrame(
    {
        "X_obs": x_obs,
        "Y_obs": y_obs,
        "Z_obs": z_obs,
        "Observed_Temperature": T_obs,
        "X_pred": [coord[0] for coord in predicted_coordinates],
        "Y_pred": [coord[1] for coord in predicted_coordinates],
        "Z_pred": [coord[2] for coord in predicted_coordinates],
        "Predicted_Temperature": T_pred,
        "T_Misfits": T_misfits,
    }
)

Dataframe["TMSE"] = np.nan
Dataframe.loc[0, "TMSE"] = Tmse

Dataframe["TMAE"] = np.nan
Dataframe.loc[0, "TMAE"] = Tmae

Dataframe["TR"] = np.nan
Dataframe.loc[0, "TR"] = Tcorrelation_coef

# Vérifier si le dossier existe, sinon le créer
if not os.path.exists("output_folder"):
    os.makedirs("output_folder")

# Enregistrer le DataFrame dans un fichier CSV
Dataframe.to_csv(os.path.join("output_folder", "data_csv"), index=False, sep=";")
