import numpy as np

import ComPASS
from ComPASS.utils.units import *


# initial reservoir pressure
pL, pR = 2.0 * bar, 1.0 * bar
# initial reservoir temperature - convert Celsius degrees to Kelvin degrees
TL, TR = (
    degC2K(30),
    degC2K(20),
)
# column permeability in m^2 (low permeability -> bigger time steps)
k_matrix = 1e-12
# column porosity
phi_matrix = 0.15
# bulk thermal conductivity in W/m/K
K_matrix = 2.0

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(1, 1, 1),
    extent=(1, 1, 1),
    origin=(0, 0, 0),
)


def dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
    return (x == x.min()) & (x == x.max())


def select_fractures():
    zc = simulation.compute_global_face_centers()[:, 2]
    return (zc < 0.1) | (zc > 0.9)


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    fracture_faces=select_fractures,
    fracture_permeability=k_matrix,
    fracture_porosity=phi_matrix,
    fracture_thermal_conductivity=K_matrix,
)

initial_state = simulation.build_state(simulation.Context.liquid, p=pR, T=TR)
simulation.all_states().set(initial_state)

current_time = simulation.standard_loop(
    initial_timestep=1,
    nitermax=2,
    output_period=1,
)
