import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.timeloops import TimeStepManager, Event

from ComPASS.utils.angular_sector import extruded_sector

# The topmost fracture is finite (cf. fracture_radius)
# The bottommost fracture is infinite and ensures recharge.

R = 1000  # radius of the angular sector
fracture_radius = 200  # radius of the fracture in meters
fracture_thickness = 0.1  # fracture thickness in meters
theta = np.pi / 6  # angle of the angular sector (in radians)
ds = 10, 100  # target edge size around the origin and at the boundary
H = 10  # reservoir thickness
nb_layers = 6  # number of horizontal layers
rw = 0.1  # well radius
k_fracture = 1e-12  # fracture permeability in m^2
k_matrix = 1e-20  # matrix permeability in m^2
omega_fracture = 0.5  # fracture porosity
omega_matrix = 0.15  # reservoir porosity
K = 2  # bulk thermal conductivity in W/m/K
pres = 10 * MPa  # reservoir pressure
Tres = degC2K(60)  # reservoir temperature
gravity = 9.81

Qinjection = 10.0  # injection flow rate (for the whole reservoir i.e. theta = 2pi)
injection_duration = 10 * hour
Tinjection = degC2K(30)
well_id = 0

#%% -------------------------------------------------------------------

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(gravity)
simulation.set_fracture_thickness(fracture_thickness)

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))
epsilon = 0.01 * ds[0]


def make_fractures():
    face_centers = simulation.compute_global_face_centers()
    dz = H / nb_layers
    # select horizontal fault axis in the middle of the simulation domain
    x, y, z = [face_centers[:, j] for j in range(3)]
    top_fracture = (x**2 + y**2 < fracture_radius**2) & (
        np.abs(z - (2 / 3) * H) < 0.25 * dz
    )
    bottom_fracture = np.abs(z - (1 / 3) * H) < 0.25 * dz
    return bottom_fracture | top_fracture


def make_injector():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = well_id
    well.operate_on_flowrate = (theta / (2 * np.pi)) * Qinjection, np.inf
    well.inject(Tinjection)
    return [well]


simulation.init(
    mesh=mesh,
    wells=make_injector,
    cell_porosity=omega_matrix,
    cell_permeability=k_matrix,
    cell_thermal_conductivity=K,
    fracture_faces=make_fractures,
    fracture_porosity=omega_fracture,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=K,
)

## ------------------------------
## First step - equilibrium state
## ------------------------------

# top nodes are set to dirichlet conditions
simulation.reset_dirichlet_nodes(lambda xyz: xyz[:, 2] > H - epsilon)
X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(X0)

simulation.close_well(well_id)

tsmger = TimeStepManager(
    initial_timestep=day,
    maximum_timestep=0.3 * year,
    increase_factor=2.0,
    decrease_factor=0.1,
)

simulation.standard_loop(
    final_time=10 * year,  # more than enough to reach pressure equilibrium
    no_output=True,
    time_step_manager=tsmger,
)

## ------------------------
## Second step - production
## ------------------------

# Dirichlet nodes are switch to boundary nodes
# they will be locked to their equilibrium values
simulation.reset_dirichlet_nodes(
    lambda xyz: np.linalg.norm(xyz[:, :2], axis=1) > R - epsilon
)

simulation.open_well(well_id)

tsmger = TimeStepManager(
    initial_timestep=1.0,
    maximum_timestep=hour,
    increase_factor=1.5,
    decrease_factor=0.5,
)


def close_well(tick):
    simulation.close_well(well_id)
    # we force a small time steps at closing time
    tsmger.current_step = 1.0


# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    initial_time=0,
    newton=newton,
    final_time=3 * injection_duration,
    output_period=hour,
    specific_outputs=[n * minute for n in [1, 2, 3, 5, 10, 15, 20, 25, 30, 45]]
    + [injection_duration + n * minute for n in [1, 2, 3, 5, 10, 15, 20, 25, 30, 45]],
    time_step_manager=tsmger,
    events=[Event(injection_duration, [close_well])],
)
