import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *

# from ComPASS.newton import Newton
# from ComPASS.linalg.factory import linear_solver
from ComPASS.timeloops import TimeStepManager
import MeshTools as MT

# import ComPASS.mpi as mpi
# from ComPASS.mpi import MPI # underlying mpi4py.MPI

from ComPASS.utils.angular_sector import extruded_sector
from my_kr import kr_functions

# from generate_mesh import R, theta, ds, H, nb_layers

R = 1000  # radius of the angular sector
theta = np.pi / 3  # angle of the angular sector (in radians)
rw = 0.1  # well radius
ds = 1, 100  # target edge size around the origin and at the boundary
H = 100  # reservoir thickness
nb_layers = 10  # number of horizontal layers
omega = 0.2  # reservoir porosity
k = 1e-14  # reservoir permeability in m^2
K = 0  # bulk thermal conductivity in W/m/K
qw = 14.0  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 90 * bar  # reservoir pressure
Tres = degC2K(300)
gravity = 9.81

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("water2ph")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(2.65e6)  # SI units J/m^3/°C

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))
# mesh_data = np.load("angular_mesh.npz")
# mesh = MT.WedgeMesh.make(mesh_data["vertices"], mesh_data["cellnodes"])

epsilon = 0.01 * ds[0]  # tolerance value to select nodes (boundary conditions...)


def make_well():
    well = simulation.create_vertical_well((0, 0), rw)
    well.operate_on_flowrate = (theta / (2 * np.pi)) * qw, 1 * bar
    well.produce()
    return [well]


simulation.init(
    mesh=mesh,
    wells=make_well,
    cell_porosity=omega,
    cell_permeability=k,
    cell_thermal_conductivity=K,
    well_model="two_phases",
)

# Set the kr functions after initialization
simulation.set_kr_functions(kr_functions)

# Set initial values
X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
# Set boundary conditions
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
simulation.reset_dirichlet_nodes((x**2 + y**2) > R**2 - epsilon)

# tsmger = TimeStepManager(
#     initial_timestep=1e-4 * day,
#     maximum_timestep=20,
#     increase_factor=1.1,
#     decrease_factor=0.5,
# )

# lsolver = linear_solver(simulation, direct=True)
# newton = Newton(simulation, 1e-5, 8, lsolver)

# simulation.standard_loop(
#     final_time=day, output_period=hour,
#     newton=newton,
#  time_step_manager=tsmger,
#  )

simulation.standard_loop(
    initial_timestep=1e-4 * day,
    final_time=day,
    output_period=0.5 * hour,
    #   time_step_manager=tsmger, well_pressure_offset=1
)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess()
