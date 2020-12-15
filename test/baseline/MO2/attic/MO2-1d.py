import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver

# import ComPASS.mpi as mpi
# from ComPASS.mpi import MPI # underlying mpi4py.MPI

from angular_sector import extruded_sector
from my_kr import kr_functions

R = 1000  # radius of the angular sector
theta = np.pi / 12  # angle of the angular sector (in radians)
rw = 0.1  # well radius
ds = 10, 100  # target edge size around the origin and at the boundary
H = 100  # reservoir thickness
nb_layers = 1  # number of horizontal layers
omega = 0.2  # reservoir porosity
k = 1e-14  # reservoir permeability in m^2
K = 0  # bulk thermal conductivity in W/m/K
qw = 14.0  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 90 * bar  # reservoir pressure
Tres = degC2K(260)  # reservoir temperature

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)
simulation.set_rock_volumetric_heat_capacity(2.65e6)  # SI units J/m^3/°C
simulation.set_kr_functions(kr_functions)

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))
epsilon = 0.01 * ds[0]  # tolerance value to select nodes (boundary conditions...)
from MeshTools import to_vtu

to_vtu(mesh, "sector")


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
)

# Set initial values
X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
# Set boundary conditions
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
simulation.reset_dirichlet_nodes((x ** 2 + y ** 2) > R ** 2 - epsilon)

simulation.standard_loop(
    initial_timestep=1e-5 * day, final_time=day, output_period=hour,
)
