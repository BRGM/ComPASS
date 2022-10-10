import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *

# from ComPASS.newton import Newton
# from ComPASS.linalg.factory import linear_solver
import MeshTools as MT

# import ComPASS.mpi as mpi
# from ComPASS.mpi import MPI # underlying mpi4py.MPI

from angular_sector import extruded_sector
from my_kr import kr_functions

# from generate_mesh import R, theta, ds, H, nb_layers

R = 1000  # radius of the angular sector
theta = np.pi / 12  # angle of the angular sector (in radians)
rw = 0.1  # well radius
ds = 1, 100  # target edge size around the origin and at the boundary
H = 1.1  # reservoir thickness
nb_layers = 11  # number of horizontal layers
omega = 0.2  # reservoir porosity
k_fracture = 3e-13  # fracture permeability in m^2
k_overburden = 3e-17  # overburden permeability in m^2 (no horizontal permeability)
K = 0  # bulk thermal conductivity in W/m/K
qw = 0.028  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
# pres = 90 * bar  # reservoir pressure
Tres = degC2K(234)

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("water2ph")
simulation.set_gravity(0)
simulation.set_rock_volumetric_heat_capacity(2.65e6)  # SI units J/m^3/°C

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))
# mesh_data = np.load("angular_mesh.npz")
# mesh = MT.WedgeMesh.make(mesh_data["vertices"], mesh_data["cellnodes"])

epsilon = 0.01 * ds[0]  # tolerance value to select nodes (boundary conditions...)


def make_well():
    well = simulation.create_vertical_well((0, 0), rw, zmin=0, zmax=0.1)
    well.operate_on_flowrate = (theta / (2 * np.pi)) * qw, 1 * bar
    well.produce()
    return [well]


def set_permeability():
    centers = simulation.compute_global_cell_centers()
    z = centers[:, 2]
    fracture_cells = z < 0.1
    overburden_cells = z > 0.1
    nb_cells = centers.shape[0]
    k = np.zeros((nb_cells, 3, 3), dtype="d")
    k[fracture_cells] = k_fracture * np.eye(3)  # isotropic
    k[overburden_cells, 2, 2] = k_overburden  # only z component
    return k


simulation.init(
    mesh=mesh,
    wells=make_well,
    cell_porosity=omega,
    cell_permeability=set_permeability,
    cell_thermal_conductivity=K,
)

# Set the kr functions after initialization
simulation.set_kr_functions(kr_functions)

# Set initial values
X0 = simulation.build_state(simulation.Context.diphasic, p=30.5 * bar, Sg=1)
simulation.all_states().set(X0)
# Set boundary conditions
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
simulation.reset_dirichlet_nodes((x**2 + y**2) > R**2 - epsilon)

simulation.standard_loop(
    fixed_timestep=1,
    nitermax=104,
    output_period=10,
)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess()
