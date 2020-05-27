import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.newton import Newton
from ComPASS.legacy_petsc import LegacyLinearSolver

# import ComPASS.mpi as mpi
# from ComPASS.mpi import MPI # underlying mpi4py.MPI

from angular_sector import extruded_sector


R = 1000  # radius of the angular sector
theta = np.pi / 6  # angle of the angular sector (in radians)
ds = 10, 100  # target edge size around the origin and at the boundary
H = 10  # reservoir thickness
nb_layers = 6  # number of horizontal layers
rw = 0.1  # well radius
omega = 0.15  # reservoir porosity
k = 1e-12  # reservoir permeability in m^2
K = 2  # bulk thermal conductivity in W/m/K
qw = 100.0  # production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 10 * MPa  # reservoir pressure
Tres = degC2K(60)  # reservoir temperature

#%% -------------------------------------------------------------------

ComPASS.set_output_directory_and_logfile(f"angular_sector_{nb_layers:02d}_layers")
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))
epsilon = 0.01 * ds[0]
from MeshTools import to_vtu

to_vtu(mesh, "sector")


def make_well():
    well = simulation.create_vertical_well((0, 0), rw)
    well.operate_on_flowrate = (theta / (2 * np.pi)) * qw, 1 * bar
    well.produce()
    return [well]


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    return np.linalg.norm(vertices[:, :2], axis=1) > R - epsilon


def build_permeability():
    centers = simulation.compute_global_cell_centers()
    result = np.zeros(centers.shape[0], dtype=np.double)
    result[:] = 0.01 * k
    zc = centers[:, 2]
    result[(zc >= H / 3) & (zc <= 2 * H / 3)] = k
    # result = np.zeros((centers.shape[0], 3, 3), dtype=np.double)
    # result[..., 0, 0] = kres
    # result[..., 1, 1] = kres
    return result


simulation.init(
    mesh=mesh,
    set_dirichlet_nodes=select_dirichlet_nodes,
    wells=make_well,
    cell_porosity=omega,
    cell_permeability=k,
    cell_thermal_conductivity=K,
)

X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
dirichlet = simulation.dirichlet_node_states().set(X0)

import ComPASS.mpi as mpi

assert mpi.communicator().size == 1, "Sequential runs only!"

times = []
timesteps = []
wellheads = []

axis_nodes = np.nonzero(np.fabs(simulation.vertices()[:, 1]) < 1e-5)[0]
node_states = simulation.node_states()
axis_vertices = np.array(simulation.vertices()[axis_nodes], copy=True)
axis_pressure = []

simulation.add_well_connections([(0, 0)])


def collect_timestep(tick):
    times.append(tick.time)
    timesteps.append(tick.latest_timestep)
    axis_pressure.append(np.array(node_states.p[axis_nodes], copy=True))
    wellhead = simulation.well_connections[0]
    wellheads.append(wellhead)
    print(
        f"Well head state: q={wellhead.mass_flowrate:.2f} p={wellhead.pressure / bar:.2f}, T={K2degC(wellhead.temperature):.2f}"
    )


# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = LegacyLinearSolver(activate_direct_solver=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    newton=newton,
    initial_timestep=1,
    final_time=year,
    output_period=year,
    iteration_callbacks=[collect_timestep],
)


def simulation_hash_code(**kwargs):
    return f"{hash(frozenset(list(kwargs.items()))) % ((sys.maxsize + 1) * 2):x}"


def dump_summary(params, objects):
    hash_code = simulation_hash_code(**params)
    with open("analysis.yaml", "a") as f:
        print(yaml.dump([{"hash": hash_code, "params": params}]), file=f)
    with open(f"results-{hash_code}", "wb") as f:
        pickle.dump(objects, f)


def pack(**kwargs):
    return kwargs


def pack_as_arrays(**kwargs):
    return {name: np.array(value) for name, value in kwargs.items()}


dump_summary(
    params=pack(qw=qw, R=R, nb_layers=nb_layers, k=k, theta=theta, ds=ds),
    objects=pack_as_arrays(
        times=times,
        timesteps=timesteps,
        axis_vertices=axis_vertices,
        axis_pressure=axis_pressure,
        pwh=[wh.pressure for wh in wellheads],
        qwh=[wh.mass_flowrate for wh in wellheads],
    ),
)
