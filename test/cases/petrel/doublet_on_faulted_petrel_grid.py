from pathlib import Path

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.utils.petrel import PetrelGrid
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager

filename = Path("Simple20x20x5_Fault.grdecl")

# fmt: off
pres = 20. * MPa            # initial reservoir pressure
Tres = degC2K( 70. )        # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 30. )  # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour      # production flowrate
rw = 0.1                    # well radius in meters
omega_reservoir = 0.15      # reservoir porosity
k_reservoir = 1E-13         # reservoir permeability in m^2
K_reservoir = 2             # bulk thermal conductivity in W/m/K
omega_fracture = 0.5        # reservoir porosity
k_fracture = 1E-11          # reservoir permeability in m^2
K_fracture = 2              # bulk thermal conductivity in W/m/K
g = 9.81                    # gravity in m/s^-2
producer_id, injector_id = 0, 1
# fmt: on

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(g)

grid = PetrelGrid(filename)
# scale mesh to realistic dimension
if grid.mesh is not None:
    grid.mesh.vertices[:, :2] *= 2000  # dx=100
    grid.mesh.vertices[:, 2] *= 200


def select_dirichlet():
    xyz = simulation.global_vertices()
    x, y, z = xyz.T
    # sort by x, then by y, then by z
    order = np.lexsort((z, y, x))
    where = np.zeros(xyz.shape[0], dtype="b")
    where[order[-1]] = True
    return where


def make_wells():
    # here well segments are hard coded
    xyz = simulation.global_vertices()
    injector = simulation.create_single_branch_well((946, 947), well_radius=rw)
    injector.id = injector_id
    injector.operate_on_flowrate = Qm, pres + 100.0 * MPa
    injector.inject(Tinjection)
    producer = simulation.create_single_branch_well((1828, 1829), well_radius=rw)
    producer.id = producer_id
    producer.operate_on_flowrate = Qm, 1.0 * bar
    producer.produce()
    return (producer, injector)


simulation.init(
    mesh=grid.mesh,
    set_dirichlet_nodes=select_dirichlet,
    wells=make_wells,
    fracture_faces=grid.fracture_faces,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    fracture_porosity=omega_fracture,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=K_fracture,
    dump_mesh_before_distribution=True,
)

# You may want to free memory on master proc on the simulation is initialized
del grid

zref = 1500
pref = 1 * bar  # zref = 3000m
rhoref = 1000.0  # kg / m3

X0 = simulation.build_state(simulation.Context.liquid, p=rhoref * g * zref, T=Tres)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

# -- Step 1: permanent state -------------------------------------------------

simulation.close_well(injector_id)
simulation.close_well(producer_id)

# search pressure equilibrium
simulation.standard_loop(
    initial_timestep=30 * day, final_time=10 * year, no_output=True,
)

# -- Step 2: doublet production -----------------------------------------------

simulation.open_well(injector_id)
simulation.open_well(producer_id)

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)
tsmger = TimeStepManager(initial_timestep=hour, increase_factor=2)

simulation.standard_loop(
    initial_time=0,
    final_time=10 * year,
    output_period=year,
    newton=newton,
    time_step_manager=tsmger,
)

simulation.postprocess()
