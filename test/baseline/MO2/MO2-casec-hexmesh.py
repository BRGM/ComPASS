import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.utils.angular_sector import extruded_sector
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

from ComPASS.timeloops import TimeStepManager

import hexmesh
from my_kr import kr_functions

R = 1000  # radius of the angular sector
theta = np.pi / 180  # angle of the angular sector (in radians) must be small
rw = 0.1  # well radius
r0 = 0.5 * rw
nr = 200  # number of radius discretization points
H = 100  # reservoir thickness
omega = 0.2  # reservoir porosity
k = 1e-14  # reservoir permeability in m^2
K = 0  # bulk thermal conductivity in W/m/K
qw = 14  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 90 * bar  # reservoir pressure
Tres = degC2K(300)
gravity = 0
output_period = 30 * minute
specific_outputs = [1, 2, 5, 10, 30, 60, 120, 180, 240, 300, 600, 900, 1200, 1500]
final_time = 10 * day
maximum_timestep = 30  # s

# -----------------------------------------------------------------------------

epsilon = 0.01 * r0  # tolerance value to select nodes (boundary conditions...)

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(2.65e6)  # SI units J/m^3/°C

a = hexmesh.geometric_progression(nr - 1, r0, R - rw)
r = [0, r0]
for _ in range(nr - 1):
    r.append(r[-1] * a)

mesh = hexmesh.make_mesh(rw + np.cumsum(r), theta, rw)


def select_well_face(centers):
    x, y, z = [centers[:, j] for j in range(3)]
    return (x**2 + y**2) <= rw**2 + epsilon


def select_dirichlet(points):
    x, y, z = [points[:, j] for j in range(3)]
    return (x**2 + y**2) >= R**2 - epsilon


simulation.init(
    mesh=mesh,
    cell_porosity=omega,
    cell_permeability=k,
    cell_thermal_conductivity=K,
)

simulation.set_kr_functions(kr_functions)

# initial values
X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)

# boundary conditions
vertices = simulation.vertices()
# assertion works for sequential script
# assert np.count_nonzero(select_dirichlet(vertices)) == 4
simulation.reset_dirichlet_nodes(select_dirichlet(vertices))

# Neumann at the well face
face_centers = simulation.face_centers()
well_face = select_well_face(face_centers)
mass_flux = qw / (2 * np.pi * rw * H)
Neumann = ComPASS.NeumannBC(-mass_flux, compute_heat_flux=True)
simulation.set_Neumann_faces(well_face, Neumann)

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 20, lsolver)

all_states = simulation.all_states()


def check_physics(tick):
    print(all_states.p.min(), all_states.T.min())
    assert np.all(all_states.p > 1e5)
    assert np.all(all_states.T > 280)
    pass


tsmger = TimeStepManager(
    initial_timestep=1,
    increase_factor=1.2,
    decrease_factor=0.8,
    maximum_timestep=maximum_timestep,
)

simulation.standard_loop(
    final_time=final_time,
    output_period=output_period,
    specific_outputs=specific_outputs,
    time_step_manager=tsmger,
    newton=newton,
    iteration_callbacks=[check_physics],
)

simulation.postprocess(time_unit="day")
