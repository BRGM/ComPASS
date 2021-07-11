import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton

import hexmesh

from my_kr import kr_functions

R = 1000  # radius of the angular sector
theta = np.pi / 180  # angle of the angular sector (in radians) must be small
rw = 0.1  # well radius
nr = 100  # number of radius discretization points
Htot = 100  # reservoir thickness
H = rw  # thickness used for simulation
omega = 0.15  # reservoir porosity
k = 24e-14  # reservoir permeability in m^2
K = 0  # bulk thermal conductivity in W/m/K
qw = (
    H / Htot
) * 16.7  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 30 * bar  # reservoir pressure
Sgres = 1.0 - 0.65  # reservoir saturation

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)
simulation.set_rock_volumetric_heat_capacity(2.0e6)  # SI units J/m^3/°C

a = hexmesh.geometric_progression(nr - 1, rw, R - rw)
r = [0, rw]
for _ in range(nr - 1):
    r.append(r[-1] * a)
mesh = hexmesh.make_mesh(rw + np.cumsum(r), theta, H)

epsilon = 0.1 * rw  # tolerance value to select nodes (boundary conditions...)

simulation.init(
    mesh=mesh, cell_porosity=omega, cell_permeability=k, cell_thermal_conductivity=K,
)

# Set the kr functions after initialization
simulation.set_kr_functions(kr_functions)

# Set initial values
X0 = simulation.build_state(simulation.Context.diphasic, p=pres, Sg=Sgres)
simulation.all_states().set(X0)
# Set boundary conditions
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
simulation.reset_dirichlet_nodes((x ** 2 + y ** 2) > R ** 2 - epsilon)

# Neumann at the well face
mass_flux = qw * (theta / (2 * np.pi))
Neumann = ComPASS.NeumannBC(-mass_flux, compute_heat_flux=True)
face_centers = simulation.face_centers()
x, y, z = [face_centers[:, j] for j in range(3)]
well_face = (x ** 2 + y ** 2) <= rw ** 2 + epsilon
well_face_nodes = simulation.facenodes(well_face)
simulation.set_Neumann_faces(well_face, Neumann)

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 20, lsolver)

simulation.standard_loop(
    initial_timestep=1e-5 * day, final_time=day, output_period=minute, newton=newton
)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess(time_unit="day")
