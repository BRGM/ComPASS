import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.utils.angular_sector import extruded_sector

R = 1000  # radius of the angular sector
theta = np.pi / 3  # angle of the angular sector (in radians)
rw = 0.1  # well radius
ds = 1, 100  # target edge size around the origin and at the boundary
H = 100  # reservoir thickness
nb_layers = 1  # number of horizontal layers
omega = 0.2  # reservoir porosity
k = 1e-12  # reservoir permeability in m^2
K = 20  # bulk thermal conductivity in W/m/K
rock_heat_capacity = 1e3  # J/kg/K
rock_density = 2500  # kg/m3
qw = 10.0  # total injection flow rate kg/s
pres = 50 * bar  # reservoir pressure
Tres = degC2K(170)
Tinj = degC2K(160)
gravity = 0

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(rock_heat_capacity * rock_density)

mesh = extruded_sector(R, theta, ds, np.tile(H / nb_layers, nb_layers))

epsilon = 0.01 * ds[0]  # tolerance value to select nodes (boundary conditions...)


def make_well():
    well = simulation.create_vertical_well((0, 0), rw)
    well.operate_on_flowrate = -(theta / (2 * np.pi)) * qw, np.inf
    well.inject(Tinj)
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
simulation.reset_dirichlet_nodes((x**2 + y**2) > R**2 - epsilon)

# final_time=1e9
timestep = 1.67e7
simulation.standard_loop(
    fixed_timestep=timestep,
    final_time=60 * timestep,
    output_period=timestep,
)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess(time_unit="second")
