import sys
import pickle
import yaml

import numpy as np

import ComPASS
from ComPASS.utils.units import *
import ComPASS.mpi as mpi
import MeshTools as MT

R = 1000  # radius of the angular sector
theta = np.pi / 12  # angle of the angular sector (in radians)
rw = 0.1  # well radius
H = 100  # reservoir thickness
omega = 0.2  # reservoir porosity
k = 1e-14  # reservoir permeability in m^2
K = 0  # bulk thermal conductivity in W/m/K
qw = 14.0  # total production flow rate (for the whole reservoir i.e. theta = 2pi)
pres = 30 * bar  # reservoir pressure
Sgres = 1.0 - 0.65  # reservoir saturation
wid = 0  # well id

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)
simulation.set_rock_volumetric_heat_capacity(2.65e6)  # SI units J/m^3/°C

vertices = np.array(
    [
        [0, 0, 0],
        [R, 0, 0],
        [R * np.cos(theta), R * np.sin(theta), 0],
        [0, 0, H],
        [R, 0, H],
        [R * np.cos(theta), R * np.sin(theta), H],
    ]
)
cellnodes = np.array([[0, 1, 2, 3, 4, 5],])
mesh = MT.WedgeMesh.make(vertices, cellnodes)

epsilon = 1e-4 * R  # tolerance value to select nodes (boundary conditions...)
angular_qw = (theta / (2 * np.pi)) * qw


def make_well():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = wid
    well.operate_on_flowrate = angular_qw, 1 * bar
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
X0 = simulation.build_state(simulation.Context.diphasic, p=pres, Sg=Sgres)
simulation.all_states().set(X0)
# Set boundary conditions
vertices = simulation.vertices()
x, y = vertices[:, 0], vertices[:, 1]
simulation.reset_dirichlet_nodes((x ** 2 + y ** 2) > R ** 2 - epsilon)

assert mpi.is_on_master_proc  # test is not to be run in parallel (one cell !)


def check(tick):
    wellhead = simulation.get_wellhead(wid)  # wellhead could be None in parallel
    assert np.allclose(angular_qw, wellhead.molar_flowrate)
    assert tick.iteration < 55  # no more than 55 iterations
    pass


simulation.standard_loop(
    initial_timestep=1e-4 * day,
    final_time=day,
    output_period=hour,
    iteration_callbacks=[check,],
)

# if necessary simulation results can be directly postprocessed here
simulation.postprocess()
