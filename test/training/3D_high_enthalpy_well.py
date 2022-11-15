#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import itertools
import numpy as np

import ComPASS
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.units import *
from MeshTools import HexMesh
from MeshTools.utils import axis_extrusion

import vtkwriters as vtkw

# fmt: off
# A vertical well in the middle of a regular square grid with x,y coordinates (0,0)
# 2nc cells along one side, 2nv vertical layers, H total thickness
# and a fracture in the middle
L = 1000               # side length of the reservoir
Lf = 600               # side length of the square fracture
nc = 10                # half the number of cells along reservoir side
H = 100                # thickness of the reservoir / heigth of the well
nv = 4                 # half the number of vertical layers
rw = 0.1               # well radius (m)
f_thickness = 0.1      # fracture thickness (m)
ptop = 4 * MPa         # oressure at the top of the reservoir
Teps = 1               #  temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 9.81
omega = 0.15           # reservoir porosity
kh = 5e-15             # reservoir horizontal permeability in m^2
kvh = 0.01             # ratio of vertical permeability to horizontal permeability - low value to prevent natural convection
kf = 1e-12             # fracture permeability
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour  # production flow rate
wid = 0                # well id - could be any number
epsilon = 1e-3         # tolerance for entity selection
# fmt: on

# -- Mesh generation --------------------------------------
assert nc > 0
assert nv > 0
x = (L / (2 * nc)) * np.arange(-nc, nc + 1)
vertices = np.array(list(itertools.product(x, x, (-H / 2,))), dtype="d")
cells = np.array([1, 0, 2 * nc + 1, 2 * nc + 2])  # one cell
cells = np.vstack([cells + k for k in range(2 * nc)])  # one column along Oy
cells = np.vstack([cells + k * (2 * nc + 1) for k in range(2 * nc)])  # all cells
# layers ticknesses
thicknesses = np.tile(H / (2 * nv), 2 * nv)
vertices, cells = axis_extrusion(vertices, cells, offsets=thicknesses)
mesh = HexMesh.make(vertices, cells)

# -- Setup simulation --------------------------------------

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)


def make_producer():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = wid
    well.operate_on_flowrate = Qm, 1.0 * bar
    well.produce()
    return [well]


simulation.set_fracture_thickness(f_thickness)


def select_fractures():
    face_centers = simulation.compute_global_face_centers()
    xc, yc, zc = face_centers.T
    return (
        (np.abs(xc) < Lf / 2 + epsilon)
        & (np.abs(yc) < Lf / 2 + epsilon)
        & (np.abs(zc) < epsilon)
    )


simulation.init(
    mesh=mesh,
    wells=make_producer,
    cell_porosity=omega,
    cell_permeability=np.diag((kh, kh, kvh * kh)),
    cell_thermal_conductivity=K,
    fracture_faces=select_fractures,
    fracture_porosity=omega,
    fracture_permeability=kf,
    fracture_thermal_conductivity=K,
    well_model="two_phases",
)

Ttop = simulation.Tsat(ptop) - Teps
print(f"temperature at reservoir top: {K2degC(Ttop):.1f} °C")
print(f"   delta pressure to boiling: {(simulation.Psat(Ttop) - ptop)/bar:.1f} bar")

## ------------------------------
## First step - equilibrium state
## ------------------------------

X0 = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
# approximate rho
rho = simulation.liquid_volumetric_mass_density(ptop, Ttop)
ztop = H / 2


def set_states(states, z):
    states.set(X0)
    states.p[:] = ptop + rho * gravity * (ztop - z)


set_states(simulation.all_states(), simulation.all_positions()[:, 2])
x, y, z = simulation.vertices().T
simulation.reset_dirichlet_nodes(np.abs(z - H / 2) < epsilon)
set_states(simulation.dirichlet_node_states(), z)
assert np.all(simulation.dirichlet_node_states().p == simulation.node_states().p)

simulation.close_well(wid)

tsmger = TimeStepManager(
    initial_timestep=1 * year,
    increase_factor=2.0,
    decrease_factor=0.2,
)

simulation.standard_loop(
    final_time=10 * year,  # more than enough to reach pressure equilibrium
    time_step_manager=tsmger,
    no_output=True,
)

## ------------------------
## Second step - production
## ------------------------

# Dirichlet nodes will be locked to their equilibrium values
simulation.reset_dirichlet_nodes(
    (np.abs(x) > L / 2 - epsilon) | (np.abs(y) > L / 2 - epsilon)
)

simulation.open_well(wid)
simulation.set_well_property(wid, imposed_flowrate=Qm)

simulation.standard_loop(
    initial_time=0,
    initial_timestep=hour,
    output_period=10 * day,
    specific_outputs=[k * day for k in range(1, 10)],
    final_time=100 * day,
)

simulation.postprocess()
