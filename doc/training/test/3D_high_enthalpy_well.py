# Documentation : https://charms.gitlabpages.inria.fr/ComPASS/

import itertools
import numpy as np

import ComPASS
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.units import *
from MeshTools import HexMesh
from MeshTools.utils import axis_extrusion

# fmt: off
# A vertical well in the middle of a regular square grid with x,y coordinates (0,0)
# 2nc cells along one side, 2nv vertical layers, H total thickness
# and a fracture in the middle
L = 1000               # side length of the reservoir
Lf = 600               # side length of the square fracture
nc = 20                # half the number of cells along reservoir side
H = 100                # thickness of the reservoir / heigth of the well
nv = 1                 # half the number of vertical layers
rw = 0.1               # well radius (m)
ptop = 4 * MPa         # pressure at the top of the reservoir
Ttop = degC2K(20)      #  temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 0
omega = 0.15           # reservoir porosity
kh = 5e-15             # reservoir horizontal permeability in m^2
kvh = 0.01             # ratio of vertical permeability to horizontal permeability - low value to prevent natural convection
kf = 1e-12             # fracture permeability
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour  # production flow rate
wid = 0                # well id - could be any number
epsilon = 1e-3
# fmt: on

# -- Mesh generation --------------------------------------
mesh = ComPASS.Grid(
    shape=(2 * nc, 2 * nc, 2 * nv),
    extent=(L, L, H),
    origin=(-L / 2, -L / 2, -H / 2),
)


# -- Setup simulation --------------------------------------

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)


simulation.init(
    mesh=mesh,
    cell_porosity=omega,
    cell_permeability=np.diag((kh, kh, kvh * kh)),
    cell_thermal_conductivity=K,
)

X0 = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(X0)
x, y, z = simulation.vertices().T
simulation.reset_dirichlet_nodes(
    (np.abs(x) > L / 2 - epsilon) | (np.abs(y) > L / 2 - epsilon)
)

outputs = np.logspace(np.log10(60), np.log10(5 * day))
outputs = [float(x) for x in outputs]  # bug workaround
simulation.standard_loop(
    initial_time=0,
    initial_timestep=30,
    specific_outputs=outputs,
    final_time=5 * day,
)

simulation.postprocess()
