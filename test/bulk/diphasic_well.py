#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.grid import on_zmax, on_vertical_boundaries
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.units import *


# fmt: off
# A vertical well in the middle of a regular grid
# with dh * dh squared cells, H total thickness and, nv vertical layers
L = 1000               # half-width of the reservoir
dh = 200               # horizontal cell size
H = 200                # thickness of the reservoir / height of the well
nv = 5                 # number of vertical layers
rw = 0.1               # well radius (m)
ptop = 4 * MPa
Teps = 1               # temperature variation temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 9.81
omega = 0.15           # reservoir porosity
kh = 5e-14             # reservoir horizontal permeability in m^2
kvh = 0.001            # ratio of vertical permeability to horizontal permeability - low value to prevent natural convection
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour
wid = 0                # well id - could be any number
# fmt: on

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)

nh = 2 * int(L / dh)
assert nh > 0
grid = ComPASS.Grid(shape=(nh, nh, nv), extent=(2 * L, 2 * L, H), origin=(-L, -L, 0),)


def make_producer():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = wid
    well.operate_on_flowrate = Qm, 1.0 * bar
    well.produce()
    return [well]


def permeability():
    k = kh * np.eye(3, dtype="d")
    k[2, 2] *= kvh
    return k


simulation.init(
    mesh=grid,
    wells=make_producer,
    cell_porosity=omega,
    cell_permeability=permeability,
    cell_thermal_conductivity=K,
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
rho = simulation.liquid_molar_density(ptop, Ttop)
ztop = grid.origin[2] + grid.extent[2]


def set_states(states, z):
    states.set(X0)
    states.p[:] = ptop + rho * gravity * (ztop - z)


set_states(simulation.node_states(), simulation.vertices()[:, 2])
set_states(simulation.cell_states(), simulation.cell_centers()[:, 2])
set_states(simulation.fracture_states(), simulation.compute_fracture_centers()[:, 2])
simulation.reset_dirichlet_nodes(on_zmax(grid))
assert np.all(simulation.dirichlet_node_states().p == simulation.node_states().p)

simulation.close_well(wid)

tsmger = TimeStepManager(
    initial_timestep=1 * year, increase_factor=2.0, decrease_factor=0.2,
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
simulation.reset_dirichlet_nodes(on_vertical_boundaries(grid))

simulation.open_well(wid)
simulation.set_well_property(wid, imposed_flowrate=Qm)

simulation.standard_loop(
    initial_time=0, initial_timestep=hour, output_period=day, final_time=30 * day,
)

assert simulation.all_states().S[:, 0].max() > 0.1
