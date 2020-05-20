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
from ComPASS.utils.units import *


# fmt: off
# A vertical well in the middle of a regular grid
# with dh * dh squared cells, H total thickness and, nv vertical layers
L = 1000               # half-width of the reservoir
H = 200                # thickness of the reservoir / height of the well
dh = 25               # horizontal cell size
Rf = 600               # horizontal fracture radius, less than L - 0.5* dh so that the fracture does not touch the boundary
zf = 0.8 * H           # horizontal fracture elevation, must be between 0 and H, and such that zf % (H / nv) == 0
nv = 50                # number of vertical layers
rw = 0.1               # well radius (m)
ptop = 4 * MPa
Teps = 1               # temperature variation temperature at the top of the reservoir is set to Tsat(ptop) - Teps
gravity = 9.81
omega_matrix = 0.15    # reservoir porosity
omega_frac = 0.5       # fracture porosity
k_matrix = 1e-16       # reservoir permeability in m^2 - low value to prevent natural convection
k_frac = 1e-11         # fractured zone permeability in m^2
th_frac = 0.1          # fractured zone thickness
K = 2                  # bulk thermal conductivity in W/m/K
Qm = 200 * ton / hour
wid = 0                # well id - could be any number
# fmt: on

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)
simulation.set_fracture_thickness(th_frac)


nh = 2 * int(L / dh)
assert nh > 0
grid = ComPASS.Grid(shape=(nh, nh, nv), extent=(2 * L, 2 * L, H), origin=(-L, -L, 0),)


def make_producer():
    well = simulation.create_vertical_well((0, 0), rw)
    well.id = wid
    well.operate_on_flowrate = Qm, 1.0 * bar
    well.produce()
    return [well]


def fracture_faces():
    fc = simulation.compute_global_face_centers()
    x, y, z = [fc[:, j] for j in range(3)]
    dz = H / nv
    assert dz > 0
    assert abs(zf % dz) < 1e-5, "fracture must be located on cell boundary"
    return (np.fabs(z - zf) < (0.25 * dz)) & (np.sqrt(x ** 2 + y ** 2) < Rf)


simulation.init(
    mesh=grid,
    wells=make_producer,
    fracture_faces=fracture_faces,
    cell_porosity=omega_matrix,
    cell_permeability=k_matrix,
    fracture_porosity=omega_frac,
    fracture_permeability=k_frac,
    cell_thermal_conductivity=K,
    fracture_thermal_conductivity=K,
)

Ttop = simulation.Tsat(ptop) - Teps
print(f"temperature at reservoir top: {K2degC(Ttop):.1f} °C")
print(f"   delta pressure to boiling: {(simulation.Psat(Ttop) - ptop)/bar:.1f} bar")

## ------------------------------
## First step - equilibrium state
## ------------------------------

simulation.reset_dirichlet_nodes(on_zmax(grid))
X0 = simulation.build_state(simulation.Context.liquid, p=ptop, T=Ttop)
simulation.all_states().set(X0)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(X0)

simulation.close_well(wid)

simulation.standard_loop(
    initial_timestep=10 * day,
    final_time=10 * year,  # more than enough to reach pressure equilibrium
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
    initial_time=0, initial_timestep=hour, output_period=day, final_time=2 * year
)
