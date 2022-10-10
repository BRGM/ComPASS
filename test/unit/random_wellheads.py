#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import copy
from collections import defaultdict

import numpy as np

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.units import *

# fmt: off
np.random.seed(12345)  # set the seed to always have the same well pattern
nb_random_wells = 10   # nb_random_wells producers + nb_random_wells injectors
pres = 15.0 * MPa      # initial reservoir pressure
Tres = degC2K(60)      # temperature of the reservoir
k_reservoir = 1e-12    # reservoir permeability in m^2
K_reservoir = 2        # bulk thermal conductivity in W/m/K
phi_reservoir = 0.15   # reservoir porosity
gravity = 0
rw = 0.1               # well radius
Qw = 100.0             # well imposed flowrate
# fmt: on

Lx, Ly, Lz = 1000.0, 1000.0, 100.0
nx, ny, nz = 40, 40, 1  # discretization

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(gravity)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
)


def make_producer(x, y, Qw):
    well = simulation.create_vertical_well((x, y), rw)
    well.operate_on_flowrate = Qw, -np.inf
    well.produce()
    return well


def make_injector(x, y, Qw):
    well = simulation.create_vertical_well((x, y), rw)
    well.operate_on_flowrate = -Qw, np.inf
    well.inject(Tres)
    return well


def set_wells(n):
    def make_wells(id_offset, factory):
        xy = np.random.random(2 * n)
        xy.shape = -1, 2
        # We do not want well to be on the grid boundaries
        # so that they do not interfere with dirichlet boundary conditions
        assert nx > 2
        dx = Lx / nx
        xy[:, 0] = xy[:, 0] * (Lx - 2 * dx) + dx
        assert ny > 2
        dy = Ly / ny
        xy[:, 1] = xy[:, 1] * (Ly - 2 * dy) + dy
        wells = []
        for wid, pos in enumerate(xy):
            wells.append(factory(pos[0], pos[1], Qw))
            wells[-1].id = id_offset + wid
        return wells

    return make_wells(0, make_producer) + make_wells(n, make_injector)


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    wells=lambda: set_wells(nb_random_wells),
    cell_porosity=phi_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
    display_well_ids=True,
)

for wid in range(2 * nb_random_wells):
    wh = simulation.get_wellhead(wid)
    assert wh is not None
    print(f"Well {wid} has pressure: {wh.pressure}")
