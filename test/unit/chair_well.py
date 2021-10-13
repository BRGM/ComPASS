#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.utils.grid import on_vertical_boundaries

Lx, Ly, Lz = 600.0, 600.0, 1000.0
Ox, Oy, Oz = 0, 0, 0
nx, ny, nz = 12, 12, 6

pres = 1e7
Tres = degC2K(60)
rw = 0.1
Qw = 1.0
dx = Lx / nx
lhb = 100.0  # length of the horizontal branch
assert abs(lhb % dx) < 1e-5
nhb = int(lhb / dx)  # length of the horizontal branch in cells

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)


def create_well():
    vertices = simulation.global_vertices()
    v = vertices.reshape((nx + 1, ny + 1, nz + 1, 3), order="F")
    indices = np.reshape(
        np.arange(vertices.shape[0]), (nx + 1, ny + 1, nz + 1), order="F"
    )
    ic, jc, kc = (nx + 1) // 2, (ny + 1) // 2, (nz + 1) // 2  # center indices
    vertical_branch = indices[ic, jc, :]
    horizontal_branch = indices[ic : ic + nhb + 1, jc, kc]
    small_vertical_branch = indices[ic + nhb, jc, : kc + 1]
    for branch in [vertical_branch, horizontal_branch, small_vertical_branch]:
        branch.shape = -1, 1
    # well segments must be oriented from wellhead downwards
    segments = np.vstack(
        [
            # nodes are from bottom to top
            np.hstack([vertical_branch[1:], vertical_branch[:-1]]),
            # order of nodes is ok
            np.hstack([horizontal_branch[:-1], horizontal_branch[1:]]),
            # nodes are from bottom to top
            np.hstack([small_vertical_branch[1:], small_vertical_branch[:-1]]),
        ]
    )
    return simulation.create_well_from_segments(segments, well_radius=rw)


def make_producer():
    wid = 0  # well id - could be any number
    mswell = create_well()
    mswell.operate_on_flowrate = Qw, bar
    mswell.produce()
    mswell.id = wid
    return [mswell]


simulation.init(
    mesh=grid,
    wells=make_producer,
    cell_porosity=0.1,
    cell_permeability=1e-12,
    cell_thermal_conductivity=2,
)

X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
simulation.reset_dirichlet_nodes(on_vertical_boundaries(grid))

simulation.standard_loop(
    initial_time=0, initial_timestep=hour, output_period=day, final_time=day,
)

simulation.postprocess(time_unit="day", convert_temperature=True)
