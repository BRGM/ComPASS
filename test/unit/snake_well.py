# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from itertools import cycle
import sys
import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop


nx = ny = nz = 4

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")
simulation.set_gravity(0)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(nx, ny, nz),
    origin=(0, 0, 0),
)

hpath = [
    (1, 1),
    (2, 1),
    (3, 1),
    (3, 2),
    (3, 3),
    (2, 3),
    (1, 3),
    (1, 2),
]


def make_well():
    nv = (nx + 1, ny + 1, nz + 1)
    vid = np.reshape(np.arange(np.prod(nv)), nv, order="F")
    wpath = []
    k = nz
    n = 0
    i0, j0 = hpath[-1]
    for i, j in cycle(hpath):
        if n == 0:
            k -= 1
            wpath.append((vid[i0, j0, k], vid[i0, j0, k + 1]))
            if k == 0:
                break
            n = 7
        n -= 1
        # order of vid is important here: always from bottom to well-head
        wpath.append((vid[i, j, k], vid[i0, j0, k]))
        i0, j0 = i, j
    vertices = simulation.global_vertices()
    for va, vb in wpath:
        print(vertices[va], "->", vertices[vb])
    well = simulation.create_well_from_segments(wpath)
    well.operate_on_flowrate = 1.0, 1.0
    well.produce()
    return [well]


simulation.init(
    mesh=grid,
    wells=make_well,
    cell_permeability=1,
    cell_porosity=0.5,
    cell_thermal_conductivity=1,
)


X0 = simulation.build_state(simulation.Context.liquid, p=MPa, T=300)
simulation.all_states().set(X0)


simulation.standard_loop(
    nitermax=0,
    fixed_timestep=1,
)

# simulation results can be directly postprocessed here
from ComPASS.postprocess import postprocess

postprocess(simulation.runtime.output_directory)
