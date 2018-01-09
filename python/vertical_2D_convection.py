#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.load_eos('water2ph')

Ttop = degC2K( 25. )
Tbot = degC2K( 350. )
Hcaprock = 0.3 * km
H = 3 * km
ptop = 1 * bar
pbot = ptop + 9.81 * 1000. * H
kcap = 1E-17 # permeability caprok
kres = 1E-12 # permeability reservoir

ncells = 31
grid = ComPASS.Grid(
    shape = (ncells, 1, ncells),
    extent = (H, H/ncells, H),
    origin = (0., 0., -H),
)

def dirichlet_boundaries():
    vertices = np.rec.array(ComPASS.global_vertices())
    on_top = (vertices.z == grid.origin[2] + grid.extent[2])
    on_bottom = (vertices.z == grid.origin[2])
    return on_top | on_bottom

def cell_permeability():
    cell_centers = ComPASS.compute_global_cell_centers()
    zc = cell_centers[:, 2]
    nbcells = cell_centers.shape[0]
    # tensor array
    permeability = np.empty((nbcells, 3, 3), dtype=np.double)
    permeability[:] = np.eye(3)
    permeability[zc >= -Hcaprock]*= kcap
    permeability[zc <  -Hcaprock]*= kres
    return permeability

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    set_dirichlet_nodes = dirichlet_boundaries,
    cells_permeability = cell_permeability
)

def set_states(states, z):
    states.context[:] = 2
    states.p[:] = ptop - ((pbot - ptop) / H) * z
    states.T[:] = Ttop - ((Tbot - Ttop) / H) * z
    states.S[:] = [0, 1]
    states.C[:] = 1.
set_states(ComPASS.dirichlet_node_states(), np.rec.array(ComPASS.vertices()).z)
set_states(ComPASS.node_states(), np.rec.array(ComPASS.vertices()).z)
set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,2])

standard_loop(initial_timestep= 30 * day, final_time = 1 * ky, output_period = 10 * year)
