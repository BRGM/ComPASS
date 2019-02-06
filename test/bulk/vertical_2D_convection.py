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
# KSP failures with iterative solver
ComPASS.activate_direct_solver = True
final_time = 20 * year
# maximum output period for small time step 1.5 year
output_period = 1.5 * year

Ttop = degC2K( 25. )
Tbot = degC2K( 280. )
H = 3 * km
Hcaprock = 0.1 * H
ptop = 1 * bar
pbot = ptop + 9.81 * 900. * H
kcap = 1E-18         # permeability caprok
kres = 1E-12         # permeability reservoir
phires = 0.15        # column porosity
Kres = 2             # bulk cell thermal conductivity W/m/K

nz = 101
grid = ComPASS.Grid(
    shape = (nz, 1, nz),
    extent = (H, H/nz, H),
    origin = (0., 0., -H),
)

def dirichlet_boundaries():
    z = ComPASS.global_vertices()[:, 2]
    on_top = (z == grid.origin[2] + grid.extent[2])
    on_bottom = (z == grid.origin[2])
    return on_top | on_bottom

def cell_permeability():
    cell_centers = ComPASS.compute_global_cell_centers()
    zc = cell_centers[:, 2]
    nbcells = cell_centers.shape[0]
    # tensor array
    identity = np.eye(3)
    permeability = np.empty((nbcells, 3, 3), dtype=np.double)
    permeability[:] = kcap * identity
    reservoir = ( zc < -Hcaprock ) & ( zc > -H + Hcaprock )
    permeability[reservoir] = kres * identity
    return permeability

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = grid,
    set_dirichlet_nodes = dirichlet_boundaries,
    cell_permeability = cell_permeability,
    cell_porosity = phires,
    cell_thermal_conductivity = Kres,
)

def set_states(states, z):
    states.context[:] = 2
    states.p[:] = ptop - ((pbot - ptop) / H) * z
    states.T[:] = Ttop - ((Tbot - Ttop) / H) * z
    states.S[:] = [0, 1]
    states.C[:] = 1.
set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:, 2])
set_states(ComPASS.node_states(), ComPASS.vertices()[:, 2])
set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:,2])

standard_loop(
    initial_timestep = 1 * year,
    final_time = final_time,
    output_period = output_period,
)
