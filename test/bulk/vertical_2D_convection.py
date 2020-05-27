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
from ComPASS.newton import Newton
from ComPASS.legacy_petsc import LegacyLinearSolver

simulation = ComPASS.load_eos("water2ph")
final_time = 20 * year
# maximum output period for small time step 1.5 year
output_period = 1.5 * year

Ttop = degC2K(25.0)
Tbot = degC2K(280.0)
H = 3 * km
Hcaprock = 0.1 * H
ptop = 1 * bar
pbot = ptop + 9.81 * 900.0 * H
kcap = 1e-18  # permeability caprok
kres = 1e-12  # permeability reservoir
phires = 0.15  # column porosity
Kres = 2  # bulk cell thermal conductivity W/m/K

nz = 101
grid = ComPASS.Grid(shape=(nz, 1, nz), extent=(H, H / nz, H), origin=(0.0, 0.0, -H),)


def dirichlet_boundaries():
    z = simulation.global_vertices()[:, 2]
    on_top = z == grid.origin[2] + grid.extent[2]
    on_bottom = z == grid.origin[2]
    return on_top | on_bottom


def cell_permeability():
    cell_centers = simulation.compute_global_cell_centers()
    zc = cell_centers[:, 2]
    nbcells = cell_centers.shape[0]
    # tensor array
    identity = np.eye(3)
    permeability = np.empty((nbcells, 3, 3), dtype=np.double)
    permeability[:] = kcap * identity
    reservoir = (zc < -Hcaprock) & (zc > -H + Hcaprock)
    permeability[reservoir] = kres * identity
    return permeability


ComPASS.set_output_directory_and_logfile(__file__)

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=dirichlet_boundaries,
    cell_permeability=cell_permeability,
    cell_porosity=phires,
    cell_thermal_conductivity=Kres,
)

X0 = simulation.build_state(simulation.Context.liquid, p=pbot, T=Tbot)
simulation.all_states().set(X0)


def apply_linear_gradients(states, z):
    states.p[:] = ptop - ((pbot - ptop) / H) * z
    states.T[:] = Ttop - ((Tbot - Ttop) / H) * z


apply_linear_gradients(simulation.all_states(), simulation.all_positions()[:, 2])
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(X0)
apply_linear_gradients(dirichlet, simulation.vertices()[:, 2])

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = LegacyLinearSolver(activate_direct_solver=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

simulation.standard_loop(
    initial_timestep=1 * year,
    final_time=final_time,
    newton=newton,
    output_period=output_period,
)
