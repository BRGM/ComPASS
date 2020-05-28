#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# This comes from:
# @techreport{Thiery2012,
# author = {Thi{\'{e}}ry, Dominique},
# file = {:D$\backslash$:/work/data/biblio/Thi{\'{e}}ry - 2012 - Validation des calculs de transport d'{\'{e}}nergie dans le logiciel MARTHE.pdf:pdf},
# institution = {BRGM},
# keywords = {charms cas test},
# mendeley-tags = {charms cas test},
# title = {{Validation des calculs de transport d'{\'{e}}nergie dans le logiciel MARTHE}},
# url = {http://www.brgm.fr/sites/default/brgm/logiciels/marthe/nt{\_}eau{\_}2011{\_}02{\_}marthe{\_}valid{\_}thermique.pdf},
# year = {2012}
# }

import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.newton import Newton
from ComPASS.legacy_petsc import LegacyLinearSolver

ComPASS.load_eos("linear_water")

p_reservoir = 0  # no flow - linear water eos
Tleft, Tright = 33.0, 5.0  # K or deg C no matter
K_reservoir = 3  # bulk thermal conductivity W/m/K
rho_reservoir = 2600.0  # rock density kg/m3
cp_reservoir = 900  # rock specific heat J/K/kg
k_reservoir = 1.0  # not relevant - dummy value
omega_reservoir = 0.5  # not relevant - conduction with rock properties

axis = 0  # axis {0: Ox, 1: Oy, 2: Oz}
L = 1000.0  # domain length (m)
ds = 10.0  # step (m)
final_time = 1e6  # seconds
nb_outputs = 20

#%% Simulation -----------------------------------------------------------------

ComPASS.set_gravity(0)
rhocp = rho_reservoir * cp_reservoir
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.volumetric_heat_capacity = rhocp
ComPASS.set_rock_volumetric_heat_capacity(rhocp)

nb_steps = int(L / ds) + 1
shape = [1,] * 3
shape[axis] = nb_steps
extent = [L / ds,] * 3
extent[axis] = L
grid = ComPASS.Grid(shape=shape, extent=extent,)


on_the_left = lambda xyz: xyz[:, axis] <= grid.origin[axis]
on_the_right = lambda xyz: xyz[:, axis] >= (grid.origin[axis] + grid.extent[axis])


def select_dirichlet_nodes():
    coords = ComPASS.global_vertices()
    return on_the_left(coords) | on_the_right(coords)


def set_boundary_conditions():
    def set_states(states, xyz):
        left = on_the_left(xyz)
        states.p[left] = p_reservoir
        states.T[left] = Tleft
        right = on_the_right(xyz)
        states.p[right] = p_reservoir
        states.T[right] = Tright
        both = left | right
        states.context[both] = 1
        states.S[both] = 1
        states.C[both] = 1.0

    set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices())


def set_initial_values():
    def set_states(states):
        states.context[:] = 1
        states.p[:] = p_reservoir
        states.T[:] = Tright
        states.S[:] = 1
        states.C[:] = 1.0

    set_states(ComPASS.node_states())
    set_states(ComPASS.cell_states())


ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_thermal_conductivity=K_reservoir,
    cell_permeability=k_reservoir,
    cell_porosity=omega_reservoir,
)

set_initial_values()
set_boundary_conditions()

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = LegacyLinearSolver(simulation, activate_direct_solver=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

output_period = final_time / nb_outputs
standard_loop(
    final_time=final_time,
    newton=newton,
    output_period=output_period,
    time_step_manager=TimeStepManager(final_time / (10 * nb_steps), output_period),
)
