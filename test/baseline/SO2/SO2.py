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
from ComPASS.utils.grid import on_xmin, on_xmax
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.properties.enthalpies import build_pure_phase_enthalpy

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

ComPASS.set_output_directory_and_logfile(__file__)

simulation = ComPASS.load_eos("linear_water")
simulation.set_gravity(0)
rhocp = rho_reservoir * cp_reservoir
simulation.set_rock_volumetric_heat_capacity(rhocp)
simulation.set_molar_enthalpy_functions(
    build_pure_phase_enthalpy(volumetric_heat_capacity=rhocp),
)

nb_steps = int(L / ds) + 1
shape = [
    1,
] * 3
shape[axis] = nb_steps
extent = [
    L / ds,
] * 3
extent[axis] = L

grid = ComPASS.Grid(
    shape=shape,
    extent=extent,
)
on_the_left = on_xmin(grid)
on_the_right = on_xmax(grid)


def select_dirichlet_nodes():
    vertices = simulation.global_vertices()
    return on_the_left(vertices) | on_the_right(vertices)


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_thermal_conductivity=K_reservoir,
    cell_permeability=k_reservoir,
    cell_porosity=omega_reservoir,
)


# homogeneous reservoir initial state
X0 = simulation.build_state(p=p_reservoir, T=Tright)
simulation.all_states().set(X0)

# Dirichlet boundary conditions
dirichlet_nodes = simulation.dirichlet_node_states()
dirichlet_nodes.set(X0)
dirichlet_nodes.T[on_the_left(simulation.vertices())] = Tleft

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solver is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

output_period = final_time / nb_outputs
simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    output_period=output_period,
    time_step_manager=TimeStepManager(final_time / (10 * nb_steps), output_period),
)

print("\nWARNING: No check against anaytical solution was performed\n")

# if necessary simulation results can be directly postprocessed here
# from ComPASS.postprocess import postprocess
# postprocess(simulation.runtime.output_directory, convert_temperature=True)
