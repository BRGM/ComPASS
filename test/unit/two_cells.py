# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop


# initial reservoir pressure
pL, pR = 2.0 * bar, 1.0 * bar
# initial reservoir temperature - convert Celsius degrees to Kelvin degrees
TL, TR = (
    degC2K(30),
    degC2K(20),
)
# column permeability in m^2 (low permeability -> bigger time steps)
k_matrix = 1e-12
# column porosity
phi_matrix = 0.15
# bulk thermal conductivity in W/m/K
K_matrix = 2.0

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(2, 1, 1), extent=(2, 1, 1), origin=(0, 0, 0),)


def dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
    return (x == x.min()) & (x == x.max())


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_dirichlet_nodes=dirichlet_nodes,
)

initial_state = simulation.build_state(simulation.Context.liquid, p=pR, T=TR)
simulation.all_states().set(initial_state)
dirichlet = simulation.dirichlet_node_states()
dirichlet.set(initial_state)
x = simulation.vertices()[:, 0]
left = x == x.min()
dirichlet.p[left] = pL
dirichlet.T[left] = TL

current_time = standard_loop(
    simulation, initial_timestep=1, nitermax=2, output_period=1,
)

# simulation results can be directly postprocessed here
# from ComPASS.postprocess import postprocess
# postprocess(simulation.runtime.output_directory)
