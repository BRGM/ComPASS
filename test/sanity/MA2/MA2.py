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
from ComPASS.timeloops import standard_loop

ComPASS.load_eos("liquid_water")

p_reservoir = 0  # no flow - linear water eos
Tleft, Tright = 33.0, 5.0  # K or deg C no matter
K_reservoir = 3  # W/m/K
k_reservoir = 1.0  # not relevant - dummy value
omega_reservoir = 0.5  # not relevant - conduction with rock properties

L = 1000.0  # domain length (m)
ds = 10.0  # step (m)
axis = 0  # axis {0: Ox, 1: Oy, 2: Oz}

#%% Simulation -----------------------------------------------------------------

ComPASS.set_gravity(0)

shape = [
    1,
] * 3
shape[axis] = int(L / ds) + 1
extent = [
    L / ds,
] * 3
extent[axis] = L
grid = ComPASS.Grid(
    shape=shape,
    extent=extent,
)


def select_dirichlet_nodes():
    x = ComPASS.global_vertices()[:, 0]
    return on_the_left(x) | on_the_right(x)


def set_boundary_conditions():
    def set_states(states, x):
        left = on_the_left(x)
        states.p[left] = pleft
        states.T[left] = Tleft
        right = on_the_right(x)
        states.p[right] = pright
        states.T[right] = Tright
        both = left | right
        states.context[both] = 1
        states.S[both] = 1
        if onecomp:
            states.C[both] = 1.0
        else:
            states.C[left] = (0, 1)
            states.C[right] = (1, 0)

    set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:, 0])


def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pright  # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.T[:] = Tright
        states.S[:] = 1
        if onecomp:
            states.C[:] = 1.0
        else:
            states.C[:] = (1, 0)

    set_states(ComPASS.node_states(), ComPASS.vertices()[:, 0])
    set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:, 0])


# %%% Simulation %%%

ComPASS.set_output_directory_and_logfile(__file__)

if onecomp:
    ComPASS.load_eos("liquid_water")
else:
    ComPASS.load_eos("water_with_tracer")

ComPASS.init(
    grid=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

set_initial_values()
set_boundary_conditions()

import ComPASS.utils.mpl_backends

# from ComPASS.utils.units import *
# from ComPASS.timeloops import standard_loop, TimeStepManager

from MA2_analytical import MA2_analytical as MA2


plt = ComPASS.utils.mpl_backends.import_pyplot(False)

analytical = MA2()

if plt:
    plt.clf()
    x = np.linspace(0, 4, 20)
    # t = (5e5, 1e6, 2e6)
    # t = (5e5, 1e6)
    t = (5e5, 1e6)
    for T in analytical(x, t):
        print(T)
        plt.plot(x, T)
    plt.savefig("MA2.png")
