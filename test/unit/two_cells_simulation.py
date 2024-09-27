import numpy as np

import ComPASS
from ComPASS.utils.units import *


def make(**init_args):
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

    simulation = ComPASS.load_physics("water2ph")

    grid = ComPASS.Grid(
        shape=(2, 1, 1),
        extent=(2, 1, 1),
        origin=(0, 0, 0),
    )

    def dirichlet_nodes():
        x = simulation.global_vertices()[:, 0]
        return (x == x.min()) & (x == x.max())

    args = {
        "mesh": grid,
        "cell_permeability": k_matrix,
        "cell_porosity": phi_matrix,
        "cell_thermal_conductivity": K_matrix,
        "set_dirichlet_nodes": dirichlet_nodes,
    }
    args.update(init_args)

    simulation.init(**args)

    initial_state = simulation.build_state(simulation.Context.liquid, p=pR, T=TR)
    simulation.all_states().set(initial_state)
    dirichlet = simulation.dirichlet_node_states()
    dirichlet.set(initial_state)
    x = simulation.vertices()[:, 0]
    left = x == x.min()
    dirichlet.p[left] = pL
    dirichlet.T[left] = TL

    return simulation
