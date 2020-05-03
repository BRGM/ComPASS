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


pleft, pright = 30 * MPa, 10 * MPa
Tleft, Tright = degC2K(60), degC2K(100)
omega_reservoir = 0.2  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K

Lx = 1000.0
nx = 100

onecomp = True

mu = 3e-4  # dynamic viscosity of pur water around 100°C (will change with temperature)
U = (k_reservoir / mu) * (pleft - pright) / Lx
print("Average Darcy velocity:", U * year, "m/year")
print(
    "                  i.e.: %.2f%%" % (100 * U * year / Lx),
    "of the simulation domain in one year.",
)
final_time = (Lx / 3) / (U / omega_reservoir)
print("Final time is set to: %.2f years" % (final_time / year))

grid = ComPASS.Grid(shape=(nx, 1, 1), extent=(Lx, Lx / nx, Lx / nx),)

on_the_left = lambda x: x <= grid.origin[0]
on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]


def select_dirichlet_nodes():
    x = simulation.global_vertices()[:, 0]
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

    set_states(simulation.dirichlet_node_states(), simulation.vertices()[:, 0])


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

    set_states(simulation.node_states(), simulation.vertices()[:, 0])
    set_states(simulation.cell_states(), simulation.compute_cell_centers()[:, 0])


# %%% Simulation %%%

ComPASS.set_output_directory_and_logfile(__file__)

if onecomp:
    simulation = ComPASS.load_eos("water2ph")
    simulation.lock_context(2)
else:
    assert False, "configuration temporarily desactivated"
    simulation = ComPASS.load_eos("water_with_tracer")

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

set_initial_values()
set_boundary_conditions()

cell_temperatures = []


def collect_node_temperature(iteration, t):
    if ComPASS.mpi.communicator().size > 1:
        if ComPASS.mpi.is_on_master_proc:
            print("WARNING - No output in parallel mode")
        return
    print("Collecting temperature at iteration", iteration)
    print("                           and time", t / year, "years")
    states = simulation.cell_states()
    cell_temperatures.append((t, np.copy(states.T)))


simulation.standard_loop(
    final_time=final_time,
    output_period=final_time / 50,
    initial_timestep=final_time / (10 * nx),
    output_callbacks=(collect_node_temperature,),
)

if ComPASS.mpi.communicator().size == 1:
    assert ComPASS.mpi.is_on_master_proc
    x = simulation.compute_cell_centers()[:, 0]
    with open(ComPASS.to_output_directory("cell_temperatures.csv"), "w") as f:
        s = ";".join(["%f" % (xi) for xi in x])
        print('"time (years)\\x";' + s, file=f)
        for tT in cell_temperatures:
            t, T = tT
            T = K2degC(T)
            print("%f;" % (t / year) + ";".join(["%f" % (Ti) for Ti in T]), file=f)
    import ComPASS.utils.mpl_backends as mpl_backends

    plt = mpl_backends.import_pyplot(False)
    if plt:
        plt.clf()
        for tT in cell_temperatures:
            t, T = tT
            T = K2degC(T)
            plt.plot(x, T)
        plt.xlabel("x in meters")
        plt.ylabel("temperature in Celsius degrees")
        plt.savefig(ComPASS.to_output_directory("cell_temperatures"))
