#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

import numpy as np

ComPASS.load_eos("linear_water")
assert not ComPASS.has_energy_transfer_enabled()

rhof = 1e3  # specific mass in kg/m^3
cpf = 4200  # specific heat in J/kg/K
rhofcpf = rhof * cpf  # volumetric heat capacity
muf = 1e-3  # viscosity Pa.s
# get_fluid_properties does not exist anymore, use set_property_functions
fluid_properties = ComPASS.get_fluid_properties()
fluid_properties.specific_mass = rhof
fluid_properties.volumetric_heat_capacity = rhofcpf
fluid_properties.dynamic_viscosity = muf

pleft, pright = 30 * MPa, 10 * MPa
Tleft, Tright = degC2K(60), degC2K(100)
omega_reservoir = 0.2  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K

Lx = 1000.0
nx = 100

mu = 3e-4  # dynamic viscosity of pur water around 100°C (will change with temperature)
U = (k_reservoir / mu) * (pleft - pright) / Lx
print("Average Darcy velocity:", U * year, "m/year")
print(
    "                  i.e.: %.2f%%" % (100 * U * year / Lx),
    "of the simulation domain in one year.",
)
final_time = (Lx / 3) / (U / omega_reservoir)
print("Final time is set to: %.2f years" % (final_time / year))

grid = ComPASS.Grid(
    shape=(nx, 1, 1),
    extent=(Lx, Lx / nx, Lx / nx),
)

on_the_left = lambda x: x <= grid.origin[0]
on_the_right = lambda x: x >= grid.origin[0] + grid.extent[0]


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
        states.C[both] = 1.0

    set_states(ComPASS.dirichlet_node_states(), ComPASS.vertices()[:, 0])


def set_initial_values():
    def set_states(states, x):
        states.context[:] = 1
        states.p[:] = pright  # pleft + (pright - pleft) * (x - grid.origin[0]) / Lx
        states.T[:] = Tright
        states.S[:] = 1
        states.C[:] = 1.0

    set_states(ComPASS.node_states(), ComPASS.vertices()[:, 0])
    set_states(ComPASS.cell_states(), ComPASS.compute_cell_centers()[:, 0])


# %%% Simulation %%%

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid=grid,
    set_dirichlet_nodes=select_dirichlet_nodes,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

set_initial_values()
set_boundary_conditions()

standard_loop(
    final_time=final_time,
    output_period=final_time / 50,
    initial_timestep=final_time / (10 * nx),
)
