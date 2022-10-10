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
from ComPASS.utils.petrel import PetrelGrid

# fmt: off
pres = 20. * MPa            # initial reservoir pressure
Tres = degC2K( 70. )        # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 30. )  # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour      # production flowrate
omega_reservoir = 0.15      # reservoir porosity
k_reservoir = 1E-12         # reservoir permeability in m^2
K_reservoir = 2             # bulk thermal conductivity in W/m/K
g = 9.81
# fmt: on

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(g)


grid = PetrelGrid("sample.grdecl")


def select_dirichlet():
    xyz = simulation.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    return (x == x.min()) | (y == y.min()) | (x == x.max()) | (y == y.max())


def make_wells():
    xyz = simulation.global_vertices()
    x, y, z = [xyz[:, i] for i in range(3)]
    Lx, Ly, Lz = (x.max() - x.min(), y.max() - y.min(), z.max() - z.min())
    interwell_distance = 0.5 * Ly
    Cx, Cy = (x.min() + 0.5 * Lx, y.min() + 0.5 * Ly)
    producer = simulation.create_vertical_well((Cx, Cy - 0.5 * interwell_distance))
    producer.operate_on_flowrate = Qm, 1.0 * bar
    producer.produce()
    injector = simulation.create_vertical_well((Cx, Cy + 0.5 * interwell_distance))
    injector.operate_on_flowrate = Qm, pres + 100.0 * MPa
    injector.inject(Tinjection)
    return (producer, injector)


simulation.init(
    mesh=grid.mesh,
    set_dirichlet_nodes=select_dirichlet,
    wells=make_wells,
    cell_porosity=omega_reservoir,
    cell_permeability=grid.permeability,
    cell_thermal_conductivity=K_reservoir,
)

zref = 3000
pref = 1 * bar  # zref = 3000m
rhoref = 1000.0  # kg / m3

X0 = simulation.build_state(simulation.Context.liquid, p=rhoref * g * zref, T=Tres)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)


def set_pressure(states, xyz):
    z = xyz[:, 2]
    states.p[:] = pref + rhoref * g * (zref - z)


for states, xyz in [
    (simulation.dirichlet_node_states(), simulation.vertices()),
    (simulation.node_states(), simulation.vertices()),
    (simulation.cell_states(), simulation.compute_cell_centers()),
]:
    set_pressure(states, xyz)

simulation.standard_loop(
    initial_timestep=30 * day,
    final_time=30 * year,
    output_period=year,
)
