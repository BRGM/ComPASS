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
from ComPASS.utils.grid import grid_center

# fmt: off
pres = 20. * MPa            # initial reservoir pressure
Tres = degC2K( 70. )        # initial reservoir temperature - convert Celsius to Kelvin degrees
Tinjection = degC2K( 30. )  # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour      # production flowrate
k_matrix = 1E-13            # matrix permeability in m^2
omega_matrix = 0.15         # matrix porosity
K_matrix = 2                # bulk thermal conductivity in W/m/K
k_fracture = 1E-12          # fracture permeability in m^2
omega_fracture = 0.5        # fracture porosity
K_fracture = 2              # bulk thermal conductivity in W/m/K
# fmt: on

Lx, Ly, Lz = 100.0, 100.0, 20.0
Ox, Oy, Oz = 0.0, -0.5 * Ly, -0.5 * Lz
nx, ny, nz = 11, 11, 4
dx = Lx / nx
dz = Lz / nz

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)


def select_fractures():
    face_centers = simulation.compute_global_face_centers()
    # select horizontal fault axis in the middle of the simulation domain i.e. z=0
    return np.abs(face_centers[:, 2]) < 0.25 * dz


def make_wells():
    interwell_distance = Lx / 3
    Cx, Cy, Cz = grid_center(grid)
    producer = simulation.create_vertical_well((Cx - 0.5 * interwell_distance, Cy))
    producer.operate_on_flowrate = Qm, 1.0 * bar
    producer.produce()
    injector = simulation.create_vertical_well((Cx + 0.5 * interwell_distance, Cy))
    injector.operate_on_flowrate = Qm, pres + 100.0 * MPa
    injector.inject(Tinjection)
    return (producer, injector)


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    wells=make_wells,
    fracture_faces=select_fractures,
    cell_porosity=omega_matrix,
    cell_permeability=k_matrix,
    cell_thermal_conductivity=K_matrix,
    fracture_porosity=omega_fracture,
    fracture_permeability=k_fracture,
    fracture_thermal_conductivity=K_fracture,
)

X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

Neumann = ComPASS.NeumannBC()
Neumann.molar_flux[:] = Qm
face_centers = simulation.face_centers()
where = (
    (np.abs(face_centers[:, 0]) < 0.25 * dx)
    & (np.abs(face_centers[:, 1]) <= 0.2 * Ly)
    & (np.abs(face_centers[:, 2]) <= 0.25 * dz)
)
left_fracture_edges = simulation.find_fracture_edges(where)
simulation.set_Neumann_fracture_edges(left_fracture_edges, Neumann)

simulation.standard_loop(
    initial_timestep=10,
    final_time=2 * day,
    output_period=4 * hour,
)
