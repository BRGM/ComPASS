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

Lx, Ly, Lz = 3000.0, 2000.0, 100.0
Ox, Oy, Oz = -1500.0, -1000.0, -1600.0
nx, ny, nz = 31, 21, 6

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

def fractures_factory(grid):
    def select_fractures():
        face_centers = simulation.compute_global_face_centers()
        dz = grid.extent[2] / grid.shape[2]
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = grid.origin[2] + 0.5 * grid.extent[2]
        return np.abs(face_centers[:, 2] - zfrac) < 0.25 * dz

    return select_fractures

def wells_factory(grid):
    def make_wells():
        interwell_distance = 1 * km
        Cx, Cy, Cz = grid_center(grid)
        producer = simulation.create_vertical_well((Cx - 0.5 * interwell_distance, Cy))
        producer.operate_on_flowrate = Qm, 1.0 * bar
        producer.produce()
        injector = simulation.create_vertical_well((Cx + 0.5 * interwell_distance, Cy))
        injector.operate_on_flowrate = -Qm, pres + 100.0 * MPa
        injector.inject(Tinjection)
        return (producer, injector)

    return make_wells

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    wells=wells_factory(grid),
    fracture_faces=fractures_factory(grid),
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

simulation.standard_loop(
    initial_timestep=day, final_time=2 * year, output_period=30 * day
)
