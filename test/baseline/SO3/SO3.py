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
pres = 80. * bar            # initial reservoir pressure
Tres = degC2K( 78. )        # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Qm = 300. * ton / hour      # production flowrate
Tinjection = degC2K( 33. )  # injection temperature - convert Celsius to Kelvin degrees
k_reservoir = 3E-12         # reservoir permeability in m^2
omega_reservoir = 0.2       # reservoir porosity
k_burden = 0                # reservoir permeability in m^2
omega_burden = 0.01         # reservoir porosity
reservoir_thickness = 40.   #
burden_thickness = 30.      #
top_model_depth = 1500      #
K_reservoir = 2             # bulk thermal conductivity in W/m/K
cp = 900.                   # J/K/kg
rhos = 2500.                # kg/m3
# fmt: on

Lx, Ly, Lz = 8000.0, 10000.0, 2 * burden_thickness + reservoir_thickness
Ox, Oy, Oz = -Lx / 2, -Ly / 2, -top_model_depth - Lz
nx, ny, nz = 80, 100, 10
bottom_reservoir = Oz + burden_thickness
top_reservoir = bottom_reservoir + reservoir_thickness

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)
simulation.set_gravity(0)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)


def make_wells():
    Cx, Cy, Cz = grid_center(grid)
    P = (Cx, Cy - 1 * km)
    I = (Cx, Cy + 1 * km)
    producer = simulation.create_vertical_well(
        P, zmin=bottom_reservoir, zmax=top_reservoir
    )
    producer.operate_on_flowrate = Qm, 1.0 * bar
    producer.produce()
    injector = simulation.create_vertical_well(
        I, zmin=bottom_reservoir, zmax=top_reservoir
    )
    injector.operate_on_flowrate = -Qm, pres + 100.0 * MPa
    injector.inject(Tinjection)
    return (producer, injector)


def set_permeability():
    xyz = simulation.compute_global_cell_centers()
    z = xyz[:, 2]
    k = np.zeros(xyz.shape[0], dtype="d")
    k[(z > bottom_reservoir) & (z < top_reservoir)] = k_reservoir
    return k


def set_porosity():
    xyz = simulation.compute_global_cell_centers()
    z = xyz[:, 2]
    omega = np.zeros(xyz.shape[0], dtype="d")
    omega[:] = omega_burden
    omega[(z > bottom_reservoir) & (z < top_reservoir)] = omega_reservoir
    return omega


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=simulation.vertical_boundaries(grid),
    # wells=make_wells,
    cell_permeability=set_permeability,
    cell_porosity=set_porosity,
    cell_thermal_conductivity=K_reservoir,
)

X0 = simulation.build_state(simulation.Context.liquid, p=pres, T=Tres)
simulation.all_states().set(X0)
simulation.dirichlet_node_states().set(X0)

simulation.standard_loop(
    initial_timestep=30 * day, final_time=30 * year, output_period=1 * year,
)

simulation.postprocess()
