#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import doublet_utils
from chessboard import chessboard


simulation = ComPASS.load_eos('water2ph')

pres = 20. * MPa                  # initial reservoir pressure
Tres = degC2K( 70. )              # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
Tinjection = degC2K( 30. )        # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour            # production flowrate
omega_reservoir = 0.15            # reservoir porosity
k_reservoir = 1E-12               # reservoir permeability in m^2
K_reservoir = 2                   # bulk thermal conductivity in W/m/K

Lx, Ly, Lz = 3000., 2000., 100.
Ox, Oy, Oz = -1500., -1000., -1600.
nx, ny, nz = 30, 20, 6

simulation.set_gravity(0)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (Ox, Oy, Oz),
)

def color_grid():
    colors = chessboard(
        grid.shape, (3, 2, 3),
        mpi.communicator().size,
        loops_order=(1, 2, 0,) # x first then z then y
    )
    return colors.ravel(order='F')

def make_wells():
    interwell_distance = 1 * km
    Cx, Cy, Cz = doublet_utils.center(grid)
    producer = doublet_utils.make_well(simulation, (Cx - 0.5 * interwell_distance, Cy))
    producer.operate_on_flowrate = Qm , 1. * bar
    producer.produce()
    injector = doublet_utils.make_well(simulation, (Cx + 0.5 * interwell_distance, Cy))
    injector.operate_on_flowrate = Qm, pres + 100. * MPa
    injector.inject(Tinjection)
    return (producer, injector)

ComPASS.set_output_directory_and_logfile(__file__)
simulation.init(
    mesh = grid,
    set_dirichlet_nodes = doublet_utils.select_boundary_factory(simulation, grid),
    wells = make_wells,
    cell_porosity = omega_reservoir,
    cell_permeability = k_reservoir,
    cell_thermal_conductivity = K_reservoir,
    mesh_parts = color_grid(),
)

doublet_utils.init_states(simulation, pres, Tres)

standard_loop(
    simulation,
    initial_timestep = 1 * day, final_time = 30 * year,
    output_period = year,
)
