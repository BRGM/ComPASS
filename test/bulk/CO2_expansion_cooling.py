#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import numpy as np
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.timestep_management import TimeStepManager


Lz = 4000.0
nz = 400
dz = Lz / nz
Lx = Ly = dz
Ox, Oy, Oz = 0.0, 0.0, -3000.0
nx = ny = 1
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12 * np.eye(3)  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity
pall = 1.0e6  # top Pressure
pbot_dir = 1.0e5  # top Pressure
Tall = 700.0  # top Temperature
Tbot_dir = 700.0  # bottom Temperature
gravity = 0.0

simulation = ComPASS.load_physics("diphasicCO2")
simulation.set_gravity(gravity)
ComPASS.set_output_directory_and_logfile(__file__)


grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)


def Dirichlet_node():
    vertices = np.rec.array(simulation.global_vertices())
    return vertices[:, 2] <= Oz


simulation.init(
    mesh=grid,
    set_dirichlet_nodes=Dirichlet_node,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)

X0 = simulation.build_state(simulation.Context.gas, p=pall, T=Tall, Cag=0.8)
simulation.all_states().set(X0)
XDir = simulation.build_state(simulation.Context.gas, p=pbot_dir, T=Tbot_dir, Cag=0.8)
simulation.dirichlet_node_states().set(XDir)

mol_enth_ex = simulation.cpp_liquid_molar_enthalpy(pall, Tall, [0.8, 0.2])

timestep = TimeStepManager(
    initial_timestep=100.0,
    minimum_timestep=1e-1,
    maximum_timestep=50.0 * year,
    increase_factor=1.2,
    decrease_factor=0.2,
)

final_time = 1000.0 * year
output_period = 0.01 * final_time


current_time = standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=timestep,
    output_period=output_period,
)

print("time after the time loop", current_time / year, "years")
