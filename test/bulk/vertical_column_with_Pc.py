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
from ComPASS.timeloops import standard_loop, TimeStepManager


p0 = 1.0 * bar  # initial reservoir pressure
T0 = degC2K(
    20.0
)  # initial reservoir temperature - convert Celsius degrees to Kelvin degrees
bottom_heat_flux = 0.08  # W/m2
k_matrix = 1e-18  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K

H = 3000.0  # column height
nx, ny, nz = 1, 1, 300  # discretization

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(10.0, 10.0, H), origin=(-5, -5, -H),)

# Cappillary pressure
Pc0 = 2.0e5
Sg0 = 1.0 - 1.0e-2
Sl0 = 1.0 - Sg0
A = -Pc0 * np.log(Sl0) - (Pc0 / Sl0) * Sg0

s


def Pc(Sg):
    if Sg < Sg0:
        return -Pc0 * np.log(1.0 - Sg)
    return Pc0 * Sg / Sl0 + A


def dPcdS(Sg):
    if Sg < S0:
        return Pc0 / (1.0 - Sg)
    return Pc0 / Sl0


simulation.set_liquid_capillary_pressure(Pc, dPcdS)


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
)

states = simulation.all_states()
states.context[:] = simulation.Context.liquid
states.p[:] = p0
states.T[:] = T0
states.S[:] = 0
states.S[:, 1] = 1
states.C[:, :, :] = 0
states.C[1, 1, :] = 1  # liquid is pure water
# Top boundary conditions
top = simulation.vertices()[:, 2] >= 0
states.S[:, 0] = Sg0
states.S[:, 1] = 1 - Sg0
states.C[0, 0, :] = 1 - Cw0
states.C[0, 1, :] = Cw0
states.C[1, 0, :] = 0.001
states.C[1, 1, :] = 1 - 0.001
# Liquid composition should be consistent with Henry law at top boundary conditions...
simulation.reset_dirichlet_nodes(top)

final_time = 1e4 * year
output_period = 1e3 * year
standard_loop(
    simulation,
    final_time=final_time,
    time_step_manager=TimeStepManager(1, 100 * year),
    output_period=output_period,
)
