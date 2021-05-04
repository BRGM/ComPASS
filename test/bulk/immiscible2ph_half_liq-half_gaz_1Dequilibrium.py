#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# use immiscible2ph EOS (air only in gaz phase, water only in liquid phase)
# Cartesian grid, 1D box with lx=1000m, 200 cells
# Homogeneous Neumann at all BC
# no gravity
# the left half is initiated with liquid,
#   the right half with linear Sg between 0 at Lx/2 and 1 at Lx

import ComPASS
import numpy as np
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.timestep_management import TimeStepManager
from ComPASS.mpi import master_print
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton


Lx = 1000.0
nx = 200
dx = Lx / nx
Lz = Ly = dx
Ox, Oy, Oz = 0.0, 0.0, 0.0
nz = ny = 1

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12 * np.eye(3)  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity : no thermal diffusion
Ptop = 1.0e5  # porous top Pressure
Tporous = 300.0  # porous Temperature (used also to init the freeflow nodes)
CpRoche = 2.0e6
pure_phase_molar_fraction = [[1.0, 0.0], [0.0, 1.0]]

simulation = ComPASS.load_eos("immiscible2ph")
ComPASS.set_output_directory_and_logfile(__file__)

simulation.set_rock_volumetric_heat_capacity(CpRoche)

liquid_context = simulation.Context.liquid
gas_context = simulation.Context.gas

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz))


if not ComPASS.mpi.is_on_master_proc:
    grid = omega_reservoir = k_reservoir = cell_thermal_cond = None

simulation.init(
    mesh=grid,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
)


def set_states(state, x):
    left = x <= Lx / 2.0
    right = np.logical_not(left)
    # liquid below 0
    state.context[left] = liquid_context
    state.p[left] = Ptop
    state.T[left] = Tporous
    state.S[left] = [0, 1]
    state.C[left] = pure_phase_molar_fraction
    # diphasic above 0
    state.context[right] = gas_context
    state.p[right] = Ptop
    state.T[right] = Tporous
    state.S[right] = [1, 0]
    state.C[right] = pure_phase_molar_fraction


def set_variable_initial_bc_values():
    set_states(simulation.node_states(), simulation.vertices()[:, 0])
    set_states(simulation.cell_states(), simulation.compute_cell_centers()[:, 0])


master_print("set initial and BC")
set_variable_initial_bc_values()

timestep = TimeStepManager(
    initial_timestep=100.0,
    minimum_timestep=1e-3,
    maximum_timestep=10.0 * year,
    increase_factor=1.2,
    decrease_factor=0.5,
)

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

final_time = 100.0 * year
output_period = 0.01 * final_time

current_time = simulation.standard_loop(
    final_time=final_time,
    time_step_manager=timestep,
    output_period=output_period,
    newton=newton,
)

print("time after the time loop", current_time / year, "years")
