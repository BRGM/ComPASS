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
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton
from ComPASS.timestep_management import TimeStepManager
import ComPASS.io.mesh as io

p0 = 4.0e6  # 40.0 * bar
patm = 1.0 * bar
T0 = 303.0
k_matrix = 5e-20  # column permeability in m^2 (low permeability -> bigger time steps)
phi_matrix = 0.15  # column porosity
K_matrix = 2.0  # bulk thermal conductivity in W/m/K
gravity = 9.81
rock_internal_energy = 2.0e6  # volumetric heat capacity J/K/m^3

H = 20.0  # column height
nx, ny, nz = 1, 1, 200  # discretization
ds = H / nz

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("diphasic")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(rock_internal_energy)

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(ds, ds, H),
    origin=(-0.5 * ds, -0.5 * ds, 0.0),
)

bottom_flag = 3
top_flag = 2


def set_physical_flags():
    nodeflags = simulation.global_nodeflags()
    nodeflags[:] = 0

    nodeflags[simulation.bottom_boundary(grid)()] = bottom_flag
    nodeflags[simulation.top_boundary(grid)()] = top_flag


def select_dirichlet_nodes():
    return simulation.top_boundary(grid)() | simulation.bottom_boundary(grid)()


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    set_global_flags=set_physical_flags,
    set_dirichlet_nodes=select_dirichlet_nodes,
)


# Set petrophyics functions after initialization
simulation.set_vanGenuchten_capillary_pressure()  # contains ox and cc
from data.van_genuchten_kr import kr_functions

simulation.set_kr_functions(kr_functions)

X0 = simulation.build_state(simulation.Context.liquid, p=p0, T=T0)
# diphasic states will call phase pressure function
# set the correct capillary pressure before using this function
Xdiph = simulation.build_state(simulation.Context.diphasic, p=patm, T=T0, Sg=0.35)


simulation.all_states().set(X0)

node_flags = simulation.nodeflags()
simulation.dirichlet_node_states().set(node_flags == bottom_flag, X0)
simulation.dirichlet_node_states().set(node_flags == top_flag, Xdiph)


def export_initial_states():
    node_states = simulation.node_states()
    cell_states = simulation.cell_states()
    petrophysics = simulation.petrophysics()

    pointdata = {
        "dirichlet pressure": simulation.pressure_dirichlet_values(),
        "dirichlet temperature": K2degC(simulation.temperature_dirichlet_values()),
        "initial pressure": node_states.p,
        "initial temperature": K2degC(node_states.T),
        "initial gas saturation": node_states.S[:, 0],
    }
    celldata = {
        "initial pressure": cell_states.p,
        "initial gas saturation": cell_states.S[:, 0],
        "initial temperature": K2degC(cell_states.T),
        "phi": petrophysics.cell_porosity,
    }
    io.write_mesh(
        simulation, "initial_states_andra.vtu", pointdata=pointdata, celldata=celldata
    )


# export_initial_states()

# linear solver tolerance must be smaller than the Newton tolerance
lsolver = linear_solver(simulation, tolerance=1e-8, legacy=False, from_options=True)
newton = Newton(simulation, 1e-7, 20, lsolver)
tsm = TimeStepManager(
    initial_timestep=20.0 * day,
    increase_factor=1.2,
    decrease_factor=0.5,
    minimum_timestep=0.1,
    maximum_timestep=10.0 * year,
)

run_loop = lambda final_time, no_output=True: simulation.standard_loop(
    initial_time=0,
    final_time=final_time,
    newton=newton,
    time_step_manager=tsm,
    no_output=no_output,
    output_every=10,
)

run_loop(
    1000 * year,
)
# simulation.postprocess()
