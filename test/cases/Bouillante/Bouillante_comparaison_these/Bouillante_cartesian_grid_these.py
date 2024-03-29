#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box similar in depth with Bouillante
# Imposed molar and heat flux at the bottom boundary

import ComPASS
import numpy as np
import MeshTools as MT
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.messages
from ComPASS.newton import Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.timestep_management import TimeStepManager
from ComPASS.utils.grid import on_zmin

Lz = 4000.0
nz = 500
dz = Lz / nz
Lx = Ly = 2 * dz
Ox, Oy, Oz = 0.0, 0.0, -3000.0
nx = ny = 2
Topz = Oz + Lz

omega_reservoir = 0.35  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2, 1D = 10^-12 m^
cell_thermal_cond = 3.0  # reservoir thermal conductivity
ptop = 0.1 * MPa  # top Pressure
Ttop = 300.0  # top Temperature
pbot = 2.5e7
Tbot = 400.0  # bottom Temperature
Tbot_input = 550.0  # bottom Temperature input
qbin = 2.0e-2  # 329./Lx/Ly
CpRoche = 2.0e6
gravity = 10.0
liq_molar_fraction_qbin_air = 0.0
liq_molar_fraction_qbin = np.array(
    [liq_molar_fraction_qbin_air, 1.0 - liq_molar_fraction_qbin_air]
)

geotherm = (Tbot - Ttop) / Lz

bot_flag = 4

simulation = ComPASS.load_physics("diphasic")
simulation.set_gravity(gravity)
simulation.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(
        shape=(nx, ny, nz),
        extent=(Lx, Ly, Lz),
        origin=(Ox, Oy, Oz),
    )

    def Dirichlet_node():
        return simulation.top_boundary(grid)()

    def set_global_flags():
        fc = np.rec.array(simulation.compute_global_face_centers())
        faceflags = simulation.global_faceflags()
        faceflags[:] = 0
        faceflags[on_zmin(grid)(fc)] = bot_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = (
        Dirichlet_node
    ) = omega_reservoir = k_reservoir = cell_thermal_cond = set_global_flags = None

simulation.init(
    mesh=grid,
    set_dirichlet_nodes=Dirichlet_node,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=cell_thermal_cond,
    set_global_flags=set_global_flags,
)


def lininterp(depths, top, gradient):
    assert np.all(depths) >= 0, "depths should be positif"
    return top + (gradient) * (depths)


def set_states(state, depths):
    state.context[:] = simulation.Context.liquid
    state.p[:] = lininterp(depths, ptop, gravity * 1000.0)
    state.T[:] = lininterp(depths, Ttop, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = [[1.0, 0.0], [0.0, 1.0]]


def set_variable_initial_bc_values():
    set_states(simulation.node_states(), Topz - simulation.vertices()[:, 2])
    set_states(simulation.cell_states(), Topz - simulation.compute_cell_centers()[:, 2])
    Xtop = simulation.build_state(simulation.Context.gas, p=ptop, T=Ttop)
    simulation.dirichlet_node_states().set(Xtop)


# at the bottom, molar and heat flux.
def set_variable_boundary_heat_flux():
    face_flags = simulation.faceflags()

    neumann_heat_faces = np.zeros(len(face_flags), dtype=bool)
    neumann_heat_faces[face_flags == bot_flag] = True

    # heat source
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = liq_molar_fraction_qbin[:] * qbin
    Neumann.heat_flux = qbin * simulation.liquid_molar_enthalpy(
        pbot, Tbot_input, liq_molar_fraction_qbin
    )
    simulation.set_Neumann_faces(neumann_heat_faces, Neumann)


set_variable_initial_bc_values()
set_variable_boundary_heat_flux()

# context = SimulationContext()
# context.abort_on_ksp_failure = False
# context.dump_system_on_ksp_failure = False
# context.abort_on_newton_failure = False

timestep = TimeStepManager(
    initial_timestep=1000.0,
    minimum_timestep=1e-8,
    maximum_timestep=10.0 * year,
    increase_factor=1.3,
    decrease_factor=0.2,
)

final_time = 300.0 * year
output_period = 0.01 * final_time

# Construct the linear solver and newton objects outside the time loop
# to set their parameters. Here direct solving is activated
lsolver = linear_solver(simulation, direct=False)
newton = Newton(simulation, 1e-8, 8, lsolver)

current_time = standard_loop(
    simulation,
    final_time=final_time,
    newton=newton,
    time_step_manager=timestep,
    # output_period=output_period,
    # specific_outputs=[1. * day], output_every=20,
)
