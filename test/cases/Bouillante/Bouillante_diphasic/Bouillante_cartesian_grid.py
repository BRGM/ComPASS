#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Cartesian grid, box similar in size with Bouillante
# Dirichlet at the top bondary with Liquid in half and gas in the other half
# Imposed molar and heat flux at the bottom boundary
import os
import ComPASS
import numpy as np
import MeshTools as MT
import MeshTools.CGALWrappers as CGAL
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.mpi
from ComPASS.timestep_management import TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton


omega_reservoir = 0.35  # reservoir porosity  .1
k_reservoir = 1e-12  # reservoir permeability in m^2
cell_thermal_cond = 2.0  # reservoir thermal conductivity
p0 = 0.1 * MPa  # top Pressure
T0 = degC2K(20)  # top Temperature
Tbot_input = 550.0  # input flux Temperature
p_bot = 5.0e8
qbin = 1.0
CpRoche = 2.0e6
gravity = 10.0
liq_molar_fraction_qbin_air = 1.0e-1
liq_molar_fraction_qbin = np.array(
    [liq_molar_fraction_qbin_air, 1.0 - liq_molar_fraction_qbin_air]
)
geotherm = 0.007


Lx, Ly, Lz = 15.0e3, 16.0e3, 11125.0
Ox, Oy, Oz = 625.0e3, 1778.0e3, -10.0e3
nx, ny, nz = 15, 16, 50
Topz = Oz + Lz

x_source, y_source, radius = (
    Ox + Lx / 2.0,
    Oy + Ly / 2.0,
    3.0e3,
)  # coordinate and radius of the source center

gas_context = 1
liquid_context = 2
diphasic_context = 3
bot_flag = 4


simulation = ComPASS.load_eos("diphasic")
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
        vertices = np.rec.array(simulation.global_vertices())
        return vertices[:, 2] >= Topz

    def set_top_flags(pts, flags):
        flags[:] = 0
        for i, pt in enumerate(pts):
            if pt[2] >= Topz:
                if pt[0] - Ox < Lx / 2.0:
                    flags[i] = liquid_context
                else:
                    flags[i] = gas_context

    def set_global_flags():
        vertices = np.rec.array(simulation.global_vertices())
        nodeflags = simulation.global_nodeflags()
        set_top_flags(vertices, nodeflags)

        face_centers = np.rec.array(simulation.compute_global_face_centers())
        faceflags = simulation.global_faceflags()
        faceflags[:] = 0
        faceflags[face_centers[:, 2] <= Oz] = bot_flag


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

sys.stdout.write("Maillage distribue" + "\n")
sys.stdout.flush()


def to_array(xyz):
    array = np.zeros((xyz.shape[0], 3))
    for i, elt in enumerate(xyz):
        array[i] = np.fromiter(elt, dtype=np.float)
    return array


def lininterp(depths, top, gradient):
    return top + (gradient) * (depths)


def inside_heat_source(pts):
    return (pts[0] - x_source) ** 2 + (pts[1] - y_source) ** 2 < radius**2


def molar_fraction_balance(Pg, T):
    Hur = 0.5
    Cwg = Hur * simulation.Psat(T) / Pg
    Cag = 1.0 - Cwg
    Ha = 1.0e8
    Cal = Cag * Pg / Ha
    Cwl = 1.0 - Cal
    return [[Cag, Cwg], [Cal, Cwl]]


def set_Dirichlet_state(state):
    node_flags = simulation.nodeflags()
    # liq
    state.context[node_flags == liquid_context] = liquid_context
    state.p[node_flags == liquid_context] = p0
    state.T[node_flags == liquid_context] = T0
    state.S[node_flags == liquid_context] = [0, 1]
    state.C[node_flags == liquid_context] = [[1.0, 0.0], [1.0e-5, 1.0 - 1.0e-5]]
    # monophasic gas
    state.context[node_flags == gas_context] = gas_context
    state.p[node_flags == gas_context] = p0
    state.T[node_flags == gas_context] = T0
    state.S[node_flags == gas_context] = [1, 0]  # gas, liquid
    state.C[node_flags == gas_context] = molar_fraction_balance(p0, T0)


def set_states(state, depths):
    state.context[:] = liquid_context
    state.p[:] = lininterp(depths, p0, 900.0 * gravity)
    state.T[:] = lininterp(depths, T0, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = [[1.0, 0.0], [1.0e-5, 1.0 - 1.0e-5]]


def set_variable_initial_bc_values():
    set_states(simulation.node_states(), Topz - simulation.vertices()[:, 2])
    set_states(simulation.cell_states(), Topz - simulation.compute_cell_centers()[:, 2])
    set_Dirichlet_state(simulation.dirichlet_node_states())


# in some part of the bottom, molar flux and energy flux.
def set_variable_boundary_heat_flux():
    face_centers = np.rec.array(simulation.face_centers())
    face_flags = simulation.faceflags()

    neumann_heat_faces = np.zeros(len(face_flags), dtype=bool)
    neumann_faces = np.zeros(len(face_flags), dtype=bool)
    for i, pts in enumerate(face_centers):
        if face_flags[i] == bot_flag and inside_heat_source(pts):
            neumann_heat_faces[i] = True
        elif face_flags[i] == bot_flag:
            neumann_faces[i] = True
    # heat source
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = liq_molar_fraction_qbin[:] * qbin
    Neumann.heat_flux = qbin * simulation.liquid_molar_enthalpy(
        p_bot, Tbot_input, liq_molar_fraction_qbin
    )
    simulation.set_Neumann_faces(neumann_heat_faces, Neumann)
    # outside heat source
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = 0.0
    Neumann.heat_flux = simulation.liquid_molar_enthalpy(p_bot, 400.0, [0.0, 1.0])
    simulation.set_Neumann_faces(neumann_faces, Neumann)


sys.stdout.write("set initial and BC" + "\n")
set_variable_initial_bc_values()
sys.stdout.write("set Neumann BC" + "\n")
set_variable_boundary_heat_flux()
sys.stdout.flush()


init_dt = 1.0e-3
max_dt = 0.7 * year
final_time = 100.0 * year
output_period = 0.01 * final_time


lsolver = linear_solver(simulation, from_options=True)
newton = Newton(simulation, 1e-5, 8, lsolver)

standard_loop(
    simulation,
    final_time=final_time,
    output_every=20,
    time_step_manager=TimeStepManager(
        initial_timestep=init_dt,
        maximum_timestep=max_dt,
    ),
    newton=newton,
    nitermax=50,
    # output_period = output_period, specific_outputs=[1. * day], output_every=20,
)
