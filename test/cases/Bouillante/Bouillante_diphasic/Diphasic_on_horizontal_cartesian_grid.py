#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#
# Test case with 2 comp, 2 phases
# Neuman flux of liquid (water + air ) from the left part of the domain
# Dirichlet on the right part of the domain

import ComPASS
import numpy as np
import MeshTools as MT
import MeshTools.CGALWrappers as CGAL
import sys
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import ComPASS.mpi


omega_reservoir = 0.35  # reservoir porosity  .1
k_reservoir = 1e-12  # reservoir permeability in m^2
cell_thermal_cond = 2.0  # reservoir thermal conductivity
p0 = 0.1 * MPa  # top Pressure
T0 = degC2K(20)  # top Temperature
CpRoche = 2.0e6
gravity = 10.0
qbin = 1.0
liq_molar_fraction_qbin_air = 1.0e-1
liq_molar_fraction_qbin = np.array(
    [liq_molar_fraction_qbin_air, 1.0 - liq_molar_fraction_qbin_air]
)
bottom_heat_flux = 0.1
geotherm = 0.02  # bottom_heat_flux / cell_thermal_cond # gradient Temperature/m


Lx, Ly, Lz = 1.0e3, 100.0, 100.0
Ox, Oy, Oz = 0.0, 0.0, 0.0
nx, ny, nz = 20, 3, 3
Topz = Oz + Lz

gas_context = 1
liquid_context = 2
diphasic_context = 3


left_flag = 2
right_flag = 3


ComPASS.load_eos("diphasic")
ComPASS.set_gravity(gravity)
ComPASS.set_rock_volumetric_heat_capacity(CpRoche)
ComPASS.set_output_directory_and_logfile(__file__)

if ComPASS.mpi.is_on_master_proc:

    grid = ComPASS.Grid(
        shape=(nx, ny, nz),
        extent=(Lx, Ly, Lz),
        origin=(Ox, Oy, Oz),
    )

    def Dirichlet_node():
        vertices = np.rec.array(ComPASS.global_vertices())
        return vertices[:, 0] >= Ox + Lx

    def set_global_flags():
        vertices = np.rec.array(ComPASS.global_vertices())
        nodeflags = ComPASS.global_nodeflags()
        nodeflags[:] = 0
        nodeflags[(vertices[:, 0] >= Ox + Lx)] = right_flag
        nodeflags[(vertices[:, 0] <= Ox)] = left_flag

        face_centers = np.rec.array(ComPASS.compute_global_face_centers())
        faceflags = ComPASS.global_faceflags()
        faceflags[:] = 0
        faceflags[(face_centers[:, 0] >= Ox + Lx)] = right_flag
        faceflags[(face_centers[:, 0] <= Ox)] = left_flag


if not ComPASS.mpi.is_on_master_proc:
    grid = (
        Dirichlet_node
    ) = fracture_faces = set_global_flags = select_global_rocktype = None

ComPASS.init(
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


def set_states(state, depths):
    # context 1:Ig; 2:Il; 3:Ig+Il
    state.context[:] = liquid_context
    state.p[:] = lininterp(depths, p0, 800 * gravity)
    state.T[:] = lininterp(depths, T0, geotherm)
    state.S[:] = [0, 1]
    state.C[:] = [[1.0, 0.0], [1.0e-3, 1.0 - 1.0e-3]]


def set_variable_initial_bc_values():
    set_states(ComPASS.node_states(), abs(ComPASS.vertices()[:, 2] - Topz))
    set_states(ComPASS.cell_states(), abs(ComPASS.compute_cell_centers()[:, 2] - Topz))
    set_states(ComPASS.dirichlet_node_states(), abs(ComPASS.vertices()[:, 2] - Topz))


# molar and energy flux on the left part
def set_variable_boundary_heat_flux():
    face_centers = np.rec.array(ComPASS.face_centers())
    face_flags = ComPASS.faceflags()

    neumann_faces = np.zeros(face_flags.size, dtype=bool)  # face_centers.size
    neumann_faces[(face_flags == left_flag)] = True
    Neumann = ComPASS.NeumannBC()
    Neumann.molar_flux[:] = (
        liq_molar_fraction_qbin[:] * qbin
    )  # two components (air, water)
    Neumann.heat_flux = qbin * ComPASS.liquid_molar_enthalpy(
        p0, T0 + 20.0, liq_molar_fraction_qbin
    )
    ComPASS.set_Neumann_faces(neumann_faces, Neumann)


sys.stdout.write("set initial and BC" + "\n")
set_variable_initial_bc_values()
sys.stdout.write("set Neumann BC" + "\n")
set_variable_boundary_heat_flux()
sys.stdout.flush()


init_dt = 0.01 * hour
final_time = 1000.0 * year
output_period = 0.01 * final_time
ComPASS.set_maximum_timestep(0.7 * year)


def ma_fonction(it_timestep, time):
    if it_timestep > 300:
        ComPASS.mpi.abort()


standard_loop(
    initial_timestep=init_dt,
    final_time=final_time,
    output_every=20,
    iteration_callbacks=[ma_fonction],
    # output_period = output_period, specific_outputs=[1. * day], output_every=20,
)
