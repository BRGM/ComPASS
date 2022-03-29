#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import MeshTools as MT
import MeshTools.GridTools as GT
import vtkwriters as vtkw

ComPASS.load_eos("diphasic")

x, y, z = 15, 1, 15
nx, ny, nz = 30, 1, 60


liquid_water_context = 2
water_density = 1.0e3
atmospheric_pressure = 1e5
temperature = 303
gravity = 10

cct_z = 1.5


grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(x, y, z),
    origin=(0, 0, 0),
)
# vertices, hexs = GT.grid2hexs((nx, ny, nz), (x, y, z))
# hexs = np.asarray(hexs, dtype=MT.idtype())
# mesh = MT.HexMesh.make(vertices, hexs)

vertices, hexs = GT.grid2hexs((nx, ny, nz), (x, y, z))
hexs = np.ascontiguousarray(hexs, dtype=MT.idtype())
mesh = MT.HexMesh.make(vertices, hexs)
vertices = MT.as_coordinate_array(mesh.vertices)


def select_global_flags():
    select_global_nodeflags()
    select_global_cellflags()


def select_global_nodeflags():
    on_zmin = lambda pts: pts[:, 2] == grid.origin[2]
    on_zmax = lambda pts: pts[:, 2] == grid.origin[2] + grid.extent[2]

    vertices = ComPASS.global_vertices().view(np.double).reshape((-1, 3))

    flags = ComPASS.global_nodeflags()
    flags[:] = 0
    flags[on_zmin(vertices)] = 1
    flags[on_zmax(vertices)] = 2


def select_global_cellflags():
    cell_centers = ComPASS.compute_global_cell_centers()
    CCT = cell_centers[:, 2] <= cct_z
    COX = cell_centers[:, 2] > cct_z

    flags = ComPASS.global_cellflags()
    flags[CCT] = 2
    flags[COX] = 1


# def select_dirichlet_nodes():
#    print('Selecting', top_nodes.shape[0], 'top nodes.')
#
#    dirichlet = np.zeros(mesh.nb_nodes(), dtype=np.bool)
#    dirichlet[top_nodes] = True
#    dirichlet[bottom_nodes] = True
#    return dirichlet


def select_dirichlet_nodes(grid):
    on_xmin = lambda pts: pts[:, 0] == grid.origin[0]
    on_xmax = lambda pts: pts[:, 0] == grid.origin[0] + grid.extent[0]
    on_ymin = lambda pts: pts[:, 1] == grid.origin[1]
    on_ymax = lambda pts: pts[:, 1] == grid.origin[1] + grid.extent[1]
    on_zmin = lambda pts: pts[:, 2] == grid.origin[2]
    on_zmax = lambda pts: pts[:, 2] == grid.origin[2] + grid.extent[2]

    def select():
        vertices = ComPASS.global_vertices().view(np.double).reshape((-1, 3))
        return on_zmin(vertices) | on_zmax(vertices)

    return select


def select_global_rocktype():
    set_cell_rocktype()
    set_frac_rocktype()


def set_frac_rocktype():
    rocktype = ComPASS.global_fracrocktype().reshape((-1, 2))
    rocktype[:] = 1


def set_cell_rocktype():
    cell_centers = ComPASS.compute_global_cell_centers()
    CCT = cell_centers[:, 2] <= cct_z
    COX = cell_centers[:, 2] > cct_z

    rocktype = ComPASS.global_cellrocktype().reshape((-1, 2))
    rocktype[CCT] = 2
    rocktype[COX] = 1


def pressure(reference_pressure, density, gravity, z):
    return reference_pressure - density * gravity * z


def pure_phase_state(context, pressure, temperature, state):
    state.context[:] = context
    state.p[:] = pressure
    state.T[:] = temperature
    state.S[:] = [0, 1]
    state.C[:] = 1.0


def set_initial_values():
    z = ComPASS.vertices().view(np.double).reshape((-1, 3))[:, -1]
    pure_phase_state(
        liquid_water_context,
        pressure(atmospheric_pressure, water_density, gravity, z),
        temperature,
        ComPASS.node_states(),
    )

    z = ComPASS.compute_face_centers()[ComPASS.frac_face_id() - 1, -1]
    pure_phase_state(
        liquid_water_context,
        pressure(atmospheric_pressure, water_density, gravity, z),
        temperature,
        ComPASS.fracture_states(),
    )

    z = ComPASS.compute_cell_centers()[:, -1]
    pure_phase_state(
        liquid_water_context,
        pressure(atmospheric_pressure, water_density, gravity, z),
        temperature,
        ComPASS.cell_states(),
    )


def set_boundary_conditions():
    nodes = ComPASS.vertices().view(np.double).reshape((-1, 3))
    dirichlet_node = ComPASS.dirichlet_node_states()

    top_nodes = np.nonzero(nodes[:, 2] == z)[0]
    liquid_water_pressure = pressure(atmospheric_pressure, water_density, gravity, z)
    liquid_water_temperature = temperature

    dirichlet_node.p[top_nodes] = liquid_water_pressure
    dirichlet_node.T[top_nodes] = liquid_water_temperature
    dirichlet_node.context[top_nodes] = liquid_water_context

    bottom_nodes = np.nonzero(nodes[:, 2] == 0)[0]
    liquid_water_pressure = pressure(atmospheric_pressure, water_density, gravity, 0)
    liquid_water_temperature = temperature

    dirichlet_node.p[bottom_nodes] = liquid_water_pressure
    dirichlet_node.T[bottom_nodes] = liquid_water_temperature
    dirichlet_node.context[bottom_nodes] = liquid_water_context

    dirichlet_node.context[:] = liquid_water_context
    dirichlet_node.S[:] = [0, 1]
    dirichlet_node.C[:] = 1


ComPASS.set_output_directory_and_logfile(__file__)


ComPASS.init(
    #    mesh = mesh,
    grid=grid,
    set_dirichlet_nodes=select_dirichlet_nodes(grid),
    set_global_rocktype=select_global_rocktype,
    set_global_flags=select_global_flags,
)
# set_initial_values()
# set_boundary_conditions()


standard_loop(initial_timestep=1000, output_every=1, final_time=200 * year)


# vtkw.write_vtu(
#    vtkw.vtu_doc(vertices, mesh.cellnodes),
#    'andra.vtu'
# )
