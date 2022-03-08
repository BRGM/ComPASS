# -*- coding: utf-8 -*-
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
import vtkwriters as vtkw

import gmsh_reader


ComPASS.load_eos("diphasic")

filename = "case-tests/andra/small_gallery.msh"

nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)


# grep 2d element
face = [elt for elt, tag in elements if type(elt) in (MT.Triangle, MT.Quad)]

face_tag = [tag for elt, tag in elements if type(elt) in (MT.Triangle, MT.Quad)]
face_tag = np.array([tag[0] for tag in face_tag])


# grep 3d element
cell = [
    elt
    for elt, tag in elements
    if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)
]

cell_tag = [
    tag
    for elt, tag in elements
    if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)
]
cell_tag = np.array([tag[0] for tag in cell_tag])


mesh = MT.HybridMesh.create(nodes, cell)


r = 2.625
y_min = -10
y_max = 10
z_min = -10
z_max = 10

gallery_flag = 1
left_flag = 2
right_flag = 2
bottom_flag = 2
top_flag = 2

eps = 1e-6
gallery_node = [y * y + z * z <= r * r + eps for x, y, z in nodes]


left_node = [abs(y - y_min) <= eps for x, y, z in nodes]
right_node = [abs(y - y_max) <= eps for x, y, z in nodes]

bottom_node = [abs(z - z_min) <= eps for x, y, z in nodes]
top_node = [abs(z - z_max) <= eps for x, y, z in nodes]


def set_physical_flags():
    cellflags = ComPASS.global_cellflags()
    cellflags[:] = cell_tag

    faceflags = ComPASS.global_faceflags()
    faceflags[:] = 0

    nodeflags = ComPASS.global_nodeflags()
    nodeflags[:] = 0

    nodeflags[gallery_node] = gallery_flag
    #    nodeflags[left_node] = left_flag
    #    nodeflags[right_node] = right_flag
    nodeflags[bottom_node] = bottom_flag
    nodeflags[top_node] = top_flag


def select_dirichlet_nodes():
    nodeflags = ComPASS.global_nodeflags()
    return nodeflags != 0


def select_global_rocktype():
    cellflags = ComPASS.global_cellflags()
    rocktype = 1 + ((cellflags - 1) % 2)

    cellrocktype = ComPASS.global_cellrocktype().reshape((-1, 2))
    cellrocktype[:] = np.stack((rocktype, rocktype), axis=-1)


ComPASS.set_output_directory_and_logfile(__file__)


ComPASS.init(
    mesh=mesh,
    set_global_rocktype=select_global_rocktype,
    set_global_flags=set_physical_flags,
    set_dirichlet_nodes=select_dirichlet_nodes,
)
# set_initial_values()
# set_boundary_conditions()


standard_loop(initial_timestep=100000, output_every=1, final_time=200 * year)
