# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import MeshTools as MT
import vtkwriters as vtkw
import gmsh_reader

filename = "case-tests/andra/andra_gallery.msh"

nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)

# select heaxedra only
hexaedra = [elt for elt, tags in elements if type(elt) is MT.Hexahedron]
mesh = MT.HexMesh.create(nodes, hexaedra)

vertices = MT.as_coordinate_array(mesh.vertices)
cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, cellnodes), filename.replace(".msh", "_hexaedra_only.vtu")
)

# filter out 2D elements
elements_3D = [
    (elt, tags)
    for elt, tags in elements
    if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)
]

mesh = MT.HybridMesh.create(nodes, [elt for elt, tags in elements_3D])
physical = np.array([tags[0] for elt, tags in elements_3D])

offsets, cellnodes = mesh.cells_nodes_as_COC()

vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        nodes,
        np.array(offsets[1:], copy=False),  # no first zero offset for wtk
        np.array(cellnodes, copy=False),
        mesh.cells_vtk_ids(),
        celldata={"physical": physical},
    ),
    filename.replace(".msh", "_3D.vtu"),
)

# Indentify rear faces and export them as 2D mesh with the physical tag
fc = MT.as_coordinate_array(mesh.face_centers())
bf = MT.as_id_array(mesh.boundary_faces())
where = fc[bf, 0] == 0  # x coordinate == 0 ?
rear_faces = bf[where]
# retrieve rear_faces tag using their cells
bc = MT.as_id_array(mesh.boundary_cells())
rear_faces_tag = physical[bc[where]]

face_nodes = mesh.connectivity.faces.nodes
rear_elements = [face_nodes[fk] for fk in rear_faces]
rear_cellnodes = [np.array(elt) for elt in rear_elements]
rear_offsets = np.cumsum([a.shape[0] for a in rear_cellnodes])
rear_cellnodes = np.hstack(rear_cellnodes)
faces_vtk_ids = mesh.faces_vtk_ids()
rear_celltypes = np.array([faces_vtk_ids[fk] for fk in rear_faces])

vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        nodes,
        rear_offsets,
        rear_cellnodes,
        rear_celltypes,
        celldata={"physical": rear_faces_tag},
    ),
    filename.replace(".msh", "_rear_face.vtu"),
)
