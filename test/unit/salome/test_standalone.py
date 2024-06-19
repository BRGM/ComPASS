# -*- coding: utf-8 -*-

import numpy as np
from ComPASS.utils import salome

# could be: salome.load(nodes_file="mynodes.txt", tets_file="mytetras.txt", groups_module="GROUPS"):
mesh, mesh_info = salome.load()

mesh_info.to_vtu_block("salome-block")
mesh_info.faces_to_multiblock("salome-surfaces")

mesh = mesh_info.mesh
flags = salome.ElementValues(
    np.zeros(mesh.nb_vertices, dtype=np.int32),
    np.zeros(mesh.nb_faces, dtype=np.int32),
    np.zeros(mesh.nb_cells, dtype=np.int32),
)
mesh_info.set_flags(flags)
masks = mesh_info.collect_masks()

# this is used to rebuild the mesh from local inforamtion
mesh_info.rebuild_from_flags(flags, masks, mesh=mesh)

mesh_info.to_vtu_block("rebuilt-block")
mesh_info.faces_to_multiblock("rebuilt-surfaces")
