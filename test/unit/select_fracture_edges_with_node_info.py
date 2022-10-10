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
import ComPASS.mpi as mpi
from ComPASS.utils.tags import retrieve_fracture_edges_with_node_tag


k = 1e-15  # dummy value
phi = 0.15  # dummy value
thermal_cond = 2.0  # dummy value

nx, ny, nz = (4,) * 3

simulation = ComPASS.load_physics("linear_water")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(nx, ny, nz), origin=(-0.5, -0.5, -0.5), extent=(1.0, 1.0, 1.0)
)


def set_node_flags():
    xyz = simulation.global_vertices()
    x, y, z = [xyz[:, j] for j in range(3)]
    flags = simulation.global_nodeflags()
    flags[(x == 0) & ((np.abs(y) == 0.5) | (np.abs(z) == 0.5))] = 1
    flags[(y == 0) & ((np.abs(x) == 0.5) | (np.abs(z) == 0.5))] = 2
    flags[(z == 0) & ((np.abs(x) == 0.5) | (np.abs(y) == 0.5))] = 3


def select_fractures():
    centers = simulation.compute_global_face_centers()
    xc, yc, zc = [centers[:, j] for j in range(3)]
    where = (xc == 0) | (yc == 0) | (zc == 0)
    print(np.sum(where), "fracture faces")
    return where


# The call to simulation.init will distribute the mesg

simulation.init(
    mesh=grid,
    cell_permeability=k,  # dummy value
    cell_porosity=phi,  # dummy value
    cell_thermal_conductivity=thermal_cond,  # dummy value
    fracture_permeability=k,  # dummy value
    fracture_porosity=phi,  # dummy value
    fracture_thermal_conductivity=thermal_cond,  # dummy value
    fracture_faces=select_fractures,
    set_global_flags=set_node_flags,
)

# The mesh is now distributed
assert simulation.mesh_is_local

retrieve_fracture_edges_with_node_tag(simulation, range(3), verbose=True)
