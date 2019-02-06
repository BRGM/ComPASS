#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS
import ComPASS.debug_utils

def test_dump_mesh():

    ComPASS.load_eos('liquid_water')
    ComPASS.set_output_directory_and_logfile(__file__)

    Lx, Ly, Lz = 3000., 2000., 100.
    Ox, Oy, Oz = -1500., -1000., -1600.
    nx, ny, nz = 9, 6, 4

    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
        origin = (Ox, Oy, Oz),
    )

    def set_dirichlet():
        vertices = ComPASS.global_vertices()
        z = vertices[:, 2]
        zmax = z.max()
        return np.nonzero(z == zmax)[0]

    def set_flags():
        vertices = ComPASS.global_vertices()
        cell_centers = ComPASS.compute_global_cell_centers()
        face_centers = ComPASS.compute_global_face_centers()
        ComPASS.global_nodeflags()[:] = np.floor(vertices[:, 0])
        ComPASS.global_faceflags()[:] = np.floor(face_centers[:, 0])
        ComPASS.global_cellflags()[:] = np.floor(cell_centers[:, 0])


    ComPASS.init(
        mesh = grid,
        cell_porosity = 0.1,
        cell_permeability = 1e-12,
        cell_thermal_conductivity = 2,
        set_global_flags = set_flags,
        set_dirichlet_nodes = set_dirichlet,
    )

    ComPASS.debug_utils_dump_mesh_info()
    ComPASS.debug_utils.extract_all_meshes()

if __name__ == '__main__':
    test_dump_mesh()
