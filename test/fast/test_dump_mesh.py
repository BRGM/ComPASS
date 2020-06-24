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

    simulation = ComPASS.load_eos("water2ph")
    ComPASS.set_output_directory_and_logfile(__file__)
    simulation.lock_context(
        2
    )  # not useful here but just as a strict replacement of liquid_water eos

    Lx, Ly, Lz = 3000.0, 2000.0, 100.0
    Ox, Oy, Oz = -1500.0, -1000.0, -1600.0
    nx, ny, nz = 9, 6, 4

    grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(Lx, Ly, Lz), origin=(Ox, Oy, Oz),)

    def set_dirichlet():
        vertices = simulation.global_vertices()
        z = vertices[:, 2]
        zmax = z.max()
        return np.nonzero(z == zmax)[0]

    def set_flags():
        vertices = simulation.global_vertices()
        cell_centers = simulation.compute_global_cell_centers()
        face_centers = simulation.compute_global_face_centers()
        simulation.global_nodeflags()[:] = np.floor(vertices[:, 0])
        simulation.global_faceflags()[:] = np.floor(face_centers[:, 0])
        simulation.global_cellflags()[:] = np.floor(cell_centers[:, 0])

    simulation.init(
        mesh=grid,
        cell_porosity=0.1,
        cell_permeability=1e-12,
        cell_thermal_conductivity=2,
        set_global_flags=set_flags,
        set_dirichlet_nodes=set_dirichlet,
    )

    simulation.debug_utils_dump_mesh_info()
    ComPASS.debug_utils.extract_all_meshes()


if __name__ == "__main__":
    test_dump_mesh()
