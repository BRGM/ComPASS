#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS


def test_access_vertices():

    ComPASS.set_output_directory_and_logfile(__file__)
    simulation = ComPASS.load_physics("water2ph")

    nx, ny, nz = 3, 2, 1

    grid = ComPASS.Grid(
        shape=(nx, ny, nz),
    )

    def dirichlet_nodes():
        vertices = simulation.global_vertices()
        newz = 10.0
        z = vertices[:, 2]
        z[z != 0] = newz
        assert np.all(
            np.logical_or(
                simulation.global_vertices()[:, 2] == 0,
                simulation.global_vertices()[:, 2] == newz,
            )
        )

    simulation.init(
        mesh=grid,
        set_dirichlet_nodes=dirichlet_nodes,
        cell_porosity=0.1,
        cell_permeability=1.0e-12,
        cell_thermal_conductivity=2.0,
    )


if __name__ == "__main__":
    test_access_vertices()
