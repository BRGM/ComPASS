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


def test_simple_tetmesh():

    omega_reservoir = 0.15  # reservoir porosity
    k_reservoir = 1e-12  # reservoir permeability in m^2
    K_reservoir = 2  # bulk thermal conductivity in W/m/K

    gridshape = (1, 1, 1)
    gridextent = (1e3, 1e3, 1e3)
    vertices, tets = GT.grid2tets(gridshape, gridextent)
    tets = [MT.Tetrahedron(MT.idarray(tet)) for tet in tets]

    mesh = MT.TetMesh.create(vertices, tets)

    vertices = mesh.vertices_array()
    zmax = vertices[:, -1].max()
    topnodes = np.nonzero(vertices[:, -1] == zmax)[0]

    simulation = ComPASS.load_physics("water2ph")

    ComPASS.set_output_directory_and_logfile(__file__)

    def select_dirichlet_nodes():
        print("Selecting", topnodes.shape[0], "top nodes.")
        on_top = np.zeros(mesh.nb_vertices, dtype=np.bool)
        on_top[topnodes] = True
        return on_top

    print("Gravity:", simulation.get_gravity())

    simulation.init(
        mesh=mesh,
        set_dirichlet_nodes=select_dirichlet_nodes,
        cell_porosity=omega_reservoir,
        cell_permeability=k_reservoir,
        cell_thermal_conductivity=K_reservoir,
    )

    def set_boundary_conditions():
        dirichlet = simulation.dirichlet_node_states()
        dirichlet.p[topnodes] = 1e5
        dirichlet.T[topnodes] = degC2K(30)
        dirichlet.context[:] = 2
        dirichlet.S[:] = [0, 1]
        dirichlet.C[:] = 1.0

    def set_initial_values():
        for state in [
            simulation.node_states(),
            simulation.fracture_states(),
            simulation.cell_states(),
        ]:
            state.context[:] = 2
            state.p[:] = 1e5
            state.T[:] = degC2K(30)
            state.S[:] = [0, 1]
            state.C[:] = 1.0

    set_boundary_conditions()
    set_initial_values()

    standard_loop(
        simulation,
        initial_timestep=1 * year,
        final_time=1e3 * year,
        output_period=1e2 * year,
    )


if __name__ == "__main__":
    test_simple_tetmesh()
