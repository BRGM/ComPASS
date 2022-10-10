# -*- coding: utf-8 -*-

import numpy as np
import ComPASS
import ComPASS.mpi as mpi
from ComPASS.utils.salome import SalomeWrapper


ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_physics("water2ph")

sw = SalomeWrapper(simulation)


def create_well():
    nodes = sw.info.well.nodes
    z = sw.info.mesh.vertices[nodes, 2]
    print(z)
    print(nodes)
    print(nodes[np.argsort(z)[::-1]])
    well = simulation.create_single_branch_well(nodes[np.argsort(z)[::-1]])
    return [well]


simulation.init(
    mesh=sw.mesh,
    wells=create_well,
    cell_permeability=1e-12,
    cell_porosity=0.5,
    cell_thermal_conductivity=2,
    fracture_thermal_conductivity=2,
    fracture_faces=lambda: sw.info.bottom.faces,
    fracture_permeability=1e-16,
    fracture_porosity=0.1,
    set_global_flags=sw.flags_setter,
)

sw.rebuild_locally()
sw.info.to_vtu_block(f"block-proc{mpi.proc_rank:04d}")
sw.info.faces_to_multiblock(f"surfaces-proc{mpi.proc_rank:04d}")
