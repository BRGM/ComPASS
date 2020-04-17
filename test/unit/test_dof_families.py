#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS import mpi

shape = nx, ny, nz = 2, 2, 1  # discretization

simulation = ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape)

simulation.init(
    mesh = grid,
    cell_permeability = 1,
    cell_porosity = 0.5,
    cell_thermal_conductivity = 1,
)

mpi.master_print("\nNode families\n")

node_family = simulation.NumNodebyProc()

print(f"Family offsets on proc {mpi.proc_rank}:", node_family.offsets())
print(f"Family on proc {mpi.proc_rank}\n{node_family.as_array()}")
print(f"Nodes on proc {mpi.proc_rank} are distributed accross {node_family.number_of_domains()} domains.")
for dk, domain_nodes in enumerate(node_family):
    # print could be used to retrieve a two colums vector
    # print(domain_nodes.as_array())
    proc = np.unique(domain_nodes.proc)
    assert len(proc)==1 # A single proc per domain
    # print(f"Domain {dk} nodes of proc {mpi.proc_rank} are held by proc {proc[0]}.")
    print(f"Domain {dk} nodes of proc {mpi.proc_rank} have local id on proc {proc[0]}:",
        domain_nodes.local_id -1 # Fortran to C indexing !!!
    )
