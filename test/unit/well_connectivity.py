#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import ComPASS.mpi as mpi

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

shape = nx, ny, nz = 4, 4, 6
grid = ComPASS.Grid(
    shape=shape,
    extent=shape,
    origin=(0, 0, 0),
)


def make_wells():
    wells = []
    for x in [1, 2]:
        well = simulation.create_vertical_well((x, 1), 0.1)
        well.operate_on_flowrate = 1.0, -np.inf
        well.produce()
        wells.append(well)
        well = simulation.create_vertical_well((x, 2), 0.1)
        well.operate_on_flowrate = 1.0, np.inf
        well.inject(293.13)
        wells.append(well)
    return wells


simulation.init(
    mesh=grid,
    wells=make_wells,
    cell_porosity=0.1,
    cell_permeability=1.0,
    cell_thermal_conductivity=1.0,
)

for well_type, info in [
    ("Producer", simulation.producers_information()),
    ("Injector", simulation.injectors_information()),
]:
    for k, well in enumerate(info):
        print(f"Proc {mpi.proc_rank}: {well_type} {k} well vertices:", well.vertices)
        print(
            f"Proc {mpi.proc_rank}: {well_type} {k} parent reservoir node:",
            well.parent_vertex,
        )
        print(
            f"Proc {mpi.proc_rank}: {well_type} {k} parent CSR offset:",
            well.parent_offset,
        )
        print(
            f"Proc {mpi.proc_rank}: {well_type} {k} parent well node:", well.parent_rank
        )
        # must always be true
        assert np.all(well.vertices[well.parent_rank] == well.parent_vertex)
