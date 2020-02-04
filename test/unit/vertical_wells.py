#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
from ComPASS.utils.wells import create_vertical_well
import ComPASS.mpi as mpi


simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(6, 6, 2), extent=(4, 4, 2), origin=(0, 0, 0),)


def make_wells():
    wells = []
    nx, ny, _ = grid.shape
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            well = create_vertical_well(simulation, (i, j), 0.1)  # dummy radius
            well.operate_on_flowrate = 1, 0.0  # dummy values
            well.produce()
            # well.inject(0)
            wells.append(well)
    for wk, well in enumerate(wells):
        well.id = wk
    return wells


simulation.init(
    mesh=grid,
    wells=make_wells,
    cell_porosity=0.1,  # dummy value
    cell_permeability=1e-12,  # dummy value
    cell_thermal_conductivity=2,  # dummy value
)


for wk, well in enumerate(simulation.producers_data()):
    print(f"well {wk} on proc {mpi.proc_rank} has id {well.id}")

print(f"Number of injectors on proc {mpi.proc_rank}: {simulation.number_of_own_injectors()} / {simulation.nb_injectors()}")
print(f"Number of producers on proc {mpi.proc_rank}: {simulation.number_of_own_producers()} / {simulation.nb_producers()}")

nx, ny, _ = grid.shape
nb_wells = (nx-2) * (ny - 2)
assert all([0 <= well.id <= nb_wells for well in simulation.producers_data()])
