#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import ComPASS.dump_wells as dw
import ComPASS.mpi as mpi


simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape=(6, 6, 2),
    extent=(4, 4, 2),
    origin=(0, 0, 0),
)


def make_wells():
    wells = []
    nx, ny, _ = grid.shape
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            well = simulation.create_vertical_well((i, j), 0.01)  # dummy radius
            well.operate_on_flowrate = 1, 0.0  # dummy values
            wells.append(well)
    for wk, well in enumerate(wells):
        well.id = wk
        if wk % 2 == 0:
            well.produce()
        else:
            well.inject(0)
    return wells


simulation.init(
    mesh=grid,
    wells=make_wells,
    cell_porosity=0.1,  # dummy value
    cell_permeability=1e-12,  # dummy value
    cell_thermal_conductivity=2,  # dummy value
)

print(
    f"Number of injectors on proc {mpi.proc_rank}: {simulation.number_of_own_injectors()} / {simulation.nb_injectors()}"
)
print(
    f"Number of producers on proc {mpi.proc_rank}: {simulation.number_of_own_producers()} / {simulation.nb_producers()}"
)


def message(well_type, data):
    for wk, well in enumerate(data):
        print(f"{well_type} well {wk} on proc {mpi.proc_rank} has id {well.id}")


message("production", simulation.producers_data())
message("injection", simulation.injectors_data())


dw.dump_producers(simulation)
dw.dump_injectors(simulation)

nx, ny, _ = grid.shape
nb_wells = (nx - 2) * (ny - 2)
assert all([0 <= well.id <= nb_wells for well in simulation.producers_data()])
