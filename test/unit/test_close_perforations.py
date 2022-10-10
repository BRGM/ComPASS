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
import ComPASS.io.mesh as io
from ComPASS.dump_wells import print_producers_stats

# A vertical well with 2x2 grid basis and nv horizontal layers over H thickness

nx = ny = nz = 10
L = 100
wid = 12  # well id

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(nx, ny, nz), extent=(L, L, L))


def make_well():
    well = simulation.create_vertical_well((0.5 * L, 0.5 * L))
    well.id = wid
    well.operate_on_flowrate = 1.0, -np.inf
    well.produce()
    return [well]


simulation.init(
    mesh=grid,
    wells=make_well,
    cell_porosity=0.5,
    cell_permeability=1.0,
    cell_thermal_conductivity=1.0,
)


simulation.close_perforations(wid, above=0.8 * L, below=0.2 * L)
print_producers_stats(simulation)
simulation.close_perforations(wid)
print_producers_stats(simulation)
