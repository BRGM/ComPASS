#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS

nx, ny, nz = 5, 5, 1

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(nx, ny, nz))

simulation.init(
    mesh=grid,
    cell_porosity=0.5,  # dummy value
    cell_permeability=1,  # dummy value
    cell_thermal_conductivity=1,  # dummy value
)

print(simulation.NumNodebyProc().as_array())
print(simulation.NumFracbyProc().as_array())
print(simulation.NumWellInjbyProc().as_array())
print(simulation.NumWellProdbyProc().as_array())
