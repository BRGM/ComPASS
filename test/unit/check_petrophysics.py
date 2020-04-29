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
from ComPASS.timeloops import standard_loop, TimeStepManager

nx, ny, nz = 2, 2, 2     # discretization

simulation = ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape = (nx, ny, nz))

def select_fractures():
    fc = simulation.compute_global_face_centers()
    zfc = fc[:, 2]
    return np.fabs(zfc - 1) < 0.1 # middle horizontal faces

simulation.init(
    mesh = grid,
    cell_permeability = 1.,              # Dummy value
    cell_porosity = 0.5,                 # Dummy value
    cell_thermal_conductivity = 1.,      # Dummy value
    fracture_permeability = 1.,          # Dummy value
    fracture_porosity = 0.5,             # Dummy value
    fracture_thermal_conductivity = 1.,  # Dummy value
    fracture_faces = select_fractures,
)

petrophysics = simulation.petrophysics()
cp = simulation.get_cell_permeability()
# Check that we have a view and not a copy
cp[:] = 3*cp
assert np.all(cp==petrophysics.cell_permeability)

# Just test accessors
petrophysics.cell_permeability
petrophysics.cell_porosity
petrophysics.cell_thermal_conductivity
petrophysics.fracture_permeability
petrophysics.fracture_porosity
petrophysics.fracture_thermal_conductivity
