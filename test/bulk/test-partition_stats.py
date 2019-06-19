#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

pres = 10. * MPa
Tres = 70. # degC

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

def fracture_factory(grid):
    def select_fracture():
        face_centers = ComPASS.compute_global_face_centers()
        dz = grid.extent[2] / grid.shape[2]
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = grid.origin[2] + 0.5 * grid.extent[2]
        return np.abs(face_centers[:, 2] - zfrac) < 0.25 * dz
    return select_fracture

ComPASS.init(
    mesh = grid,
    fracture_faces = fracture_factory(grid),
    cell_porosity = 0.5,                 # dummy value, no simulation
    cell_permeability = 1E-12,           # dummy value, no simulation
    cell_thermal_conductivity = 2,       # dummy value, no simulation
    fracture_porosity = 0.5,             # dummy value, no simulation
    fracture_permeability = 1E-12,       # dummy value, no simulation
    fracture_thermal_conductivity = 0.5, # dummy value, no simulation
)

@ComPASS.mpi.on_master_proc
def display_own_elements_number():
    for s in ['cells', 'faces', 'nodes', 'fractures']:
        parts = getattr(ComPASS.kernel, 'nb_%s_own' % s)()
        print('%s own:' % s, parts, 'total:', np.sum(parts))

display_own_elements_number()
