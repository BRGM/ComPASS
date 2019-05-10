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

def test_face_centers():

    ComPASS.load_eos('water2ph')
    ComPASS.set_output_directory_and_logfile(__file__)

    Lx, Ly, Lz = 100., 100., 100.
    nx, ny, nz = 10, 10, 10

    grid = ComPASS.Grid(
        shape = (nx, ny, nz),
        extent = (Lx, Ly, Lz),
    )

    def fractures():
        centers = ComPASS.compute_global_face_centers()
        selection = centers[:, 0] == 0.5 * Lx
        return selection

    ComPASS.init(
        mesh = grid,
        fracture_faces = fractures,
        cell_porosity = 0.5,               # dummy value, no simulation
        cell_permeability = 1E-12,         # dummy value, no simulation
        cell_thermal_conductivity = 2,     # dummy value, no simulation
        fracture_porosity = 0.5,           # dummy value, no simulation
        fracture_permeability = 1E-10,     # dummy value, no simulation
        fracture_thermal_conductivity = 2, # dummy value, no simulation
    )

    #print('Old computation')
    #print(ComPASS.old_compute_face_centers())
    #print('New computation')
    #print(ComPASS.compute_face_centers())
    assert np.all((len(ComPASS.old_compute_face_centers())==0 and len(ComPASS.compute_face_centers())==0) or
                   ComPASS.old_compute_face_centers()==ComPASS.compute_face_centers())
    #print('Fractures old computation')
    #print(ComPASS.old_compute_fracture_centers())
    #print('Fractures new computation')
    #print(ComPASS.compute_fracture_centers())
    assert np.all((len(ComPASS.old_compute_fracture_centers())==0 and len(ComPASS.compute_fracture_centers())==0) or
                   ComPASS.old_compute_face_centers()==ComPASS.compute_face_centers())


if __name__ == '__main__':
    test_face_centers()