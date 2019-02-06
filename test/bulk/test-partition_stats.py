#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

ComPASS.load_eos('water2ph')

pres = 10. * MPa
Tres = 70. # degC

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = grid,
    cell_porosity = 0.5,               # dummy value, no simulation
    cell_permeability = 1E-12,         # dummy value, no simulation
    cell_thermal_conductivity = 2,     # dummy value, no simulation
)

@ComPASS.mpi.on_master_proc
def display_own_elements_number():
    print("cells own:", ComPASS.kernel.nb_cells_own())
    print("faces own:", ComPASS.kernel.nb_faces_own())
    print("nodes own:", ComPASS.kernel.nb_nodes_own())
    print("fractures own:", ComPASS.kernel.nb_fractures_own())

display_own_elements_number()
