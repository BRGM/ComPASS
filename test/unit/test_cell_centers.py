#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

nx, ny, nz = 1, 1, 1

grid = ComPASS.Grid(shape = (nx, ny, nz))

ComPASS.init_and_load_mesh(grid)

print(ComPASS.compute_global_cell_centers())
print(ComPASS.compute_global_face_centers())