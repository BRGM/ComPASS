#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import ComPASS  # needed for cpp wrappers

def create_vertical_well(xy, well_radius = None, zmin=None, zmax=None):
    x, y, z = ComPASS.coordinates(ComPASS.global_vertices())
    x_well = x[np.argmin(np.abs(x - xy[0]))]
    y_well = y[np.argmin(np.abs(y - xy[1]))]
    selection = (x == x_well) & (y == y_well)
    if zmin is not None:
        selection&= (z >= zmin)
    if zmax is not None:
        selection&= (z <= zmax)
    well_nodes = np.nonzero(selection)[0]
    # CHECKME: What is the expected order for well nodes?
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    well = ComPASS.Well()
    if well_radius is None:
        well_radius = 0.1
    well.geometry.radius = well_radius
    segments = np.transpose(np.vstack([well_nodes[1:], well_nodes[:-1]]))
    well.geometry.add_segments(segments + 1) # Fortran indices start at 1
    return well
