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

def interwell_distance(grid):
    return grid.extent[0]/3

def center(grid):
    return tuple(grid.origin[i] + 0.5 * grid.extent[i] for i in range(3))

def select_boundary_factory(grid):
    on_xmin = lambda pts: pts[:,0] == grid.origin[0]
    on_xmax = lambda pts: pts[:,0] == grid.origin[0] + grid.extent[0]
    on_ymin = lambda pts: pts[:,1] == grid.origin[1]
    on_ymax = lambda pts: pts[:,1] == grid.origin[1] + grid.extent[1]
    def select():
        vertices = ComPASS.global_vertices().view(np.double).reshape((-1, 3))
        return on_xmin(vertices) | on_xmax(vertices) | on_ymin(vertices) | on_ymax(vertices)
    return select

def init_states(p, T):
    def set_states(states):
        states.context[:] = 2
        states.p[:] = p
        states.T[:] = T
        states.S[:] = [0, 1]
        states.C[:] = 1.
    for states in [ComPASS.dirichlet_node_states(),
                  ComPASS.node_states(),
                  ComPASS.fracture_states(),
                  ComPASS.cell_states()]:
        set_states(states)

def make_well(xy, well_radius = None):
    vertices = np.rec.array(ComPASS.global_vertices())
    x, y, z = vertices.x, vertices.y, vertices.z
    x_well = x[np.argmin(np.abs(x - xy[0]))]
    y_well = y[np.argmin(np.abs(y - xy[1]))]
    well_nodes = np.nonzero((x == x_well) & (y == y_well))[0]
    # CHECKME: What is the expected order for well nodes?
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    well = ComPASS.Well()
    if well_radius is None:
        well_radius = 0.1
    well.geometry.radius = well_radius
    segments = np.transpose(np.vstack([well_nodes[1:], well_nodes[:-1]]))
    well.geometry.add_segments(segments + 1) # Fortran indices start at 1
    return well

def make_wells_factory(grid, interwell=None):
    def make_wells(interwell=interwell):
        if interwell is None:
            interwell = interwell_distance(grid)
        Ox, Oy = center(grid)[:2]
        producer = make_well((Ox - 0.5 * interwell, Oy))
        #producer.operate_on_flowrate = 300, 1E5
        producer.operate_on_pressure = 10E6, 100 * ton / hour
        producer.produce()
        injector = make_well((Ox + 0.5 * interwell, Oy))
        #injector.operate_on_flowrate = 300, 30E6
        injector.operate_on_pressure = 30E6, 100 * ton / hour
        injector.inject(degC2K(30))
        return [producer, injector]
    return make_wells
