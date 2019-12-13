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
from ComPASS.utils.wells import create_vertical_well


def interwell_distance(grid):
    return grid.extent[0]/3

def center(grid):
    return tuple(grid.origin[i] + 0.5 * grid.extent[i] for i in range(3))

def select_boundary_factory(simulation, grid):
    on_xmin = lambda pts: pts[:,0] == grid.origin[0]
    on_xmax = lambda pts: pts[:,0] == grid.origin[0] + grid.extent[0]
    on_ymin = lambda pts: pts[:,1] == grid.origin[1]
    on_ymax = lambda pts: pts[:,1] == grid.origin[1] + grid.extent[1]
    def select():
        vertices = simulation.global_vertices().view(np.double).reshape((-1, 3))
        return on_xmin(vertices) | on_xmax(vertices) | on_ymin(vertices) | on_ymax(vertices)
    return select

def init_states(simulation, p, T):
    def set_states(states):
        states.context[:] = simulation.Context.liquid
        states.p[:] = p
        states.T[:] = T
        states.S[:, simulation.phase_index(simulation.Phase.gas)] = 0
        states.S[:, simulation.phase_index(simulation.Phase.liquid)] = 1
        states.C[:] = 1
    for states in [simulation.dirichlet_node_states(),
                  simulation.node_states(),
                  simulation.fracture_states(),
                  simulation.cell_states()]:
        set_states(states)

def make_well(simulation, xy, well_radius = None):
    return create_vertical_well(simulation, xy, well_radius)

def make_wells_factory(simulation, grid, interwell=None):
    def make_wells(interwell=interwell):
        if interwell is None:
            interwell = interwell_distance(grid)
        Ox, Oy = center(grid)[:2]
        producer = make_well(simulation, (Ox - 0.5 * interwell, Oy))
        #producer.operate_on_flowrate = 300, 1E5
        producer.operate_on_pressure = 10E6, 100 * ton / hour
        producer.produce()
        injector = make_well(simulation, (Ox + 0.5 * interwell, Oy))
        #injector.operate_on_flowrate = 300, 30E6
        injector.operate_on_pressure = 30E6, 100 * ton / hour
        injector.inject(degC2K(30))
        return [producer, injector]
    return make_wells
