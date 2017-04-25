import numpy as np
import ComPASS
from ComPASS.utils.units import *

def make_well(xy):
    vertices = ComPASS.get_vertices()
    x, y, z = (vertices[:, col] for col in range(3))
    x_well = x[np.argmin(np.abs(x - xy[0]))]
    y_well = y[np.argmin(np.abs(y - xy[1]))]
    well_nodes = np.nonzero((x == x_well) & (y == y_well))[0]
    # CHECKME: What is the expected order for well nodes?
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    well = ComPASS.Well()
    well.geometry.radius = 0.1
    segments = np.transpose(np.vstack([well_nodes[:-1], well_nodes[1:]]))
    well.geometry.add_segments(segments + 1) # Fortran indices start at 1
    return well

def make_wells_factory(grid):
    def make_wells():
        interwell_distance = grid.extent[0]/3
        producer = make_well((-0.5 * interwell_distance, 0))
        #producer.operate_on_flowrate = 300, 1E5
        producer.operate_on_pressure = 10E6, 100 * ton / hour
        producer.produce()
        injector = make_well((0.5 * interwell_distance, 0))
        #injector.operate_on_flowrate = 300, 30E6
        injector.operate_on_pressure = 30E6, 100 * ton / hour
        injector.inject(degC2K(30))
        return [producer, injector]
    return make_wells