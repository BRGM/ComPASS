#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
from .._kernel import common_wrapper


def create_vertical_well(simulation, xy, well_radius=None, zmin=None, zmax=None):
    """
    :param simulation: simulation object, the method can also be accessed 
    through a fake method as `simulation.create_vertical_well`
    :param xy: the 2D coordinates (X, Y) of the well location
    :return: The data of the well which as `wid` id.
    """
    x, y, z = simulation.coordinates(simulation.global_vertices())
    x_well = x[np.argmin(np.abs(x - xy[0]))]
    y_well = y[np.argmin(np.abs(y - xy[1]))]
    selection = (x == x_well) & (y == y_well)
    if zmin is not None:
        selection &= z >= zmin
    if zmax is not None:
        selection &= z <= zmax
    well_nodes = np.nonzero(selection)[0]
    # CHECKME: What is the expected order for well nodes?
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    well = common_wrapper.Well()
    if well_radius is None:
        well_radius = 0.1
    well.geometry.radius = well_radius
    segments = np.transpose(np.vstack([well_nodes[1:], well_nodes[:-1]]))
    well.geometry.add_segments(segments + 1)  # Fortran indices start at 1
    return well


def get_well_data(simulation, wid):
    """
    :param simulation: simulation object, the method can also be accessed 
    through a fake method as `simulation.get_well_data`
    :param wid: well unique id
    :return: The data of the well which as `wid` id.
    """
    for wells in [simulation.injectors_data(), simulation.producers_data()]:
        for data in wells:
            if data.id == wid:
                return data


# WARNING: in parallel we must modify both own and ghost wells
def set_well_property(simulation, wid, **kwargs):
    """
    :param simulation: simulation object, the method can also be accessed 
    through a fake method (cf. example below)
    :param wid: well unique id
    Select data of the well which as `wid` id and set every property
    according to the kwargs dictionnart items 

    :Example:

    simulation.set_well_property(imposed_flowrate=100., injection_temperature=300.)
    """
    data = get_well_data(simulation, wid)
    if data is not None:
        for name, value in kwargs.items():
            setattr(data, name, value)
