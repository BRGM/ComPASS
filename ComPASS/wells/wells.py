#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
from .._kernel import common_wrapper
from .. import mpi
from ..timeloops import Event
from ..mpi import master_print
from ..dump_wells import _wells_info
from ..exceptions import CompassException


def create_well_from_segments(simulation, segments, well_radius=None):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.create_well_from_segments`
    :param segments: a sequence of pair vertices id oritented from wellhead downwards
    :param well_radius: the well radius in meters (used to compute Peaceman well indices - defaults to 0.1 m)
    """
    well = common_wrapper.Well()
    well.geometry.radius = well_radius or 0.1
    segments = np.ascontiguousarray(segments)
    segments.shape = -1, 2
    well.geometry.add_segments(segments + 1)  # Fortran indices start at 1
    return well


def create_single_branch_well(simulation, nodes, well_radius=None):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.create_single_branch_well`
    :param segments: a sequence of vertices id describing the well from top to bottom
    :param well_radius: the well radius in meters (used to compute Peaceman well indices - defaults to 0.1 m)
    """
    nodes = np.asarray(nodes)
    assert nodes.ndim == 1
    segments = np.transpose(np.vstack([nodes[:-1], nodes[1:]]))
    return create_well_from_segments(simulation, segments, well_radius)


def create_vertical_well(simulation, xy, well_radius=None, zmin=None, zmax=None):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.create_vertical_well`

    :param xy: the 2D coordinates (X, Y) of the well location
    :param well_radius: the well radius in meters (used to compute Peaceman well indices - defaults to 0.1 m)
    :param zmin: only vertices above zmin will be considered (ignored if None, which is the default)
    :param zmax: only vertices above zmin will be considered (ignored if None, which is the default)

    :return: the created well
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
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    # reorder nodes from top to bottom
    return create_single_branch_well(simulation, well_nodes[::-1], well_radius)


def _get_well_data(simulation, wells_data, wid):
    for k, data in enumerate(wells_data):
        if data.id == wid:
            return k, data


def _all_wells_data(simulation, own_only):
    return [simulation.injectors_data(own_only), simulation.producers_data(own_only)]


def get_well_data(simulation, wid, own_only=False):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.get_well_data`

    :param wid: well unique id
    :param own_only: will look amongst own wells only
                     (well whose unknowns are managed by this proc)

    :return: The data of the well which as `wid` id.
    """
    for wells_data in _all_wells_data(simulation, own_only):
        data = _get_well_data(simulation, wells_data, wid)
        if data is not None:
            return data[1]


def get_well_perforations_state(simulation, wid, own_only=False):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.get_well_data`

    :param wid: well unique id
    :param own_only: will look amongst own wells only
                     (well whose unknowns are managed by this proc)

    :return: The perforations of the well which as `wid` id.
    """
    for getdata, getperf in [
        (simulation.injectors_data, simulation.injector_perforations),
        (simulation.producers_data, simulation.producer_perforations),
    ]:
        data = _get_well_data(simulation, getdata(own_only), wid)
        if data is not None:
            wk, _ = data
            return getperf(wk)


def get_wellhead(simulation, wid, own_only=False):
    """
    :param simulation: simulation object, the method can also be accessed
                       through a fake method as `simulation.get_well_data`

    :param wid: well unique id
    :param own_only: will look amongst own wells only
                     (well whose unknowns are managed by this proc)

    :return: The wellhead of the well which as `wid` id.
    """
    perfs = get_well_perforations_state(simulation, wid, own_only)
    if perfs is not None:
        return perfs.wellhead


# WARNING: in parallel we must modify both own and ghost wells
def set_well_property(simulation, wid, verbose=False, **kwargs):
    """
    Select data of the well which as `wid` id and set every property
    according to the `kwargs` dictionnary items

    :param simulation: simulation object, the method can also be accessed
                       through a fake method (cf. example below)

    :param wid: well unique id

    :Example:

    .. highlight:: python
    .. code-block:: python

        simulation.set_well_property(wid, imposed_flowrate=100., injection_temperature=300.)
    """
    data = get_well_data(simulation, wid)
    if data is not None:
        for name, value in kwargs.items():
            if verbose:
                print(
                    f"Setting {name} for well {wid} to {value} on proc {mpi.proc_rank}"
                )
            setattr(data, name, value)


# WARNING: in parallel we must modify both own and ghost wells
def close_well(simulation, wid):
    """
    Close the well which as `wid` id.

    :param simulation: simulation object, the method can also be accessed
                       through a fake method (cf. example below)

    :param wid: well unique id

    :Example:

    .. highlight:: python
    .. code-block:: python

        simulation.close_well(wid)
    """
    data = get_well_data(simulation, wid)
    if data is not None:
        # FIXME: to be generalized
        assert (
            data.operating_code in "cf"
        ), f"Only wells operating on flowrate can be closed (found operating code {data.operating_code})."
        data.close()


# WARNING: in parallel we must modify both own and ghost wells
def open_well(simulation, wid):
    """
    Open the well which as `wid` id and set its operating mode to *flowrate*.

    :param simulation: simulation object, the method can also be accessed
                       through a fake method (cf. example below)

    :param wid: well unique id

    :Example:

    .. highlight:: python
    .. code-block:: python

        simulation.open_well(wid)
    """
    data = get_well_data(simulation, wid)
    if data is not None:
        assert (
            data.operating_code == "c"
        ), f"Only closed wells can be opened and set operated on flowrate (found operating code {data.operating_code})."
        data.open()


def _make_close_well_action(simulation, wid, verbose):
    def action(_):
        if verbose:
            master_print(f"WELL OPERATION: closing well {wid}")
        close_well(simulation, wid)

    return action


def _make_open_well_action(simulation, wid, verbose):
    def action(_):
        if verbose:
            master_print(f"WELL OPERATION: opening well {wid}")
        open_well(simulation, wid)

    return action


# OPTIMIZE: we could by pass the search for the well at each loop...
#           and register actions only on concerned procs
# WARNING: in parallel we must modify both own and ghost wells
def well_production_history(simulation, wid, history, verbose=True):
    """
    Set well production history.
    """

    def make_set_flowrate_action(simulation, wid, qw, verbose):
        def action(_):
            if verbose:
                master_print(
                    f"WELL OPERATION: setting well {wid} production flowrate to {qw}"
                )
            data = get_well_data(simulation, wid)
            if data is not None:
                if data.is_closed:
                    data.open()
                assert data.operating_code == "f"
                assert qw > 0
                data.imposed_flowrate = qw

        return action

    close_this_well = _make_close_well_action(simulation, wid, verbose)
    open_this_well = _make_open_well_action(simulation, wid, verbose)
    events = []
    for t, qw in history:
        actions = []
        if qw == 0:
            actions.append(close_this_well)
        else:
            assert qw > 0, "Production flow rate should be positive or null."
            actions.append(make_set_flowrate_action(simulation, wid, qw, verbose))
        events.append(Event(t, actions))
    return events


# OPTIMIZE: we could by pass the search for the well at each loop...
#           and register actions only on concerned procs
# WARNING: in parallel we must modify both own and ghost wells
def well_injection_history(simulation, wid, history, verbose=True):
    """
    Set well injection history.
    """

    def make_set_flowrate_action(simulation, wid, qw, T, verbose):
        def action(_):
            if verbose:
                master_print(
                    f"WELL OPERATION: setting well {wid} injection flowrate to {qw} with temperature {T}"
                )
            data = get_well_data(simulation, wid)
            if data is not None:
                if data.is_closed:
                    data.open()
                assert data.operating_code == "f"
                assert qw < 0 and T > 273.15  # 0°C
                data.imposed_flowrate = qw
                data.injection_temperature = T

        return action

    close_this_well = _make_close_well_action(simulation, wid, verbose)
    open_this_well = _make_open_well_action(simulation, wid, verbose)
    events = []
    for t, qw, T in history:
        actions = []
        if qw == 0:
            actions.append(close_this_well)
        else:
            actions.append(
                make_set_flowrate_action(simulation, wid, -abs(qw), T, verbose)
            )
        events.append(Event(t, actions))
    return events


def close_perforations(simulation, wid, above=None, below=None):
    """
    Close some of the well perforations by setting their Peaceman well indices to 0.
    (This is a rough approximation).

    :param simulation: simulation object, the method can also be accessed
                       through a fake method (cf. example below)

    :param wid: well unique id
    :param above: optional, will close all perforations with z above `above`
    :param below: optional, will close all perforations with z below `below`

    If neither `above` nor `below` are specified all perforations are closed.

    :Example:

    .. highlight:: python
    .. code-block:: python

        simulation.close_perforations(wid, above=z_top_reservoir)
    """
    well = None
    for well_type in ["injection", "production"]:
        pack = _wells_info(simulation, well_type)
        if pack.nb > 0:
            for info, data in zip(pack.information, pack.data):
                if data.id == wid:
                    well = info
                    break
        if well is not None:
            break
    if well is None:
        return
    vertices = simulation.vertices()[well.vertices]
    z = vertices[:, 2]

    def close_perfs(mask):
        well.well_index_Darcy[mask] = 0
        well.well_index_Fourier[mask] = 0

    if above is not None:
        close_perfs(z >= float(above))
    if below is not None:
        close_perfs(z <= float(below))
    if above is None and below is None:
        close_perfs(np.ones(z.shape, dtype=np.bool))


def set_well_model(simulation, well_model):
    if simulation.well_model is not None:
        raise CompassException("You cannot change well model at runtime.")
    if not (well_model == "single_phase" or well_model == "two_phases"):
        raise CompassException("Well model must be single_phase or two_phases.")
    simulation.well_model = well_model
