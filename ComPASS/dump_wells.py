from collections import namedtuple
from pathlib import Path
import numpy as np
from ComPASS.utils.units import *
import ComPASS.mpi as mpi
import MeshTools.vtkwriters as vtkw


WellInfo = namedtuple("WellInfo", ["type", "nb_own", "nb", "information", "data"])


def _injectors_info(simulation):
    well_type = "injection"
    n = simulation.nb_injectors()
    if n == 0:
        return WellInfo(well_type, None, 0, None, None)
    return WellInfo(
        well_type,
        simulation.number_of_own_injectors(),
        n,
        simulation.injectors_information(),
        simulation.injectors_data(),
    )


def _producers_info(simulation):
    well_type = "production"
    n = simulation.nb_producers()
    if n == 0:
        return WellInfo(well_type, None, 0, None, None)
    return WellInfo(
        well_type,
        simulation.number_of_own_producers(),
        n,
        simulation.producers_information(),
        simulation.producers_data(),
    )


def _wells_info(simulation, well_type):
    assert well_type in ["injection", "production"]
    if well_type == "injection":
        return _injectors_info(simulation)
    return _producers_info(simulation)


def _check_own_wells(wells, info):
    assert info.nb == wells.nb_wells
    nb_own_wells = info.nb_own
    if nb_own_wells is None:
        nb_own_wells = wells.nb_wells
    return nb_own_wells


def _any_wells(info):
    if info.nb == 0:
        print()
        print(f"No {info.type} wells.")
        print()
        return False
    return True


def _print_well_stats(simulation, well_type):
    info = _wells_info(simulation, well_type)
    if not _any_wells(info):
        return
    wells = info.information
    nb_own_wells = _check_own_wells(wells, info)
    well_id = [data.id for data, _ in zip(info.data, range(nb_own_wells))]
    print()
    print(f"Statistics for {info.nb_own} {info.type} wells out of {info.nb}.")
    for wk, well in enumerate(wells):
        if wk >= nb_own_wells:
            break
        print(
            f"Well {wk} with id {well_id[wk]} has {well.nb_perforations} perforations."
        )
        print("vertices:", well.vertices)
        print("parent (vertex id):", well.parent_vertex)
        print("parent (rank):", well.parent_rank)
        print("pressure:", well.pressure)
        print("temperature:", K2degC(well.temperature))
        print("pressure_drop:", well.pressure_drop)
        print("density:", well.density)
        print("Darcy WI:", well.well_index_Darcy)
        print("Fourier WI:", well.well_index_Fourier)
    print()


def _well_vtu_filename(wellid, tag=None):
    if tag is None:
        tag = ""
    else:
        tag = f"_{tag}"
    return f"well_id{wellid:04d}{tag}.vtu"


def _dump_wells(simulation, well_type, outputdir=None, tag=None, verbose=False):
    if outputdir is None:
        outputdir = Path(".")
    else:
        outputdir = Path(outputdir)
    assert outputdir.exists() and outputdir.is_dir()
    info = _wells_info(simulation, well_type)
    if not _any_wells(info):
        return
    wells = info.information
    nb_own_wells = _check_own_wells(wells, info)
    if verbose:
        print(f"Dumping {info.nb_own} {info.type} wells out of {info.nb}.")
    wells_data = [data for data, _ in zip(info.data, range(nb_own_wells))]
    vertices = simulation.vertices()
    node_states = simulation.node_states()
    if len(set([data.id for data in wells_data])) != nb_own_wells:
        print("!!!")
        print("!!! WARNING: you must use unique well id to discriminate output files")
        print("!!!")
    for wk, well in enumerate(wells):
        if wk >= nb_own_wells:
            break
        # We want to keep only well vertices not all revervoir vertices
        well_vertices = well.vertices
        assert well_vertices.shape == np.unique(well_vertices).shape
        remap = np.zeros(vertices.shape[0], dtype=well_vertices.dtype)
        remap[well_vertices] = 1
        remap = np.cumsum(remap) - 1

        vtkw.write_vtu(
            vtkw.vtu_doc(
                vertices[well_vertices],
                np.hstack(
                    [
                        np.reshape(remap[well_vertices[:-1]], (-1, 1)),
                        np.reshape(remap[well.parent_vertex], (-1, 1)),
                    ]
                ),
                pointdata={
                    name: np.ascontiguousarray(a)
                    for name, a in [
                        ("reservoir vertices id", well_vertices),
                        ("well pressure", well.pressure),
                        ("well temperature", K2degC(well.temperature)),
                        ("well pressure drop", well.pressure_drop),
                        ("well density", well.density),
                        ("Darcy WI", well.well_index_Darcy),
                        ("Fourier WI", well.well_index_Fourier),
                        ("reservoir pressure", node_states.p[well_vertices]),
                        ("reservoir temperature", K2degC(node_states.T[well_vertices])),
                        ("well saturation pressure", simulation.Psat(well.temperature)),
                        (
                            "inflow",
                            np.where(
                                well.pressure < node_states.p[well_vertices],
                                well.well_index_Darcy
                                * (node_states.p[well_vertices] - well.pressure),
                                0,
                            ),
                        ),
                    ]
                },
            ),
            str(outputdir / _well_vtu_filename(wells_data[wk].id, tag)),
        )


def dump_injectors(simulation, outputdir=None, tag=None, verbose=False):
    _dump_wells(simulation, "injection", outputdir, tag, verbose)


def dump_producers(simulation, outputdir=None, tag=None, verbose=False):
    _dump_wells(simulation, "production", outputdir, tag, verbose)


def dump_all_wells(simulation, outputdir=None, tag=None, verbose=False):
    for well_type in ("injection", "production"):
        _dump_wells(simulation, well_type, outputdir, tag, verbose)


def print_injectors_stats(simulation):
    _print_well_stats(simulation, "injection")


def print_producers_stats(simulation):
    _print_well_stats(simulation, "production")


def collect_well_ids(simulation):
    well_ids = []
    for well_type in ("injection", "production"):
        info = _wells_info(simulation, well_type)
        if info.nb == 0:
            continue
        wells = info.information
        nb_own_wells = _check_own_wells(wells, info)
        well_ids += [data.id for data, _ in zip(info.data, range(nb_own_wells))]
    return well_ids
