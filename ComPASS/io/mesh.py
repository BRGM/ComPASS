#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from pathlib import Path
import numpy as np
from MeshTools import vtkwriters as vtkw
from .. import mpi


def cell_description(cell_connectivity):
    # copy is important not to modify ComPASS connectivity array
    cellnodes = np.copy(cell_connectivity.contiguous_content())
    cellnodes -= 1  # Fortran indexing -> C indexing
    cellnodes_offsets = cell_connectivity.offsets()
    return cellnodes, cellnodes_offsets


def mesh_description(simulation, nb_own_cells=None):
    mesh_type = ""  # local mesh
    if not simulation.mesh_is_local:
        mesh_type = "global_"
    vertices = getattr(simulation, f"{mesh_type}vertices")()
    connectivity = getattr(simulation, f"get_{mesh_type}connectivity")()
    cellnodes, offsets = cell_description(connectivity.NodebyCell)
    celltypes = getattr(simulation, f"{mesh_type}celltypes")()
    if nb_own_cells is not None:
        offsets = offsets[: nb_own_cells + 1]
        cellnodes = cellnodes[: offsets[-1]]
        celltypes = celltypes[:nb_own_cells]
    return vertices, offsets[1:], cellnodes, celltypes


def write_vtu_mesh(simulation, filename, pointdata, celldata, ofmt):
    nb_own_cells = simulation.number_of_own_cells()
    vertices, offsets, cellnodes, celltypes = mesh_description(simulation, nb_own_cells)
    assert all([a.shape[0] == vertices.shape[0] for a in pointdata.values()])
    assert all([a.shape[0] >= nb_own_cells for a in celldata.values()])
    celldata = {name: a[:nb_own_cells] for name, a in celldata.items()}
    vtkw.write_vtu(
        vtkw.vtu_doc_from_COC(
            vertices,
            offsets,
            cellnodes,
            celltypes,
            pointdata=pointdata,
            celldata=celldata,
            ofmt=ofmt,
        ),
        filename,
    )
    return vertices.dtype


def proc_filename(basename, procid):
    if mpi.communicator().size == 1:
        return basename
    ndigits = int(np.log10(mpi.communicator().size)) + 1
    ptag = f"{procid}".rjust(ndigits, "0")
    return f"{basename}_{ptag}.vtu"


def collect_dtypes(data):
    return {
        name: a.dtype if a.ndim == 1 else (a.dtype, a.shape[1]) 
        for name, a in data.items()
    }


def create_vtu_directory(parent=Path(".")):
    parallel = mpi.communicator().size > 1
    if parallel:
        vtu_directory = parent / "vtu"
        vtu_directory.mkdir(parents=True, exist_ok=True)
    mpi.synchronize()
    return vtu_directory if parallel else parent


def write_mesh(simulation, basename, pointdata={}, celldata={}, ofmt="binary"):
    """
    Write mesh data from a simulation as paraview vtu (single processor) or pvtu (multi processor) file.
    If there are more than one core in the communicator (multi processor case) the vtu files indexed
    in the pvtu files will be written in a vtu directory created alongside the pvtu file.

    :param simulation: a simulation object
    :param basename: the output filename (paraview)
    :param pointdata: a dictionnary of point based properties
    :param celldata: a dictionnary of cell based properties
    :param oftm: output format can be 'ascii' or 'binary' (the default)
    """
    basename = Path(basename)
    vtu_directory = create_vtu_directory(basename.parent)
    filename = vtu_directory / proc_filename(basename.name, mpi.proc_rank)
    # The following function ensures that properties have the right format
    # no copy is made if the underlying array already matches the requirements
    def check_properties(propdict):
        # possibly convert some values into numpy arrays
        propdict = {k: np.ascontiguousarray(v) for k, v in propdict.items()}
        # vtkwriter does not handle arrays of boolean so we convert them to numpy.int8
        propdict = {
            k: np.asarray(v, dtype=np.int8 if v.dtype == np.bool else v.dtype)
            for k, v in propdict.items()
        }
        return propdict

    pointdata = check_properties(pointdata)
    celldata = check_properties(celldata)
    vertices_type = write_vtu_mesh(
        simulation, str(filename), pointdata=pointdata, celldata=celldata, ofmt=ofmt
    )
    nbprocs = mpi.communicator().size
    # in the multi-processor case, the master proc writes a pvtu header file
    if nbprocs > 1 and mpi.is_on_master_proc:
        vtkw.write_pvtu(
            vtkw.pvtu_doc(
                vertices_type,
                [
                    str(
                        (vtu_directory / proc_filename(basename.name, pk)).relative_to(
                            basename.parent
                        )
                    )
                    for pk in range(nbprocs)
                ],
                pointdata_types=collect_dtypes(pointdata),
                celldata_types=collect_dtypes(celldata),
            ),
            str(basename),
        )
