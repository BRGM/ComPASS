#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from pathlib import Path
import numpy as np
import ComPASS
from .. import mpi
from MeshTools import vtkwriters as vtkw


def cell_description(cell_connectivity):
    # copy is important not to modify ComPASS connectivity array
    cellnodes = np.copy(cell_connectivity.contiguous_content())
    cellnodes -= 1  # Fortran indexing -> C indexing
    cellnodes_offsets = cell_connectivity.offsets()
    return cellnodes, cellnodes_offsets


def mesh_description(nb_own_cells=None):
    mesh_type = ""  # local mesh
    if not ComPASS.mesh_is_local:
        mesh_type = "global_"
    vertices = getattr(ComPASS, f"{mesh_type}vertices")()
    connectivity = getattr(ComPASS, f"get_{mesh_type}connectivity")()
    cellnodes, offsets = cell_description(connectivity.NodebyCell)
    celltypes = getattr(ComPASS, f"{mesh_type}celltypes")()
    if nb_own_cells is not None:
        offsets = offsets[: nb_own_cells + 1]
        cellnodes = cellnodes[: offsets[-1]]
        celltypes = celltypes[:nb_own_cells]
    return vertices, offsets[1:], cellnodes, celltypes


def write_vtu_mesh(filename, pointdata, celldata, ofmt):
    nb_own_cells = ComPASS.number_of_own_cells()
    vertices, offsets, cellnodes, celltypes = mesh_description(nb_own_cells)
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


def write_mesh(basename, pointdata={}, celldata={}, ofmt="binary", ascontiguousarray=True):
    """

    ascontiguousarray: (default True)
        if True, convert the inputs arrays from `pointdata` and `celldata` as contiguous array
        using np.ascontiguousarray()

    """
    basename = Path(basename)
    vtu_directory = create_vtu_directory(basename.parent)
    filename = vtu_directory / proc_filename(basename, mpi.proc_rank)
    if ascontiguousarray:
        pointdata = {k: np.ascontiguousarray(v) for k, v in pointdata.items()}
        celldata = {k: np.ascontiguousarray(v) for k, v in celldata.items()}
    vertices_type = write_vtu_mesh(
        str(filename), pointdata=pointdata, celldata=celldata, ofmt=ofmt
    )
    nbprocs = mpi.communicator().size
    if nbprocs > 1 and mpi.is_on_master_proc:
        vtkw.write_pvtu(
            vtkw.pvtu_doc(
                vertices_type,
                [
                    str(vtu_directory / proc_filename(basename, pk))
                    for pk in range(nbprocs)
                ],
                pointdata_types=collect_dtypes(pointdata),
                celldata_types=collect_dtypes(celldata),
            ),
            str(basename),
        )