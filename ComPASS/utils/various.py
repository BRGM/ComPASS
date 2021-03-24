from itertools import product
from pathlib import Path

import numpy as np

from .. import mpi


def enum_to_list(enum):
    return [
        name
        for name, _ in sorted(enum.__members__.items(), key=lambda item: int(item[1]))
    ]


def phases(simulation):
    return enum_to_list(simulation.Phase)


def components(simulation):
    return enum_to_list(simulation.Component)


def contexts(simulation):
    return enum_to_list(simulation.Context)


def states_locations(simulation):
    return [
        ("node", simulation.node_states()),
        ("cell", simulation.cell_states()),
        ("fracture", simulation.fracture_states()),
    ]


def mass_fluxes_locations(simulation):
    cell_fluxes, fracture_fluxes = simulation.mass_fluxes()
    return [
        ("cell", cell_fluxes),
        ("fracture", fracture_fluxes),
    ]


def tensor_coordinates(tensor, name, diagonal_only=False):
    tensor = np.asarray(tensor)
    assert tensor.ndim == 2 or tensor.ndim == 3, "wrong tensor array dimension"
    assert tensor.shape[-1] == tensor.shape[-2]
    dim = tensor.shape[-1]
    assert dim <= 3, "dimension should be 3 at max"
    if diagonal_only:
        return {
            f"{name}{si}{si}": tensor[..., i, i] for i, si in enumerate("xyz"[:dim])
        }
    # intended to be written to vtu as an array with numberofcomponents = 9
    # can be visualized in paraview with the following procedure
    # - create a Cell Centers filter as a child of the whole grid or a filter thereof
    # - create a Tensor Glyph filter as child of the Cell Centers
    # - in the Tensor Glyph filter's properties, select the property among the Tensors chooser
    # - adjust radius
    if tensor.ndim == 3:
        tensor = tensor.reshape([tensor.shape[0], tensor.shape[1] * tensor.shape[2]])
    return {name: tensor}


def _reload(simulation, snapshot, old_style):
    sep = "_" if old_style else " "
    states_locations = simulation.states_locations()
    phases = simulation.phases()
    if len(phases) == 1:
        phases = ["fluid"]
    components = simulation.components()
    for location, states in states_locations:
        if not old_style:
            states.context[:] = snapshot[f"{location} context"]
        else:
            states.context[:] = -1
        states.p[:] = snapshot[f"{location}{sep}pressure"]
        states.T[:] = snapshot[f"{location}{sep}temperature"]
        if len(phases) > 1:
            for phk, phase in enumerate(phases):
                states.S[:, phk] = snapshot[f"{location}{sep}{phase} saturation"]
        else:
            states.S.fill(1)
        if len(components) > 1:
            for ci, comp in enumerate(components):
                for phk, phase in enumerate(phases):
                    if old_style:
                        name = f"{location}_{comp} in {phase}"
                        if location == "fracture":  # this was a bug...
                            name = f"fracture_comp_{comp} in {phase}"
                    else:
                        name = f"{location} {comp} fraction in {phase}"
                    states.C[:, phk, ci] = snapshot[name]
        else:
            states.C.fill(1)


def reload_snapshot(
    simulation,
    path,
    iteration=None,
    verbose=True,
    old_style=False,
    reset_dirichlet=True,
):
    """
    This will reload a previous simulation state from snapshot outputs.
    The method can also be used as a *fake* simulation method: `simulation.reload_snapshot(path, iteration...)`.

    .. warning::
        The mesh and its partition must be exactly the same.
        Do not forget to reset Dirichlet conditions if necessary.

    .. warning::
        Using old style output physical context is not reloaded and will be set to -1.
        This is made on purpose to block flash, set context to the appropriate value...

    :param simulation: the simulation *object*
    :param path: the path to the ComPASS output directory that will be used to reload simulation state
    :param iteration: the ouput to reload (must be present in `path/snapshots` file)
                      if None (the default) the latest output from snapshots will be reloaded.
    :param verbose: if True will display a few information on master
                    proc about the reloaded snapshot (defaults to True).
    :param old_style: use old style output (defaults to False)
    :param reset_dirichlet: reset dirichlet node values (defaults to True)
    :return: the physical time in seconds of the reloaded snapshot
    """
    snapdir = Path(path)
    snapshot_info = None
    if mpi.is_on_master_proc:
        snapfile = snapdir / "snapshots"
        assert snapdir.is_dir(), f"Could not find snapshot directory: {str(snapdir)}"
        assert snapfile.is_file(), f"Could not find snapshot file: {str(snapfile)}"
        iterations = np.loadtxt(snapfile, usecols=(0,), dtype="i")
        times = np.loadtxt(snapfile, usecols=(1,), dtype="d")
        if iteration is None:
            assert len(iterations) > 0, "Snapshots file is empty..."
            index = len(iterations) - 1
        else:
            index = np.nonzero(iterations == iteration)[0]
            assert (
                len(index) != 0
            ), f"Could not find iteration {iteration} in snapshots (cf. {str(snapfile)})."
            assert (
                len(index) < 2
            ), f"Found several times iteration {iteration} in snapshots (cf. {str(snapfile)})."
            index = index[0]
        snapshot_info = (index, times[index])
    index, t = mpi.communicator().bcast(snapshot_info, root=mpi.master_proc_rank)
    snapshot = np.load(
        snapdir / "states" / f"state_{index:05d}_proc_{mpi.proc_rank:06d}.npz"
    )
    _reload(simulation, snapshot, old_style)
    if verbose:
        mpi.master_print(
            f"Reloaded snapshot {iteration} from {str(snapdir)} directory corresponding to time {t}"
        )
    if old_style:
        mpi.master_print(
            f"WARNING: Using old style output physical context is not reloaded and will be set to -1"
            f"WARNING: This is made on purpose to block flash"
            f"WARNING: Set context to the appropriate value..."
        )
    if reset_dirichlet:
        simulation.reset_dirichlet_nodes_states()
        if verbose:
            mpi.master_print(
                f"Dirichlet node states were updated according to reloaded snapshot."
            )
    return t


if __name__ == "__main__":
    T = np.random.random((2, 3, 3))
    print(tensor_coordinates(T, "T"))
