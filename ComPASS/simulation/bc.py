import numpy as np

from .. import mpi
from .._kernel import get_kernel


def _process_dirichlet_flags(dirichlet_flags, location1, location2):
    assert location1 is None or location2 is None
    # Node information is reset first ('i' flag)
    dirichlet_flags[:] = ord("i")
    if location1 is not None:
        dirichlet_flags[location1] = ord("d")
    elif location2 is not None:
        dirichlet_flags[location2] = ord("d")


def _process_dirichlet_nodes(
    info,
    dirichlet_nodes,
    pressure_dirichlet_nodes,
    temperature_dirichlet_nodes,
    energy_transfer=True,
):
    _process_dirichlet_flags(info.pressure, dirichlet_nodes, pressure_dirichlet_nodes)
    if energy_transfer:
        _process_dirichlet_flags(
            info.temperature, dirichlet_nodes, temperature_dirichlet_nodes
        )
    else:
        if temperature_dirichlet_nodes is not None:
            messages.error(
                "You are setting temperature dirichlet nodes without energy transfer"
            )


def set_global_dirichlet_nodes(
    simulation, both=None, pressure=None, temperature=None,
):
    assert mpi.is_on_master_proc
    kernel = get_kernel()
    info = np.rec.array(simulation.global_node_info(), copy=False)
    _process_dirichlet_nodes(
        info, both, pressure, temperature, simulation.has_energy_transfer_enabled(),
    )
    kernel.global_mesh_count_dirichlet_nodes()


def clear_dirichlet_nodes(simulation, update_scheme=True):
    def clear(a):
        a[:] = ord("i")

    info = np.rec.array(simulation.node_info(), copy=False)
    clear(info.pressure)
    if simulation.has_energy_transfer_enabled():
        clear(info.temperature)
    if update_scheme:
        simulation.scheme.compute_volumes()


def set_dirichlet_nodes(
    simulation, both=None, pressure=None, temperature=None,
):
    info = np.rec.array(simulation.node_info(), copy=False)
    _process_dirichlet_nodes(
        info, both, pressure, temperature, simulation.has_energy_transfer_enabled(),
    )


def reset_dirichlet_nodes_states(simulation):
    """
    Copy all dirichlet node states from the current node states.

    :param simulation: the simulation object
    """
    simulation.dirichlet_node_states().copy(simulation.node_states())


def reset_dirichlet_nodes(
    simulation,
    both_selection=lambda pts: None,
    pressure_selection=lambda pts: None,
    temperature_selection=lambda pts: None,
):
    """
    Select new Dirichlet nodes using mask functions
    (must be True where the node is a Dirichlet node)
    or a list of nodes that will become Dirichlet nodes.
    All dirichlet node states will be copied from the
    current node states.

    :param simulation: the simulation object
    :param both_selection: a function that will be called on vertices coordinates
                            to select both pressure and temperature dirichlet nodes
                            (default to no selection)
    :param pressure_selection: a function that will be called on vertices coordinates
                            to select both pressure dirichlet nodes (default to no selection)
    :param temperature_selection: a function that will be called on vertices coordinates
                            to select temperature dirichlet nodes (default to no selection)
    """
    clear_dirichlet_nodes(simulation, update_scheme=False)
    reset_dirichlet_nodes_states(simulation)

    def apply_if_callable(f):
        if callable(f):
            return f(simulation.vertices())
        return f

    set_dirichlet_nodes(
        simulation,
        both=apply_if_callable(both_selection),
        pressure=apply_if_callable(pressure_selection),
        temperature=apply_if_callable(temperature_selection),
    )
    simulation.scheme.compute_volumes()


def set_freeflow_faces(simulation, faces):
    """
    Select new freeflow faces using mask functions
    (must be True where the face is a *freeflow* face)
    or a sequence of faces id that will become freeflow faces.

    :param simulation: the simulation object
    :param faces: a sequence with boolean values or faces id
    """
    assert (
        simulation.mesh_is_local
    ), "Freeflow faces must be set once the mesh is distributed."
    kernel = get_kernel()
    faces = np.asarray(faces)
    assert faces.ndim == 1
    nf = simulation.number_of_faces()
    if faces.dtype == np.bool:
        assert (
            faces.size == nf
        ), "There should has many mask values has local mesh faces."
        faces = np.nonzero(faces)[0]
    assert np.all(faces >= 0) and np.all(faces < nf)
    kernel.set_freeflow_faces(faces + 1)  # C -> Fortran indexing


def reset_freeflow_faces(simulation, faces):
    """
    Delete all freeflow faces and select new ones using mask functions
    (must be True where the face is a *freeflow* face)
    or a sequence of faces id that will become freeflow faces.

    :param simulation: the simulation object
    :param faces: a sequence with boolean values or faces id
    """
    kernel.clear_freeflow_faces()
    set_freeflow_faces(simulation, faces)
