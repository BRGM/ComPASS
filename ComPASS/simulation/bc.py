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


def clear_dirichlet_nodes(simulation):
    def clear(a):
        a[:] = ord("i")

    info = np.rec.array(simulation.node_info(), copy=False)
    clear(info.pressure)
    if simulation.has_energy_transfer_enabled():
        clear(info.temperature)


def set_dirichlet_nodes(
    simulation, both=None, pressure=None, temperature=None,
):
    info = np.rec.array(simulation.node_info(), copy=False)
    _process_dirichlet_nodes(
        info, both, pressure, temperature, simulation.has_energy_transfer_enabled(),
    )


def reset_dirichlet_nodes(
    simulation,
    both_selection=lambda pts: None,
    pressure_selection=lambda pts: None,
    temperature_selection=lambda pts: None,
):
    clear_dirichlet_nodes(simulation)
    vertices = simulation.vertices()
    set_dirichlet_nodes(
        simulation,
        both=both_selection(vertices),
        pressure=pressure_selection(vertices),
        temperature=temperature_selection(vertices),
    )
