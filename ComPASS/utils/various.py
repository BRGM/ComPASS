from itertools import product
import numpy as np


def enum_to_list(enum):
    return [
        name
        for name, _ in sorted(enum.__members__.items(), key=lambda item: int(item[1]))
    ]


def phases(simulation):
    return enum_to_list(simulation.Phase)


def components(simulation):
    return enum_to_list(simulation.Component)


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
    return {
        f"{name}{si}{sj}": tensor[..., i, j]
        for (i, si), (j, sj) in product(enumerate("xyz"[:dim]), enumerate("xyz"[:dim]))
    }


if __name__ == "__main__":
    T = np.random.random((2, 3, 3))
    print(tensor_coordinates(T, "T"))
