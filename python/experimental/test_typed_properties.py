import numpy as np
from zones import ZoneManager
from properties import Property, TypedItems
from data_converters import SubArray, Vector, Tensor


class TypedProperty(TypedItems, Property):
    pass


def some_zones():
    m = ZoneManager(10)
    empty = m.from_indices([])
    full = ~empty
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))
    return empty, full, a, b


def test_scalar():
    empty, full, a, b = some_zones()
    p = TypedProperty(np.float64)
    p[...] = 0
    p[a] = 1
    p[b] = 2
    assert np.allclose(p[full], [2, 1, 2, 1, 2, 0, 2, 0, 2, 0])


def test_vector():
    empty, full, a, b = some_zones()
    p = TypedProperty(Vector(2))
    p[...] = [0, 0]
    p[a] = [1, 2]
    p[b] = [10, 20]
    assert np.allclose(p[full], [
        [10, 20], [1, 2], [10, 20], [1, 2], [10, 20],
        [0, 0], [10, 20], [0, 0], [10, 20], [0, 0]
    ])


def test_vector():
    empty, full, a, b = some_zones()
    p = TypedProperty(Tensor(2))
    p[...] = 0
    p[a] = 1
    p[b] = [[10, 20], [30, 40]]
    assert np.allclose(p[full], [
        [[10, 20], [30, 40]],
        [[1, 0], [0, 1]],
        [[10, 20], [30, 40]],
        [[1, 0], [0, 1]],
        [[10, 20], [30, 40]],
        [[0, 0], [0, 0]],
        [[10, 20], [30, 40]],
        [[0, 0], [0, 0]],
        [[10, 20], [30, 40]],
        [[0, 0], [0, 0]],
    ])
