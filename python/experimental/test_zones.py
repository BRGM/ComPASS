import pytest

from numpy import allclose
from zones import ZoneManager


def test_from_indices():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))

    assert len(a) == 5
    assert len(b) == 5

    assert allclose(a.content, range(5))
    assert allclose(b.content, range(0, 10, 2))


def test_from_mask():
    m = ZoneManager(10)
    a = m.from_mask([1] * 5 + [0] * 5)
    b = m.from_mask([1, 0] * 5)

    assert allclose(a.content, range(5))
    assert allclose(b.content, range(0, 10, 2))


def test_set_operations():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))

    assert allclose((~a).content, range(5, 10))
    assert allclose((~b).content, range(1, 10, 2))
    assert allclose((a & b).content, [0, 2, 4])
    assert allclose((b & a).content, [0, 2, 4])
    assert allclose((a | b).content, [0, 1, 2, 3, 4, 6, 8])
    assert allclose((b | a).content, [0, 1, 2, 3, 4, 6, 8])
    assert allclose((a ^ b).content, [1, 3, 6, 8])
    assert allclose((b ^ a).content, [1, 3, 6, 8])
    assert allclose((a - b).content, [1, 3])
    assert allclose((b - a).content, [6, 8])

    assert not ~(a | ~a)
    assert not a & ~a


def test_set_comparison():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))

    assert not a.isdisjoint(b)
    assert not b.isdisjoint(a)

    assert a.isdisjoint(b - a)
    assert b.isdisjoint(a - b)

    assert bool(a)
    assert bool(b)
    assert not bool(a - a)

    assert a == a
    assert b == b

    assert a != b
    assert b != a

    assert not (a <= b)
    assert not (b <= a)

    assert a < (a | b) > b

    assert a != 235434
    assert not (a == 235434)


def test_multiple_managers():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))
    mm = ZoneManager(10)
    aa = mm.from_indices(range(5))
    bb = mm.from_indices(range(0, 10, 2))

    assert a != aa
    assert not a == aa

    with pytest.raises(ValueError):
        a & bb


def test_reinit_partiton():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))

    assert len(a) == 5
    assert len(b) == 5

    assert allclose(a.content, range(5))
    assert allclose(b.content, range(0, 10, 2))

    new_partition = [m.size - 1 - p for p in m.partition]
    new_zones = m._zones
    m.reinit(new_partition, new_zones)

    assert allclose(a.content, [5, 6, 7, 8, 9])
    assert allclose(b.content, [1, 3, 5, 7, 9])


def test_reinit_partiton_and_zones():
    m = ZoneManager(10)
    a = m.from_indices(range(5))
    b = m.from_indices(range(0, 10, 2))

    assert len(a) == 5
    assert len(b) == 5

    assert allclose(a.content, range(5))
    assert allclose(b.content, range(0, 10, 2))

    # très mal construit mais on garde que < 5
    new_partition = m.partition[:-2]
    new_zones = {
        id: frozenset(i for i in parts if i < 2) for id, parts in m._zones.items()
    }
    m.reinit(new_partition, new_zones)

    assert allclose(a.content, range(5))
    assert allclose(b.content, [range(0, 5, 2)])
