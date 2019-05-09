import pytest

import numpy as np

from zones import ZoneManager
from properties import Property, Initializer, Descriptor, Setter, Getter


class Base_Setter(Initializer, Setter):
    pass


class Setter_Getter(Initializer, Setter, Getter):
    pass


def some_zones():
    mng = ZoneManager(10)
    a = mng.from_indices(range(0, 5, 1))
    b = mng.from_indices(range(0, 10, 2))
    c = mng.from_indices(range(0, 10, 5))
    return a, b, c


@pytest.mark.parametrize("Property", [Initializer, Property])
def test_base_init(Property):
    p = Property()
    assert p._default is None
    assert p._data == {}


@pytest.mark.parametrize("Property", [Base_Setter, Property])
def test_set_default(Property):
    p = Property()
    p[...] = 99
    assert p._default == 99
    assert p._data == {}


@pytest.mark.parametrize("Property", [Base_Setter, Property])
def test_set_one_value(Property):
    a, b, c = some_zones()
    val = np.arange(a.manager.size) + 100
    p = Property()
    p[a] = val[a.content]
    assert p._default is None
    assert p._data.keys() == {a}
    assert p._data[a][0] == a
    assert all(p._data[a][1] == val[a.content])


@pytest.mark.parametrize("Property", [Base_Setter, Property])
def test_set_multiple_values(Property):
    a, b, c = some_zones()
    val = np.arange(a.manager.size) + 100
    p = Property()
    p[a] = val[a.content]
    p[b] = val[b.content]
    assert p._default is None
    assert p._data.keys() == {a - b, b}
    assert p._data[a - b][0] == a
    assert all(p._data[a - b][1] == val[a.content])
    assert p._data[b][0] == b
    assert all(p._data[b][1] == val[b.content])


@pytest.mark.parametrize("Property", [Base_Setter, Property])
def test_set_multiple_zone_managers(Property):
    a, b, c = some_zones()
    aa, bb, cc = some_zones()
    val = np.arange(a.manager.size) + 100
    p = Property()
    p[a] = val[a.content]
    p[b] = val[b.content]
    p[aa] = val[a.content] + 10
    p[bb] = val[b.content] + 10
    assert p._default is None
    assert p._data.keys() == {a - b, b, aa - bb, bb}
    assert p._data[a - b][0] == a
    assert all(p._data[a - b][1] == val[a.content])
    assert p._data[b][0] == b
    assert all(p._data[b][1] == val[b.content])
    assert p._data[aa - bb][0] == aa
    assert all(p._data[aa - bb][1] == val[aa.content] + 10)
    assert p._data[bb][0] == bb
    assert all(p._data[bb][1] == val[bb.content] + 10)


@pytest.mark.parametrize("Property", [Setter_Getter, Property])
def test_get_array_default(Property):
    a, b, c = some_zones()
    p = Property()
    p[...] = 99
    arr = p.get_array(a)
    assert len(a) == len(arr) and all(arr == 99)
    arr = p.get_array(b)
    assert len(b) == len(arr) and all(arr == 99)
    arr = p.get_array(c)
    assert len(c) == len(arr) and all(arr == 99)


@pytest.mark.parametrize("Property", [Setter_Getter, Property])
def test_get_array(Property):
    a, b, c = some_zones()
    val = np.arange(a.manager.size) + 100
    p = Property()
    p[a] = val[a.content]
    arr = p.get_array(a & b)
    assert all(arr == val[(a & b).content])


@pytest.mark.parametrize("Property", [Setter_Getter, Property])
def test_get_array_multiple_manager(Property):
    a, b, c = some_zones()
    aa, bb, cc = some_zones()
    val = np.arange(a.manager.size) + 100
    vval = np.arange(aa.manager.size) + 200
    p = Property()
    p[a] = val[a.content]
    p[aa] = vval[aa.content]
    arr = p.get_array(a & b)
    assert all(arr == val[(a & b).content])
    arr = p.get_array(aa & bb)
    assert all(arr == vval[(aa & bb).content])
