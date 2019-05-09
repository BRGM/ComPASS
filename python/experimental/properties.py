import numpy as np
from zones import Zone


class Initializer:
    """ initialize internal state.
    """

    def __init__(self):
        self._data = {}
        self._default = None

    def reinit(self, default, data):
        self._data.clear()
        self._default = None
        self[...] = default
        for zone, value in data:
            self[zone] = value


class Descriptor:
    """ descriptor object.
    `obj.attr = value` delegates to `obj.attr[...] = value`
    """

    def _new(self):
        return type(self)()

    def __get__(self, obj, owner):
        if obj is None:
            return self
        name = next(
            name
            for cls in owner.mro()
            for name, val in vars(cls).items()
            if val is self
        )
        if name in obj.__dict__:
            return obj.__dict__[name]
        attr = obj.__dict__[name] = self._new()
        return attr

    def __set__(self, obj, value):
        attr = self.__get__(obj, type(obj))
        attr[...] = value


class Setter:
    """ set item with zone objects as keys.

    `...` (ellipsis) is also a valid key, setting value to ellipsis key
    clears internal state and sets the default value.

    """

    def __setitem__(self, key, value):
        if key is Ellipsis:
            self._default = value
            self._data.clear()
        else:
            self._add_zone_value(key, value)

    def _add_zone_value(self, zone, value):
        ""
        assert isinstance(zone, Zone)
        if not zone:
            return
        zones = [z for z in self._data if z.manager is zone.manager]
        for old_zone in zones:
            new_zone = old_zone - zone
            zval = self._data.pop(old_zone)
            if new_zone:
                self._data[new_zone] = zval
        self._data[zone] = (zone, value)


class MultiKeySetter(Setter):
    """ set item with multiple keys.

    Setting value on two zones:
    >>> p[a, b] = value

    Other usage:
    >>> z = (reservoir_nodes, reservoir_faces, reservoir_cells)
    >>> porosity[z] = 0.3
    >>> permeability[z] = 1e-6
    >>> z = (caprock_nodes, caprock_faces, caprock_cells)
    >>> porosity[z] = 0.1
    >>> permeability[z] = 1e-9

    """

    def __setitem__(self, key, value):
        if isinstance(key, tuple):
            for single_key in key:
                super().__setitem__(single_key, value)
        else:
            super().__setitem__(key, value)


class Getter:
    """ get values from specified zone.
    """

    # NOTE  On pourrait faire du __getitem__ si le dtype étant connu par
    #       ailleurs.

    def iter_data(self, manager):
        return ((k, v) for k, v in self._data.items() if k.manager is manager)


class TypedItems(Setter, Getter):
    """
    """

    def __init__(self, type):
        super().__init__()
        self.type = type

    def __setitem__(self, key, value):
        value = self.type(value)
        super().__setitem__(key, value)

    def __getitem__(self, key):
        dtype = np.dtype(self.type)
        res = np.empty(len(key), dtype)
        key_content = key.content
        parts = key.manager.get_parts
        p_remains = parts(key.id)
        for zone, (orig_zone, value) in self.iter_data(key.manager):
            p_zone = parts(zone.id) & parts(key.id)
            if not p_zone:
                continue
            p_remains -= p_zone
            zone_content = zone.content
            ii = np.searchsorted(key_content, zone_content)
            value_as_array = (
                isinstance(value, np.ndarray)
                and value.ndim > dtype.ndim
            )
            if value_as_array:
                jj = np.searchsorted(orig_zone.content, zone_content)
                res[ii] = value[jj]
            else:
                res[ii] = value
        if p_remains:
            if self._default is not None:
                remains = key.manager.from_parts(p_remains)
                ii = np.searchsorted(key_content, remains.content)
                res[ii] = self._default
            else:
                pass
                # TODO  que faire si il reste des trous ?
                #       en particulier si pas de default ?
                #       idée: init prop._default = np.nan
        return res


class Property(Initializer, Descriptor, Setter, Getter):
    pass
