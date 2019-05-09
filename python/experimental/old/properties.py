
import copy
import numpy as np
from index_zoning import ZoneManager, Zone, Zict, ZoneMap
from distribute import ZoneArrayItem


class PropertyBase():
    def __init__(self, doc=None):
        self._data = dict()
        self._everywhere = None
        self.__doc__ = doc

    @property
    def everywhere(self):
        return self._get_everywhere()

    @everywhere.setter
    def everywhere(self, value):
        self._set_everywhere(value)

    def _get_everywhere(self):
        return self._everywhere

    def _set_everywhere(self, value):
        self._everywhere = value
        for manager in self._data:
            self[manager['whole']] = value

    def _get_data(self, manager):
        try:
            return self._data[manager]
        except KeyError:
            pass
        assert isinstance(manager, ZoneManager)
        # data = Zict()
        data = ZoneMap(manager)
        if self.everywhere is not None:
            # data[manager['whole']] = manager['whole'], self.everywhere
            data.set(manager['whole'], (manager['whole'], self.everywhere))
        self._data[manager] = data
        return data

    def _get_zone(self, key):
        if isinstance(key, ZoneManager):
            return key['whole']
        assert isinstance(key, Zone)
        return key

    def __setitem__(self, key, value):
        if key is Ellipsis:
            self.everywhere = value
            return
        key = self._get_zone(key)
        data = self._get_data(key.manager)
        # data[key] = value
        data.set(key, value)

    def __getitem__(self, key):
        raise NotImplementedError

    def _new(self):
        return PropertyBase(doc=self.__doc__)

    def __get__(self, obj, owner):
        if obj is None:
            return self
        name = next(name for cls in type(obj).mro()
                         for name, val in vars(cls).items()
                         if val is self)
        if name in obj.__dict__:
            return obj.__dict__[name]
        attr = obj.__dict__[name] = self._new()
        return attr

    def __set__(self, obj, value):
        attr = self.__get__(obj, type(obj))
        attr.everywhere = value



class ArrayProperty(PropertyBase):
    def __init__(self, type, doc=None):
        """

        type must:
            - be callable and return conversion of an object or raise exception
            - be convertible to dtype: np.dtype(type)
        """
        self.type = type
        super().__init__(doc=doc)

    @property
    def dtype(self):
        return np.dtype(self.type)

    def is_array(self, value, convert=False):
        value = self.type(value) if convert else value
        return value.shape != self.dtype.shape

    def _set_everywhere(self, value):
        value = self.type(value)
        # if value.shape != self.dtype.shape:
        if self.is_array(value):
            raise ValueError('value shape must be %s, not %s'
                             % (self.dtype.shape, value.shape))
        super()._set_everywhere(value)

    def __setitem__(self, key, value):
        value = self.type(value)
        key = self._get_zone(key)
        arr_shape = (len(key),) + self.dtype.shape
        # if value.shape != self.dtype.shape and value.shape != arr_shape:
        if self.is_array(value) and value.shape != arr_shape:
            raise ValueError("value shape must be %s or %s, not %s"
                             % (self.dtype.shape, arr_shape, value.shape))
        item = ZoneArrayItem(key, value, self.is_array(value))
        super().__setitem__(key, item)

    def __getitem__(self, key):
        key = self._get_zone(key)
        data = self._get_data(key.manager)
        if key.difference(*data.keys()):
            raise KeyError('data not available for the entire zone')
        idx = key.content()
        res = np.empty(len(key), self.dtype)
        for zone, (orig_zone, value) in data.items():
            zone = zone & key
            if not zone:
                continue
            ii = np.searchsorted(idx, zone.content())
            # if value.shape != self.dtype.shape:
            if self.is_array(value):
                jj = np.searchsorted(orig_zone.content(), zone.content())
                res[ii] = value[jj]
            else:
                res[ii] = value
        return res

    def _new(self):
        return ArrayProperty(self.type, doc=self.__doc__)


#############################################################################
#############################################################################


class Vector():
    def __init__(self, size, item_dtype=np.float64):
        assert isinstance(size, int)
        self._size = size
        self._item_dtype = np.dtype(item_dtype)

    def __repr__(self):
        return '%s(%s, %s)' % (type(self).__name__, self.size, self.item_dtype)

    @property
    def item_dtype(self):
        return self._item_dtype

    @property
    def size(self):
        return self._size

    @property
    def dtype(self):
        return np.dtype((self.item_dtype, self.size))

    def __call__(self, value):
        err = ValueError("could not understand input as %s data" % self)
        value = np.asarray(value, dtype=self.item_dtype)
        dtype = self.dtype
        if value.ndim == 1:
            if value.shape != dtype.shape:
                raise err
            return value
        elif value.ndim == 2:
            if value.shape[-dtype.ndim:] == dtype.shape:
                return value
            return self(np.moveaxis(value, -1, 0))
        raise err


class Tensor(Vector):
    @property
    def dtype(self):
        return np.dtype((self.item_dtype, (self.size, self.size)))

    def __call__(self, value):
        err = ValueError("could not understand input as %s data" % self)
        value = np.asarray(value, dtype=self.item_dtype)
        dtype = self.dtype
        if value.ndim == 0:
            res = np.zeros(dtype.shape, dtype=self.item_dtype)
            ii = np.arange(self.size)
            res[ii, ii] = value
            return res
        elif value.ndim == 1:
            pass
            res = np.zeros(value.shape+dtype.shape, dtype=self.item_dtype)
            for i in range(self.size):
                res[:, i, i] = value
            return res
        elif value.ndim == 2:
            if value.shape != dtype.shape:
                raise err
            return value
        elif value.ndim == 3:
            if value.shape[-dtype.ndim:] == dtype.shape:
                return value
            return self(np.moveaxis(value, -1, 0))
        raise err
