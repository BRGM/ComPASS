import functools
import numpy as np


def flatten_dtype(dtype):
    " return the corresponding flat sub-array dtype "
    dtype = np.dtype(dtype)
    shape = ()
    while True:
        sdt = dtype.subdtype
        if sdt is None:
            return np.dtype((dtype, shape))
        dtype, sub_shape = sdt
        shape += sub_shape


class SubArray:
    """
    """

    def __init__(self, dtype):
        self._dtype = flatten_dtype(dtype)

    @property
    def dtype(self):
        return self._dtype

    def __call__(self, value):
        dtype = self.dtype
        value = np.asarray(value, dtype=dtype.base)
        if value.shape[-dtype.ndim :] != dtype.shape:
            raise ValueError(f"shape {value.shape} not compatible with {dtype.shape}")
        return value

    def __repr__(self):
        return f"{type(self).__name__}<SubArray of dtype {self.dtype}>"


class Vector(SubArray):
    def __init__(self, dim, dtype=None):
        super().__init__((dtype, dim))


class Tensor(SubArray):
    def __init__(self, dim, dtype=None):
        super().__init__((dtype, (dim, dim)))

    def __call__(self, value):
        dtype = self.dtype
        shape = np.shape(value)
        if shape[-dtype.ndim :] != dtype.shape:
            new_value = np.zeros(shape + dtype.shape, dtype.base)
            for i in range(dtype.shape[-1]):
                new_value[..., i, i] = value
            value = new_value
        return super().__call__(value)


@functools.lru_cache()
def vector(size, dtype=None):
    """

    >>> my_type = vector(2)
    >>> vector_arr = np.empty(5, my_type)
    >>> vector_arr = my_type(data)

    """

    dtype = flatten_dtype((dtype, size))

    def convert(value):
        value = np.asarray(value, dtype=dtype.base)
        if value.ndim == dtype.ndim:
            if value.shape != dtype.shape:
                raise ValueError(f"can't convert data in vector({size})")
            return value
        elif value.shape[-dtype.ndim :] == dtype.shape:
            return value
            # return convert(np.moveaxis(value, -1, 0))
        raise ValueError(f"can't convert data in vector({size})")

    convert.dtype = dtype
    return convert
