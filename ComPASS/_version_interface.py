import packaging.version


def to_version(obj):
    return packaging.version.parse(str(obj))


def comp_op(str_op):
    return lambda self, other: getattr(self._version, str_op)(to_version(other))


class VersionInterface:
    def __init__(self, version):
        self._version = to_version(version)

    def __repr__(self):
        return str(self._version)

    __eq__ = comp_op("__eq__")
    __lt__ = comp_op("__lt__")
    __le__ = comp_op("__le__")
    __gt__ = comp_op("__gt__")
    __ge__ = comp_op("__ge__")
