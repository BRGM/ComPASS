"""

Case connait les property, les locations (zone mng) et d'autres trucs


"""


class Case:
    def distribute(self, *args, **kwds):
        data = self._prepare_distribute(*args)
        self._distribute_mesh(data, *args, **kwds)
        self._distribute_properties(data, *args, **kwds)
        self._distribute_zones(data, *args, **kwds)
