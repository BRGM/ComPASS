from itertools import chain
from collections import defaultdict
import weakref
import uuid
import numpy as np
from mpi4py import MPI


class Zone:
    """Base class for zone objects"""

    __slots__ = "__id__", "__weakref__"

    def __new__(cls, id):
        if cls is not Zone:
            return super(Zone, cls).__new__(cls)
        raise TypeError("can't instantiate Zone class")

    def __init__(self, id):
        assert id in self.manager.ids
        self.__id__ = id

    def __repr__(self):
        return f"<Zone({self.id!r}) of {self.manager!r}>"

    def __hash__(self):
        return hash(self.id)

    @property
    def id(self):
        return self.__id__

    @property
    def manager(self):
        raise NotImplementedError()


class BaseManager:
    """"""

    def __init__(self, size):
        # indices and partition
        self._size = size
        self._partition = (np.arange(size),)
        # zones : mapping id -> [partition_index, ...]
        self._zones = {}
        # zone class
        parents, namespace = [], {}
        self._prepare_cls(parents, namespace)
        self._zone_cls = type("<mng>.Zone", tuple(parents), namespace)

    def reinit(self, partition, zones):
        partition = tuple(partition)
        zones = dict(zones)
        size = sum(len(p) for p in partition)
        assert size - 1 == max(np.max(p) for p in partition)
        assert 0 == min(np.min(p) for p in partition)
        assert zones.keys() == self._zones.keys()
        self._size = size
        self._partition = partition
        self._zones = zones

    def _prepare_cls(self, parents, namespace):
        parents.insert(0, Zone)
        namespace.setdefault("__slots__", ())
        namespace["manager"] = property(lambda s: self)
        namespace["content"] = property(lambda s: self.get_content(s.id))

    @property
    def size(self):
        return self._size

    @property
    def Zone(self):
        return self._zone_cls

    @property
    def ids(self):
        return self._zones.keys()

    @property
    def partition(self):
        return self._partition

    def from_id(self, id):
        assert id in self._zones
        return self._zone_cls(id)

    def from_parts(self, parts):
        # TODO  validate input
        parts = frozenset(parts)
        old_id = self._id_from_parts(parts)
        if old_id is not None:
            return self._zone_cls(old_id)
        new_id = self._new_id()
        self._set_zone(new_id, parts)
        return self._zone_cls(new_id)

    def _new_id(self):
        return str(uuid.uuid4())

    def _id_from_parts(self, parts):
        for id, p in self._zones.items():
            if p == parts:
                return id
        return None

    def get_content(self, zone):
        chunks = [self.partition[i] for i in self._zones[zone]]
        if not chunks:
            return np.array([], dtype=int)
        res = np.concatenate(chunks)
        res.sort()
        return res

    def get_parts(self, zone_id):
        return self._zones[zone_id]

    def _update_partition(self, new_partition):
        # TODO  ? handle callbacks ?
        print("_update_partition")
        self._partition = tuple(new_partition)

    def _set_zone(self, id, parts):
        # TODO  ? handle callbacks ?
        print("_set_zone")
        self._zones[id] = frozenset(parts)


class Factory(BaseManager):
    """ """

    def from_mask(self, mask):
        assert len(mask) == self.size
        [indices] = np.where(mask)
        return self.from_indices(indices)

    def from_indices(self, indices):
        # TODO  check indices if needed
        indices = np.array(indices, dtype=int, copy=False)
        # build new partition, mapping with old partition and zone ids
        new_partition = list()
        parts_map = list()
        ids = list()
        count = 0
        for part in self.partition:
            local_map = list()
            intersect = np.intersect1d(part, indices)
            difference = np.setdiff1d(part, indices)
            if intersect.size != 0:
                intersect.setflags(write=False)
                new_partition.append(intersect)
                local_map.append(count)
                ids.append(count)
                count += 1
            if difference.size != 0:
                difference.setflags(write=False)
                new_partition.append(difference)
                local_map.append(count)
                count += 1
            assert local_map, "I don't know how this is possible"
            parts_map.append(local_map)
        # store new partition
        self._update_partition(new_partition)
        # update all zone parts, must be done before creating new zone
        old_parts = list(self._zones.items())
        for id, parts in old_parts:
            self._set_zone(id, chain.from_iterable(parts_map[i] for i in parts))
        # build zone
        return self.from_parts(ids)


class Calculator(BaseManager):
    """ """

    def _prepare_cls(self, parents, namespace):
        super()._prepare_cls(parents, namespace)

        EMPTY = object()

        def wrap(meth, default=EMPTY):
            # TODO  improve signature, doc, ...
            #       ? property for 1 arg methods ?
            if default is EMPTY:

                def wrapper(*zones):
                    if len({z.manager for z in zones}) != 1:
                        raise ValueError("zones must have the same manager")
                    args = (z.id for z in zones)
                    return meth(self, *args)

            else:

                def wrapper(*zones):
                    try:
                        args = [z.id for z in zones]
                    except AttributeError:
                        return default
                    if len({z.manager for z in zones}) != 1:
                        return default
                    return meth(self, *args)

            return wrapper

        for name, meth in vars(Calculator).items():
            if not name.startswith("zone_"):
                continue
            name = name[len("zone_") :]
            if name.endswith("_op"):
                name = f"__{name[:-len('_op')]}__"
            namespace[name] = wrap(meth)
        namespace["__hash__"] = lambda self: Zone.__hash__(self)
        namespace["__eq__"] = wrap(Calculator.zone_eq_op, False)
        namespace["__ne__"] = wrap(Calculator.zone_ne_op, True)

    def _bool_all(self, value):
        return bool(value)

    _bool_any = _bool_all

    def zone_isdisjoint(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1.isdisjoint(p2))

    def zone_len_op(self, zone):
        parts = self._zones[zone]
        return sum(len(self.partition[i]) for i in parts)

    def zone_bool_op(self, zone):
        parts = self._zones[zone]
        return self._bool_any(parts)

    def zone_or_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self.from_parts(p1 | p2)

    def zone_and_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self.from_parts(p1 & p2)

    def zone_sub_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self.from_parts(p1 - p2)

    def zone_xor_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self.from_parts(p1 ^ p2)

    def zone_invert_op(self, zone):
        parts = self._zones[zone]
        whole = frozenset(range(len(self.partition)))
        return self.from_parts(whole - parts)

    def zone_eq_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1 == p2)

    def zone_ne_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_any(p1 != p2)

    def zone_ge_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1 >= p2)

    def zone_le_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1 <= p2)

    def zone_gt_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1 > p2)

    def zone_lt_op(self, zone_1, zone_2):
        p1 = self._zones[zone_1]
        p2 = self._zones[zone_2]
        return self._bool_all(p1 < p2)


class Parallel:
    "add parallel (MPI) support"

    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.comm = MPI.COMM_WORLD

    def _new_id(self):
        "bcast"
        comm = self.comm
        id = super()._new_id() if comm.rank == 0 else None
        return comm.bcast(id, root=0)

    def _id_from_parts(self, parts):
        "gather, bcast"
        comm = self.comm
        ids = {id for id, p in self._zones.items() if p == parts}
        ids = comm.gather(ids, root=0)
        id = None
        if comm.rank == 0:
            ids = set.intersection(*ids)
            if ids:
                [id] = ids
        return comm.bcast(id, root=0)

    def _bool_all(self, value):
        "allreduce (all)"
        comm = self.comm
        return comm.allreduce(value, MPI.LAND)

    def _bool_any(self, value):
        "allreduce (any)"
        comm = self.comm
        return comm.allreduce(value, MPI.LOR)


class ZoneManager(Factory, Calculator):
    pass
