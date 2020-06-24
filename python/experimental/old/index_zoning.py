""" Tools for managing subsets of indices
"""

__all__ = ["ZoneManager", "Zone", "Zict"]


import numpy as np
import weakref

# TODO

#   Zone.__new__(block_ids) delegue tout au manager. [attention a la recursion]


class Zone(object):
    """
    Zone objects represent subsets of indices enhanced with set algebra
    operations.

    Set algebra operations are the same as provided by the built-in set
    class. One additional operation is provided, `zone.complement()`,
    which returns a zone representing the complementary set of the zone.

    A zone object is always built through its manager, ie. a
    ZoneManager instance. Directly instanciating the Zone class will
    raise an error.

    Internally, a zone object does not hold the list of indices it
    represents. Two methods give acces to this list:

        - zone.content()

          creates and returns a new array contening the indices.

        - zone.chunked_content()

          made to avoid to allocate and copy arrays. It returns a tuple
          of references to arrays of indices held by the manager.

    zone.content() creates and returns a new array contening the indices.
    zone.chunked_content() is made to avoid to allocate and copy arrays.
    It returns a tuple of references to arrays of indices held by the
    manager.

    See also: ZoneManager
    """

    def __new__(cls):
        if cls is Zone:
            raise TypeError(
                "Can't directly instantiate Zone class: "
                "see ZoneManager to produce Zone instances."
            )
        return super(Zone, cls).__new__(cls)

    def __init__(self):
        self._content_ref = None  # used by self.content()
        self._block_ids = None  # defined by ZoneManager._create_new_zone()

    def __bool__(self):
        return bool(self._block_ids)

    __nonzero__ = __bool__

    def __len__(self):
        "len(self) -> size of subset of indices represented by the zone"
        partition = self.manager.partition
        return sum(len(partition[i]) for i in self._block_ids)

    @property
    def manager(self):
        "Return the manager of instance class."
        return self.__class__._manager

    def content(self):
        """Return the sorted array contening all the identifiers of the zone.

        Remarks:
            This method create a new array. To avoid copying arrays, see the
            zone.chunked_content() method.
            Nevertheless, it will keep an internal wekref to the result. Thus
            the array will not be rebuilt if a strong reference already exists
            elsewhere.
        """
        if (self._content_ref is None) or (self._content_ref() is None):
            content = self._build_content()
            content.setflags(write=False)
            self._content_ref = weakref.ref(content)
        return self._content_ref()

    def _build_content(self):
        "build the sorted array contening all the identifiers of the zone"
        chunks = self.chunked_content()
        if chunks:
            res = np.concatenate(chunks)
            res.sort()
            return res
        else:  # concatenate can't manage void lists
            return np.array([], dtype=self.manager.dtype)

    def chunked_content(self):
        """Return a tuple of arrays contening the identifiers of the zone.
        Returned arrays are disjoints.

        Returned arrays are references to internal arrays of the manager.
        Thus, unless manager must change its internal partition (eg. if
        manager.build_from_indices() is called), the memory footprint is null.
        """
        return tuple(self.manager._partition[i] for i in self._block_ids)

    def _valid_others(self, *others):
        "raise an error if others have not the same type as self."
        if any(not isinstance(o, self.__class__) for o in others):
            msg = "Set algebra avaible only for zones with same manager."
            raise TypeError(msg)

    ## Methods for set algebra
    def complement(self):
        "Return the complement of the zone."
        ids = set(range(len(self.manager._partition))) - self._block_ids
        return self.manager._create_new_zone(ids)

    def __eq__(self, other):
        "self == other -> True if zones are the sames."
        self._valid_others(other)
        return self._block_ids == other._block_ids

    def __ne__(self, other):
        "self != other -> False if zones are the sames."
        self._valid_others(other)
        return self._block_ids != other._block_ids

    def union(self, *others):
        "Return the union of two or more zones as a new zone."
        self._valid_others(*others)
        ids = self._block_ids.union(*(o._block_ids for o in others))
        return self.manager._create_new_zone(ids)

    def __or__(self, other):
        "self | other <=> self.union(other)"
        return self.union(other)

    def intersection(self, *others):
        "Return the intersection of two or more zones as a new zone."
        self._valid_others(*others)
        ids = self._block_ids.intersection(*(o._block_ids for o in others))
        return self.manager._create_new_zone(ids)

    def __and__(self, other):
        "self & other <=> self.intersection(other)"
        return self.intersection(other)

    def difference(self, *others):
        "Return the difference of two or more zones as a new zone."
        self._valid_others(*others)
        ids = self._block_ids.difference(*(o._block_ids for o in others))
        return self.manager._create_new_zone(ids)

    def __sub__(self, other):
        "self - other <=> self.difference(other)"
        return self.difference(other)

    def symmetric_difference(self, other):
        "Return the symmetric difference of two zones as a new zone."
        self._valid_others(other)
        ids = self._block_ids.symmetric_difference(other._block_ids)
        return self.manager._create_new_zone(ids)

    def __xor__(self, other):
        "self ^ other <=> self.symmetric_difference(other)"
        return self.symmetric_difference(other)

    def issubset(self, other):
        "Return True if other contains self."
        self._valid_others(other)
        return self._block_ids.issubset(other._block_ids)

    def issuperset(self, other):
        "Return True if self contains other."
        self._valid_others(other)
        return self._block_ids.issuperset(other._block_ids)

    def isdisjoint(self, other):
        "Return True if two zones have a null intersection."
        self._valid_others(other)
        return self._block_ids.isdisjoint(other._block_ids)

    def __ge__(self, other):
        "self >= other <=> self.issuperset(other)"
        self._valid_others(other)
        return self._block_ids >= other._block_ids

    def __le__(self, other):
        "self <= other <=> self.issubset(other)"
        self._valid_others(other)
        return self._block_ids <= other._block_ids

    def __gt__(self, other):
        "self > other <=> self.issuperset(other) and self != other"
        self._valid_others(other)
        return self._block_ids > other._block_ids

    def __lt__(self, other):
        "self < other <=> self.issubset(other) and self != other"
        self._valid_others(other)
        return self._block_ids < other._block_ids


class ZoneManager(object):
    """ ZoneManager(size, check_consistency_default)

    A ZoneManager object works on a list of indices of length 'size',
    ie. equivalent to 'range(size)'. It produces Zone objects that
    represent subset of this list of indices.

    The method `manager.register(zone, name)` allow later accesses of zone
    through name in a read-only dictionnary way: `manager[name] -> zone`

    The methods keys(), values(), items(), iterkeys(), itervalues() and
    iteritems() are also avaible as for the built-in dict class.

    See also: Zone
    """

    # TODO: est-il pertinent d'ajouter une methode 'repack_partition()'
    #       dont le but serait de passer sur la partition maximal respectant
    #       self._zones ?

    def __init__(self, size, check_consistency_default=True):
        # TODO: décider la valeur par défaut de check_consistency_default
        # find optimal dtype
        try:
            uint_type = np.min_scalar_type(size)
        except AttributeError:  # numpy version < 1.6
            for uint_type in (np.uint8, np.uint16, np.uint32, np.uint64):
                if size == uint_type(size):
                    break
            else:
                raise TypeError("Invalid value for size (=%s)" % size)
        # create attributes
        self._dtype = np.dtype(uint_type)
        self._size = size
        # all instances of zones that the manager will manage will be of
        # instances of a specific subclass of Zone that is created here
        self._zone_cls = type(
            "Zone(%s)" % (self), (Zone,), {"_manager": self, "__doc__": Zone.__doc__}
        )
        self._zones = weakref.WeakValueDictionary()
        self._zone_counter = 0
        self._partition_callbacks = weakref.WeakKeyDictionary()
        self._namedzones = dict()
        self.check_consistency_default = check_consistency_default
        # create indices
        indices = np.arange(size, dtype=self._dtype)
        # init partition and the first zone 'whole'
        self._partition = (indices,)
        # create a first zone which is the whole universe
        self.register(self._create_new_zone([0]), "whole")

    def _find_zone_by_block_ids(self, ids):
        "return the managed zone matching `ids` else None"
        for zone in self._zones.values():
            if zone._block_ids == ids:
                return zone
        return None

    def _create_new_zone(self, ids):
        """Create, register and returns the new zone defined by given ids.
        As manager first checks if the corresponding zone has already been
        defined, no duplicate well be created."""
        ids = frozenset(ids)
        zone = self._find_zone_by_block_ids(ids)
        if zone is None:
            # create new zone
            zone = self._zone_cls()
            zone._block_ids = ids
            self._zones[self._zone_counter] = zone
            self._zone_counter += 1
        return zone

    def register(self, zone, name):
        "Register zone as name."
        if name in self._namedzones:
            raise KeyError("name (=%s) is already used." % name)
        if not isinstance(zone, self._zone_cls):
            raise ValueError("Wrong zone type.")
        self._namedzones[name] = zone

    def iterblocks(self):
        """Returns an iterator over blocks constitutiong the partition.
        When the block does not exist as a managed zone, a new managed zone
        will be created."""
        return (self._create_new_zone([i]) for i in range(len(self._partition)))

    @property
    def size(self):
        "The size of index list managed. (maximum index = size - 1)"
        return self._size

    @property
    def dtype(self):
        "The dtype used to store indices."
        return self._dtype

    @property
    def partition(self):
        "The current partition of whole indices."
        return self._partition

    def __getitem__(self, name):
        "manager[name] -> Return the zone registred as name."
        return self._namedzones[name]

    def __delitem__(self, name):
        "del manager[name] -> Unregistred the zone registred as name."
        del self._namedzones[name]

    def __getattr__(self, attr):
        """ Delegation of read-only dict behavior to internal data

        Methods accessible here are:
            'keys', 'values', 'items', 'iterkeys', 'itervalues',
            'iteritems'

        See the documentation of dict for further informations about
        this methods.
        """
        # Delegation of read-only dict behavior to _namedzones
        if attr in {"keys", "values", "items", "iterkeys", "itervalues", "iteritems"}:
            return getattr(self._namedzones, attr)
        else:
            raise AttributeError(
                "'%s' object has no attribute '%s'" % (self.__class__.__name__, attr)
            )

    def check_consistency(self, indices, flag):
        """ Check if indices are in the range of manager's indices.
        Raise an error if not.

        If flag is False, checking is skipped.

        If flag is None, check_consistency_default attribute indicates
        whether checking is skipped or not.
        """
        if not ((flag is None and self.check_consistency_default) or flag):
            return
        if len(indices) == 0:
            return
        if np.min(indices) < 0 or np.max(indices) >= self._size:
            raise ValueError("indices must be >= 0 and < manager.size")

    def build_from_mask(self, mask, name=None, check_consistency=False):
        """ Return a zone corresponding to mask.

        mask must 1d array of boolean the size of the manager.

        If name is specified, the zone is registred with this name.

        If check_consistency is specified, it override the value of
        attribute check_consistency_default.
        """
        assert len(mask) == self._size
        [indices] = np.where(mask)
        return self.build_from_indices(indices, name, check_consistency)

    def build_from_indices(self, indices, name=None, check_consistency=None):
        """ Return a zone corresponding to indices.

        If name is specified, the zone is registred with this name.

        If check_consistency is specified, it override the value of
        attribute check_consistency_default.
        """
        # check indices if needed
        self.check_consistency(indices, check_consistency)
        # build new partition, mapping with old partition and zone ids
        new_partition = list()
        map_ = list()
        ids = list()
        count = 0
        for part in self._partition:
            local_map = list()
            intersect = np.intersect1d(part, indices)
            difference = np.setdiff1d(part, indices)
            # TODO: is the following more efficient (cpu? memory?)
            # footprint = np.in1d(part, indices, assume_unique=True)
            # intersection = part[footprint]
            # difference = part[np.logical_not(footprint)]
            #
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
            assert local_map  # logic wants it is always true, but who knows.
            map_.append(local_map)
        # store new partition
        self._partition = tuple(new_partition)
        # update ids of all zones (must be done before creation of new zone)
        for zone in self._zones.values():
            new_block_ids = list()
            for i in zone._block_ids:
                new_block_ids.extend(map_[i])
            zone._block_ids = frozenset(new_block_ids)
        # call the partition callbacks
        self._invoke_partition_callbacks()
        zone = self._create_new_zone(ids)
        if name is not None:
            self.register(zone, name)
        return zone

    def add_partition_callback(self, func):
        """ add a callback invoked at new partition creation.

        WARNING: advanced usage only

        See also: remove_partition_callback()
        """
        if not callable(func):
            raise TypeError
        other = getattr(func, "__self__", None)
        if other is None:
            self._partition_callbacks[func] = None
        else:
            self._partition_callbacks[other] = func.__name__

    def remove_partition_callback(self, func):
        """ remove a callback invoked at new partition creation.

        WARNING: advanced usage only

        See also: add_partition_callback()
        """
        if not callable(func):
            raise TypeError
        other = getattr(func, "__self__", None)
        if other is None:
            self._partition_callbacks.pop(func, None)
        elif self._partition_callbacks.get(other, None) == func.__name__:
            self._partition_callbacks.pop(other, None)

    def _invoke_partition_callbacks(self):
        """ invoke the partition callbacks.

        ONLY USED BY : build_from_indices()

        See also: add_partition_callback(), remove_partition_callback()
        """
        for key, val in list(self._partition_callbacks.items()):
            (key if val is None else getattr(key, val))()


# TODO: zict[manager] devrait être équivalent à zict[manager['whole']]
#       aussi, manager['whole'] c'est nul, il faudrait un attribut


class Zict:
    """
    Zict objects are dict-like that only accept Zone objects as keys.

    A zict behaves like a built-in dict except that its keys
    follow set comparison (intersection, inclusion, ...) rather than
    hash comparison.

    The main consequences of this behavior are:

      - 'D[z1]' is licit if 'D[z2]' is too and 'z1.issubset(z2) == True'

      - if 'z1.isdisjoint(z2) == False' then: in the sequence
        > D[z1] = v1
        > D[z2] = v2
        z1 becomes an illicit key but (z1-z2) still is licit.

      - ...

    See also: Zone, ZoneManager
    """

    def __init__(self, E=(), **F):
        "Create an empty zict=Zict() then set zict.update(E,**F)"
        self._manager = None
        self._data = list()
        self.update(E, **F)

    @property
    def manager(self):
        "The manager of zones currently used as keys."
        return self._manager

    def items(self):
        "iterator on items"
        return iter(self._data)

    def keys(self):
        "iterator on keys"
        return iter(k for k, v in self._data)

    def values(self):
        "iterator on values"
        return iter(v for k, v in self._data)

    def get(self, key, default=None):
        "Return zict[key] if key in zict, else default."
        try:
            return self[key]
        except KeyError:
            return default

    def setdefault(self, key, default=None):
        """ Return zict.get(k, default),
            also set zict[k]= default if k not in zict"""
        try:
            return self[key]
        except KeyError:
            self[key] = default
            return default

    def pop(self, key, *default):
        """ Remove specified key and return the corresponding value.

        If key is not found, default is returned if given, otherwise
        KeyError is raised.
        """
        if len(default) > 1:
            TypeError("pop() expected at most 3 arguments, got %i" % (len(default) + 2))
        if default:
            value = self.get(key, default[0])
            try:
                del self[key]
            except KeyError:
                pass
            return value
        else:
            try:
                value = self[key]
                del self[key]
            except KeyError:
                raise
            return value

    def popitem(self):
        """
        D.popitem() -> (k, v), remove and return some (key, value) pair
        as a 2-tuple; but raise KeyError if D is empty.
        """
        if not self:
            raise KeyError("Zict object is empty")
        res = self._data.pop()
        if not self:
            self._manager = None
        return res

    def clear(self, zone=None):
        """Remove all subset of a given zone Z. The default is to clear the
        complete Zict instance. This works even if Z is not in Zict."""
        if self._manager is None:
            return
        if zone is None:
            zone = self._manager["whole"]
        keys = list(self.keys())
        for k in keys:
            if k <= zone:
                del self[k]
        if not self:
            self._manager = None

    def update(self, E, **F):
        """
        D.update(E, **F) -> None.  Update D from dict/iterable E and F.
        If E has a .keys() method, does:     for k in E: D[k] = E[k]
        If E lacks .keys() method, does:     for (k, v) in E: D[k] = v
        In either case, this is followed by: for k in F: D[k] = F[k]

        !!! WARNING !!!
            The result may strongly depends on the traversal order used
            for each arguments, especially when keys (ie. zones) are
            not dijoints. But this order can not always be predicted.
        """
        if hasattr(E, "keys"):
            if hasattr(E, "items"):
                for (k, v) in E.items():
                    self[k] = v
            else:
                for k in E:
                    self[k] = E[k]
        else:
            for (k, v) in E:
                self[k] = v

    def __iter__(self):
        return self.keys()

    def __bool__(self):
        return bool(self._data)

    def __contains__(self, key):
        "key in zict  -> True if 'key' is a licit key, else False."
        try:
            self[key]
        except KeyError:
            return False
        return True

    def _setitem(self, key, value):
        """ __setitem__() implementation

        Note: the separation with `__setitem__` avoid unwanted recursions.
        """
        # key must be Zone instance
        if not isinstance(key, Zone):
            raise KeyError(
                "key of %s instance must be a %s instance"
                % (self.__class__.__name__, Zone.__name__)
            )
        # ignore empty zone
        if not key:
            return
        # merge into key all keys mapping to value
        key = key.union(*(k for k, v in self.items() if v is value))
        # clean old keys intersecting key
        for z_old in [z for z in self if not z.isdisjoint(key)]:
            z_out = z_old - key
            self._setitem(z_out, self._popitem(z_old))
        # set item
        self._data.append((key, value))
        # if not yet manager, store it
        if self._manager is None:
            self._manager = key.manager

    def _getitem(self, key):
        for k, v in self._data:
            if k == key:
                return v
        raise KeyError()

    def _popitem(self, key):
        for i, (k, v) in enumerate(self._data):
            if k == key:
                del self._data[i]
                return v
        raise KeyError()

    def __setitem__(self, key, value):
        """zict[key] = value
        If key is an empty zone, affectation is ignored.
        """
        self._setitem(key, value)

    def __getitem__(self, key):
        """zict[key]
        If key is an empty zone, None is returned.
        """
        try:
            # search true key
            value = self._getitem(key)
        except KeyError:
            if not key:  # empty zone -> None
                return None
            # else search sub-key
            for value in (v for z, v in self.items() if key <= z):
                return value
        else:
            return value
        # nothing found -> error
        raise KeyError(key)

    def __delitem__(self, key):
        """del zict[key]
        Valid key may be subsets/subzones of actual internal dictionnary keys.
        Deleting such a subset will preserve the complementary subset and the
        associated value.

        You cannot delete zone providing a key that would be a superset of
        zone. Nevertheless you can directly change the value of two zones
        using their union as a new key.
        """
        if key not in self:
            raise KeyError(key)
        # object() compared (is/==) to anything will always be False,
        # thus associating 'key' to object() cleans the footprint of key.
        self[key] = object()
        # Now, 'key' is an actual key and can be deleted.
        self._popitem(key)
        # if self become empty, forget manager
        if not self:
            self._manager = None


class ZoneMap:
    """
    """

    def __init__(self, zone_manager):
        self._manager = zone_manager
        self._data = ()

    @property
    def manager(self):
        return self._manager

    def get(self, zone):
        ""
        if not zone:
            return None
        for part, value in self._data:
            if zone <= part:
                return value
        raise KeyError(zone)

    def set(self, zone, value):
        ""
        if not zone:
            return
        zone = zone.union(*(k for k, v in self._data if v is value))  # utile ?
        old_data = ((z - zone, v) for z, v in self._data)
        old_data = ((z, v) for z, v in old_data if z)
        self._data = ((zone, value), *old_data)

    def keys(self):
        return (z for z, _ in self._data)

    def items(self):
        return iter(self._data)

    def clean(self):
        self._data = tuple((z, v) for z, v in self._data if z)
