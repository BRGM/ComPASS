"""

ArrayProperty stocke ses morceaux de array ainsi:
    prop._data[zone.manager][zone] -> (orig_zone, array)

Comment garder synchronisé (orig_zone, array) lors de manager._reduce_domain(z, idx) ?
    ni manager ni zone ne connait prop
    peut-on faire ça par callback ?

    mettre à jour array n'est possible que s'il ne change pas de taille.
    (orig_zone, array) -> ZoneArrayItem(orig_zone, array) ?
    permet de remplacer array par un nouvel objet si besoin.
    le besoin de synchro dépend de la seule existence de l'item.


"""

import weakref
import numpy as np
from index_zoning import ZoneManager


class ZoneArrayItem():
    def __init__(self, zone, array, is_array=False):
        self.zone = zone
        self.array = array
        if is_array:
            # zone.manager.add_reduce_callback(zone, self._reduce_array_cb())
            add = getattr(zone.manager, 'add_reduce_callback', None)
            if add:
                add(zone, self._reduce_array_cb())

    def __iter__(self):
        return iter((self.zone, self.array))

    def _reduce_array_cb(self):
        ref = weakref.ref(self)

        def cb(new_indices, ref=ref):
            """ change array par array[new_indices]

            return True si doit être nettoyé

            """
            self = ref()
            if self is None:
                return True
            self.array = self.array[new_indices]

        return cb


class ZoneManager_bis(ZoneManager):
    """

    reduce_callback :
      - fonction qui accepte un array d'indices new_indices en argument
      - appelé quand la zone associée est réduite
      - new_indices décrit les indices de zone.content() conservés par la
        réduction ainsi que leur nouvel ordre
      - reduce_callback doit mettre à jour les données concernées
      - retourne un booléen :
          True -> oublier le callback
          False -> le callback sera rappelé à la prochaine réduction

    """
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self._reduce_callbacks = {}

    @property
    def void(self):
        return self.build_from_indices(())

    def add_reduce_callback(self, zone, cb):
        assert zone.manager is self
        # self._reduce_callbacks.setdefault(zone, set()).add(cb)
        self._reduce_callbacks.setdefault(id(zone), set()).add(cb)

    def get_reduce_callbacks(self, zone):
        return self._reduce_callbacks.get(id(zone), set())

    def _reduce_domain(self, indices):
        indices = np.array(indices, copy=False)
        order = np.argsort(indices)
        indices = indices[order]
        new_domain = self.build_from_indices(indices)

    # def _reduce_domain(self, new_domain, new_content):
        """

        change l'état du manager:
         - seul les indices désignés par zone sont gardés
         - ces indices sont renumérotés selon new_content

        """
        # assert isinstance(new_domain, self._zone_cls)

        ids = set(new_domain._block_ids)
        new_ids = {n: i for i, n in enumerate(sorted(ids))}

        for zone in self._zones.values():
            # update les dépendences à zone via callback
            callbacks = self.get_reduce_callbacks(zone)
            if callbacks:
                _, jj, ii = np.intersect1d(
                    indices, zone.content(),
                    return_indices=True,
                    assume_unique=True,
                )
                ii = ii[np.argsort(order[jj])]
                to_remove = [cb for cb in callbacks if cb(ii)]
                callbacks.difference_update(to_remove)

            # reduction de la définition de zone
            block_ids = ids & zone._block_ids
            block_ids = map(new_ids.get, block_ids)
            zone._block_ids = frozenset(block_ids)
            zone._content_ref = None

        # update self._partition
        self._size = len(indices)
        partition = [self._partition[i] for i in ids]
        partition = [
            order[np.intersect1d(p, indices, return_indices=1)[2]]
            for p in partition
        ]
        # new_domain.content()[new_content]
        self._partition = tuple(partition)


class GrandManitou:
    def variables(self):
        """ renvoie le namespace des propriétés gérées

        certaines propriétés (ex locations) sont pré-remplies

        l'objet est aussi iterable:
            list(namespace) -> [(name, property), ...]

        """

    def zones(self):
        """ renvoie le namespace des domaines (zones) gérées
        """

    def distribute(self, *args):
        # args: partition des cells, nodes, faces, ...

        local_indices = ...  # dépend de MPI proc number
        self._reduce_domain(local_indices)
        ...

    def _reduce_domain(self, local_indices):
        assert isinstance(local_indices, ArrayProperty)
        managers = set(local_indices._data)
        properties = {p: n for n, p in self.variables()}
        assert managers.issubset(z.manager for _, z in self.zones())
        assert all(managers.issuperset(prop._data) for prop in properties)

        domains = {}
        for _, zone in self.zones():
            domains.setdefault(zone.manager, zone)
            domains[zone.manager] |= zone

        for manager, domain in domains.items():
            parts = local_indices._get_data(manager).keys()
            new_domain = manager.void.union(*parts)
            assert new_domain <= domain
            # manager._reduce_domain(new_domain, local_indices[new_domain])  # FIXME
            manager._reduce_domain(new_domain, local_indices)

#         # si besoin, nettoyer les ZoneMap
#         # un callback entre ZoneMap et ZoneManager est-il mieux ?
#         for prop in properties:
#             prop.clean()
