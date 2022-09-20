import numpy as np

# access to underlying mpi4py
from .mpi import MPI as mpi

# A transitory structure to hold various information
class DistributedSystem:
    @property
    def nb_comp_thermique(self):
        kernel = self.kernel
        assert kernel is not None
        NbCompThermique = kernel.model_number_of_components()
        if kernel.has_energy_transfer_enabled:
            NbCompThermique += 1
        return NbCompThermique

    def nb_item_rows(self, attr=None):
        NbCompThermique = self.nb_comp_thermique
        kernel = self.kernel
        part = kernel.part_info()
        if attr is None:
            attr = "nb"
        else:
            attr = "nb_" + attr
        # Rather be a loop over a list here...
        yield getattr(part.nodes, attr) * NbCompThermique
        yield getattr(part.fractures, attr) * NbCompThermique
        yield getattr(part.injectors, attr)
        yield getattr(part.producers, attr)
        yield getattr(part.mswell_nodes, attr) * NbCompThermique

    def items_and_rows(self):
        kernel = self.kernel
        NbCompThermique = self.nb_comp_thermique
        part = kernel.part_info()
        # Rather be a loop over a list here...
        yield part.nodes, NbCompThermique
        yield part.fractures, NbCompThermique
        yield part.injectors, 1
        yield part.producers, 1
        yield part.mswell_nodes, NbCompThermique

    def __init__(self, kernel):
        assert kernel is not None
        self.kernel = kernel

    @property
    def local_nb_rows(self):
        return np.sum([n for n in self.nb_item_rows()])

    @property
    def local_nb_cols(self):
        return np.sum([n for n in self.nb_item_rows("owns")])

    @property
    def global_nb_rows(self):
        return self.global_system_size[0]

    @property
    def global_nb_cols(self):
        return self.global_system_size[1]

    @property
    def local_system_size(self):
        return self.local_nb_rows, self.local_nb_cols

    @property
    def global_system_size(self):
        local_size = np.array(self.local_system_size, dtype=np.int)
        global_size = np.zeros(local_size.shape, dtype=local_size.dtype)
        mpi.COMM_WORLD.Allreduce(local_size, global_size, mpi.SUM)
        return global_size[0], global_size[1]  # convert to tuple

    def colnum(self):
        ColNum = np.zeros(self.local_nb_rows, dtype=np.int32)
        self.kernel.SyncPetsc_colnum(ColNum)
        ColNum -= 1  # Fortran -> C indexing
        return ColNum
