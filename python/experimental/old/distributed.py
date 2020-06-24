from mpi4py import MPI
from index_zoning import ZoneManager as SeqZoneManager


# class Distributed():
#     @property
#     def __comm__(self):
#         return MPI.COMM_WORLD
#     @property
#     def __rank__(self):
#         return self.__comm__.Get_rank()


class ZoneManager(SeqZoneManager):
    def __init__(self, size, check_consistency_default=True, comm=None):
        super().__init__(size, check_consistency_default)
        if comm is None:
            comm = MPI.COMM_WORLD
        self._comm = comm

    @property
    def comm(self):
        return self._comm

    def distribute(self, indices_by_rank):
        comm = self.comm
        rank = comm.Get_rank()
        print(rank)
        for zid, zone in self._zones.items():
            all_zid = comm.gather(zid, root=0)
            assert all_zid is None or set(all_zid) == {zid}
            if rank == 0:
                self._distribute_master(indices_by_rank, zid, zone)
            else:
                self._distribute_slave(zid, zone)

    def _distribute_master(self, indices_by_rank, zid, zone):
        pass

    def _distribute_slave(self, zid, zone):
        pass
