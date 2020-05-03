from collections import namedtuple, defaultdict
from itertools import chain

from .. import mpi
from .wells import get_well_data


WellFlowrates = namedtuple(
    "WellFlowrates", ["mass_flowrate", "energy_flowrates"], defaults=(None,) * 2
)


class WellDataConnection:
    def __init__(self, wid, source, mpi_tag=None, dest=None, data=None):
        self.well_id = wid
        self.source = source
        self.well_data = data
        # print(
        #     f"New connection {(wid, source, mpi_tag, dest)} on proc {mpi.proc_rank}"
        #     f"{'' if self.well_data is None else ' with data'}"
        # )
        assert not (self.well_data is None and source == mpi.proc_rank)
        assert dest is None or self.well_data is not None
        self.dest = []
        if dest is not None:
            self.dest = [proc for proc in dest if proc != mpi.proc_rank]
        self.mpi_tag = mpi_tag
        self._value = WellFlowrates()

    @property
    def value(self):
        return self._value

    def send(self):
        data = self.well_data
        if data is None:
            return
        flowrates = WellFlowrates(
            data.actual_mass_flowrate, data.actual_energy_flowrate
        )
        self._value = flowrates
        comm = mpi.communicator()
        for proc in self.dest:
            assert self.mpi_tag is not None
            comm.send(flowrates, dest=proc, tag=self.mpi_tag)

    def recv(self):
        if self.well_data is None:
            comm = mpi.communicator()
            assert self.mpi_tag is not None
            self._value = comm.recv(source=self.source, tag=self.mpi_tag)


class WellDataConnections:
    def __init__(self):
        self.wells = {}

    def add(self, wid, source, mpi_tag=None, dest=None, data=None):
        assert wid not in self.wells
        assert mpi_tag is None or mpi_tag not in [
            val.mpi_tag for val in self.wells.values()
        ]
        self.wells[wid] = WellDataConnection(wid, source, mpi_tag, dest, data)

    def synchronize(self):
        # WARNING: the two loops must NOT be merged
        #          FIRST we send all data
        #          SECOND we receive all data
        for well in self.wells.values():
            well.send()
        for well in self.wells.values():
            well.recv()

    def __getitem__(self, wid):
        return self.wells[wid].value


def _map_well_pairs(well_pairs, well_data_provider):
    comm = mpi.communicator()
    wells = [], [] # owns, ghosts
    for wid in set(chain(*well_pairs)):
        # FIXME: we could have an alternative way to find if a well is seen by this proc
        for k, own_only in enumerate((True, False)):
            data = well_data_provider(wid, own_only)
            if data is not None:
                wells[k].append(wid)
                break
    # print(f"Wells collected on proc {mpi.proc_rank}, owns: {wells[0]}, ghosts: {wells[1]}")
    return comm.gather(wells, root=mpi.master_proc_rank)


def _define_connections(well_pairs, well_data_provider):
    well_connections_tag = 12000  # FIXME: magic number
    mpi_tag = lambda well_id: well_connections_tag + well_id + 1
    comm = mpi.communicator()
    master = mpi.master_proc_rank
    wells_map = _map_well_pairs(well_pairs, well_data_provider)
    if mpi.is_on_master_proc:
        all_connections = [list() for _ in range(comm.size)]
        owning_proc = {}
        sent_wells = defaultdict(list)
        for source_well, target_well in well_pairs:
            need_sync = False
            for proc, (owns, ghosts) in enumerate(wells_map):
                all_wells = owns + ghosts
                if target_well in all_wells:
                    sent_wells[source_well].append(proc)
                if source_well in owns:
                    assert source_well not in owning_proc
                    owning_proc[source_well] = proc
        for well, dest in sent_wells.items():
            source_proc = owning_proc[well]
            tag = mpi_tag(well)
            dest = [proc for proc in dest if proc != source_proc]
            all_connections[source_proc].append((well, source_proc, tag, dest))
            for target_proc in dest:
                all_connections[target_proc].append((well, source_proc, tag))
        for proc, connections in enumerate(all_connections):
            if proc == master:
                continue
            comm.send(connections, dest=proc, tag=well_connections_tag)
        return all_connections[master]
    return comm.recv(source=master, tag=well_connections_tag)

def add_well_connections(well_pairs, connections, well_data_provider):
    for connection in _define_connections(well_pairs, well_data_provider):
        connections.add(*connection, data=well_data_provider(connection[0]))
