from collections import namedtuple, defaultdict
from itertools import chain
from functools import partial

from .. import mpi


WellHead = namedtuple(
    "WellHead",
    ["molar_flowrate", "energy_flowrate", "pressure", "temperature"],
    defaults=(None,) * 4,
)


class WellDataConnection:
    def __init__(self, wid, source, mpi_tag=None, dest=None, data=None, provider=None):
        assert data is None or provider is None
        self.well_id = wid
        self.source = source
        self._well_data = provider or (lambda: data)
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
        self._value = WellHead()

    @property
    def well_data(self):
        return self._well_data()

    @property
    def value(self):
        return self._value

    def send(self):
        data = self.well_data
        if data is None:
            return
        wellhead = WellHead(
            data.molar_flowrate,
            data.energy_flowrate,
            data.pressure,
            data.temperature,
        )
        self._value = wellhead
        comm = mpi.communicator()
        for proc in self.dest:
            assert self.mpi_tag is not None
            comm.send(wellhead, dest=proc, tag=self.mpi_tag)

    def recv(self):
        if self.well_data is None:
            comm = mpi.communicator()
            assert self.mpi_tag is not None
            self._value = comm.recv(source=self.source, tag=self.mpi_tag)


class WellDataConnections:
    def __init__(self):
        self.wells = {}

    def add(self, wid, source, mpi_tag=None, dest=None, data=None, provider=None):
        assert wid not in self.wells
        assert mpi_tag is None or mpi_tag not in [
            val.mpi_tag for val in self.wells.values()
        ]
        self.wells[wid] = WellDataConnection(
            wid, source, mpi_tag, dest, data=data, provider=provider
        )

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


def _map_wells_distribution(wells, well_data_provider):
    distribution = [], []  # owns, ghosts
    for wid in set(wells):
        # FIXME: we could have an alternative way to find if a well is seen by this proc
        for k, own_only in enumerate((True, False)):
            data = well_data_provider(wid, own_only)
            if data is not None:
                distribution[k].append(wid)
                break
    # print(f"Wells collected on proc {mpi.proc_rank}, owns: {wells[0]}, ghosts: {wells[1]}")
    return mpi.communicator().gather(distribution, root=mpi.master_proc_rank)


def _mpi_well_connections_tag():
    return 12000  # FIXME: magic number


def _mpi_well_tag(well_id):
    return _mpi_well_connections_tag() + well_id + 1


def _send_connections(owning_proc, sent_wells):
    assert mpi.is_on_master_proc
    master = mpi.master_proc_rank
    comm = mpi.communicator()
    all_connections = [list() for _ in range(comm.size)]
    for well, dest in sent_wells.items():
        source_proc = owning_proc[well]
        tag = _mpi_well_tag(well)
        dest = [proc for proc in dest if proc != source_proc]
        all_connections[source_proc].append((well, source_proc, tag, dest))
        for target_proc in dest:
            all_connections[target_proc].append((well, source_proc, tag))
    for proc, connections in enumerate(all_connections):
        if proc == master:
            continue
        comm.send(connections, dest=proc, tag=_mpi_well_connections_tag())
    return all_connections[master]


def _define_connections(well_data_provider, well_pairs=None, proc_requests=None):
    """
    :param well_pairs: a sequence of well pair to be chained
    :param proc_requests: a sequence of pair (proc, list of wells to make available)
    """
    wells = set()
    if well_pairs is not None:
        wells.update(chain(*well_pairs))
    else:
        well_pairs = []
    if proc_requests is not None:
        for _, well_list in proc_requests:
            wells.update(well_list)
    else:
        proc_requests = []
    # the following line MUST be executed by all procs
    wells_map = _map_wells_distribution(wells, well_data_provider)
    if mpi.is_on_master_proc:
        owning_proc = {}
        for proc, (owns, _) in enumerate(wells_map):
            for well in owns:
                owning_proc[well] = proc
        sent = defaultdict(set)
        for proc, (owns, ghosts) in enumerate(wells_map):
            all_wells = owns + ghosts
            for source_well, target_well in well_pairs:
                if target_well in all_wells:
                    sent[source_well].add(proc)
        for target_proc, well_list in proc_requests:
            for well in well_list:
                sent[well].add(target_proc)
        return _send_connections(owning_proc, sent)
    return mpi.communicator().recv(
        source=mpi.master_proc_rank, tag=_mpi_well_connections_tag()
    )


def add_well_connections(
    connections, well_data_provider, well_pairs=None, proc_requests=None
):
    """
    :param well_data_provider: a function that will consume a well id and return some compliant data
    :param well_pairs: a sequence of well pair to be chained
    :param proc_requests: a sequence of pair (proc, list of wells to make available)
    """
    new_connections = _define_connections(
        well_data_provider, well_pairs=well_pairs, proc_requests=proc_requests
    )
    for connection in new_connections:
        # FIXME: the following does not work if well_data_provider returns a copy
        #        instead of a reference to the actual data
        # connections.add(*connection, data=well_data_provider(connection[0]))
        connections.add(
            *connection, provider=partial(well_data_provider, connection[0])
        )
