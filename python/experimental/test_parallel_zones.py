import numpy as np
from mpi4py import MPI
from zones import ZoneManager, Parallel

import pytest


class ParallelManager(Parallel, ZoneManager):
    pass


def test_one_index_par_proc():
    comm = MPI.COMM_WORLD

    mng = ParallelManager(0)
    mng.reinit([np.array([0])], {})
    val = np.array([comm.rank])

    a = mng.from_mask(val % 2)
    b = mng.from_mask(val % 3)
    c = mng.from_mask(val % 4)

    assert len(a) in (0, 1)
    assert len(b) in (0, 1)

    expected = [0] if comm.rank % 2 else []
    assert list(a.content) == expected

    expected = [0] if comm.rank % 3 else []
    assert list(b.content) == expected

    expected = [0] if comm.rank % 2 and comm.rank % 3 else []
    assert list((a & b).content) == expected

    expected = [0] if comm.rank % 2 or comm.rank % 3 else []
    assert list((a | b).content) == expected

    expected = [0] if comm.rank % 2 and not (comm.rank % 3) else []
    assert list((a - b).content) == expected

    assert a <= c
    assert a & b <= b
