#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# We must load mpi4py first so that MPI is  initialized before calling PETSC_Initialize
from mpi4py import MPI

proc_rank = MPI.COMM_WORLD.rank
master_proc_rank = 0
is_on_master_proc = proc_rank == master_proc_rank


def communicator():
    return MPI.COMM_WORLD


def on_master_proc(f):
    def call(*args, **kwargs):
        if is_on_master_proc:
            return f(*args, **kwargs)

    return call


def master_print(*args, **kwargs):
    if is_on_master_proc:
        print(*args, **kwargs)


def synchronize():
    MPI.COMM_WORLD.Barrier()  # wait for every process to synchronize


def abort():
    MPI.COMM_WORLD.Abort()


if __name__ == "__main__":

    @on_master_proc
    def f():
        return 1

    print("on proc", proc_rank, "f()=", f())
