# We must load mpi4py first so that MPI is  initialized before calling PETSC_Initialize
from mpi4py import MPI

proc_rank = MPI.COMM_WORLD.rank
is_on_master_proc = proc_rank==0

def on_master_proc(f):
    def call(*args, **kwargs):
        if is_on_master_proc:
            f(*args, **kwargs)
    return call

def synchronize():
    MPI.COMM_WORLD.Barrier() # wait for every process to synchronize
