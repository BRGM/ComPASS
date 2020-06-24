import sys
from mpi4py import MPI

n = 4

# comm = MPI.COMM_WORLD

comm = MPI.COMM_SELF.Spawn(sys.executable, args=[__file__], maxprocs=n)

rank = comm.Get_rank()

data = comm.gather(rank, root=0)
print(rank, data)
