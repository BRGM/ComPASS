# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 18:40:41 2017

@author: lopez
"""

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
#    data = {'a': 7, 'b': 3.14}
    data = 'ma jolie chaine'
    req = comm.isend(data, dest=1, tag=11)
    req.wait()
    data = None
elif rank == 1:
    data = None
    req = comm.irecv(source=0, tag=11)
    data = req.wait()

print('On', rank, 'data =', data)    
