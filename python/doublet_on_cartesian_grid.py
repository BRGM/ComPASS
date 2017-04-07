import os
from mpi4py import MPI
import numpy as np
import ComPASS

comm = MPI.COMM_WORLD

head, tail = os.path.split(ComPASS.__file__)
compassdir, tail = os.path.split(head)
assert tail=='python'
dirname = os.path.join(compassdir,'tests')
meshdir = os.path.join(dirname, 'meshes')

meshfile = os.path.join(meshdir, 'cartesian.msh')

outputdir = os.path.abspath('./mytest')
if not os.path.exists(outputdir):
  os.makedirs(outputdir)

logfile = os.path.join(outputdir, 'mytest.log')

origin = (-1500., -1000., -1600.)
extent = (3000., 2000., 100.)
shape = (30, 20, 10)

# The following calls are equivalent to ComPASS.init(meshfile, logfile, outputdir)
ComPASS.init_warmup(logfile)
if comm.rank==0:
  ComPASS.init_build_grid(*origin, *extent, *shape)
ComPASS.init_phase2(outputdir)

ComPASS.main_loop(0, outputdir)
ComPASS.finalize()
