import os
from mpi4py import MPI
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

ComPASS.init(meshfile, logfile, outputdir)
ComPASS.main_loop(0, outputdir)
# ComPASS.get_vertices()
ComPASS.finalize()
