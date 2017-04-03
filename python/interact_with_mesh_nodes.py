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

# ComPASS.init(meshfile, logfile, outputdir)
ComPASS.init_up_to_mesh(meshfile, logfile, outputdir)
# ComPASS.main_loop(0, outputdir)
# FIXME: This would have to be renamed global mesh vertices
if comm.rank==0:
  vertices = np.array(ComPASS.get_vertices(), copy = False)
  # Vertices array is in "Fortran order"
  print(vertices.shape)
  print(vertices[:,:5].transpose())
  print('Bounding box min corner:', vertices.min(axis=1))
  print('Bounding box max corner:', vertices.max(axis=1))
# The following deallocate global mesh
ComPASS.init_phase2(outputdir)
ComPASS.finalize()
