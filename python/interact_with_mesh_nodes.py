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

ComPASS.init_warmup_and_read_mesh(meshfile, logfile)

# FIXME: This would have to be renamed global mesh vertices
if comm.rank==0:
  vertices = np.array(ComPASS.get_vertices(), copy = False)
  print(vertices.shape)
  print(vertices[:5])
  print('Bounding box min corner:', vertices.min(axis=0))
  print('Bounding box max corner:', vertices.max(axis=0))
  con = ComPASS.get_connectivity()
  #for cell in con.NodebyCell:
  #  print(np.array(cell, copy=False))
  boundary_faces = np.array([fi for fi, face in enumerate(con.CellbyFace) if len(face)==1])
  boundary_faces_centers = np.array([
    vertices[
      np.array(con.NodebyFace[fi], copy=False) - 1 # fortran indexes start at 1
    ].mean(axis=0)
    for fi in boundary_faces
  ])
  print(boundary_faces_centers)
  
# FIXME: The following line is mandatory as it allocates structures
# that will be deallocated afterwards
ComPASS.init_phase2(outputdir)
# The following deallocate global mesh
ComPASS.finalize()
