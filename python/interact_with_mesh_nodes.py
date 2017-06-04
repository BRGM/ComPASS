import os
from mpi4py import MPI
import numpy as np
import ComPASS

comm = MPI.COMM_WORLD

compass_package_dir, tail = os.path.split(ComPASS.__file__)
compass_python_dir, tail = os.path.split(compass_package_dir)
outputdir = os.path.join(compass_python_dir, 'tests',
                         'output-' + os.path.splitext(os.path.basename(__file__))[0])
outputdir = os.path.abspath(outputdir)

# master proc manages directory creation
if comm.rank==0:
    if not os.path.exists(outputdir):
      os.makedirs(outputdir)
comm.Barrier() # wait for every process to synchronize

logfile = os.path.join(outputdir, 'mytest.log')

#ComPASS.init_warmup_and_read_mesh(meshfile, logfile)
ComPASS.init_warmup(logfile)
if comm.rank==0:
    origin = (-1500., -1000., -1600.)
    extent = (3000., 2000., 100.)
    shape = (3, 2, 1)
    wells = [] # no wells
    ComPASS.build_grid(shape = shape, origin = origin, extent = extent)
    ComPASS.set_well_geometries(wells)
    ComPASS.global_mesh_make_post_read()
    ComPASS.set_well_data(wells)
    ComPASS.compute_well_indices()

# FIXME: This would have to be renamed global mesh vertices
if comm.rank==0:
  vertices = ComPASS.global_vertices()
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

