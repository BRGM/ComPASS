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

# Simulation domain
origin = (-1500., -1000., -1600.)
extent = (3000., 2000., 100.)
shape = (3, 2, 1)

ComPASS.init_warmup(logfile)
if comm.rank==0:
    ComPASS.build_grid(shape = shape, origin = origin, extent = extent)
    ComPASS.set_well_geometries([])
    ComPASS.global_mesh_mesh_bounding_box()
    ComPASS.global_mesh_compute_all_connectivies()
    fractures_centers = ComPASS.compute_face_centers()
    yfc = fractures_centers[:, 1]
    print("Setting", np.sum(yfc==0), "fracture faces.")
    idfaces = ComPASS.get_id_faces()
    assert fractures_centers.shape[0]==idfaces.shape[0]
    # this will set directly fracture faces as the memory array is the one used by Fortran
    idfaces[yfc==0] = -2
    ComPASS.global_mesh_set_frac() # this will collect faces with flag -2 as fracture faces
    ComPASS.global_mesh_node_of_frac()
    ComPASS.global_mesh_set_dir_BC()
    ComPASS.global_mesh_frac_by_node()
    ComPASS.global_mesh_make_post_read_set_poroperm()
    fracperm = ComPASS.get_fracture_permeability()
    print("Number of faces:", idfaces.shape[0])
    print("Number of fractures faces:", np.sum(yfc==0))
    assert fracperm.shape[0]==idfaces.shape[0]
    ComPASS.global_mesh_make_post_read_well_connectivity_and_ip()
    ComPASS.set_well_data([])
    ComPASS.compute_well_indices()

ComPASS.init_phase2(outputdir)

#ComPASS.main_loop(0, outputdir)
ComPASS.finalize()
