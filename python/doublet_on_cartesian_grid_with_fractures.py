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
shape = (31, 21, 11)

# units
tons_per_hour = 1E3 / 3600 # 1000 kg / 3600 s

def make_well(xy):
    vertices = ComPASS.global_vertices()
    x, y, z = (vertices[:, col] for col in range(3))
    x_well = x[np.argmin(np.abs(x - xy[0]))]
    y_well = y[np.argmin(np.abs(y - xy[1]))]
    well_nodes = np.nonzero((x == x_well) & (y == y_well))[0]
    # CHECKME: What is the expected order for well nodes?
    well_nodes = well_nodes[np.argsort(z[well_nodes])]
    well = ComPASS.Well()
    well.geometry.radius = 0.1
    segments = np.transpose(np.vstack([well_nodes[:-1], well_nodes[1:]]))
    well.geometry.add_segments(segments + 1) # Fortran indices start at 1
    return well

def make_wells(interwell_distance=None):
    if interwell_distance and comm.rank==0:
        producer = make_well((-0.5 * interwell_distance, 0))
        #producer.operate_on_flowrate = 300, 1E5
        producer.operate_on_pressure = 10E6, 100 * tons_per_hour
        producer.produce()
        injector = make_well((0.5 * interwell_distance, 0))
        #injector.operate_on_flowrate = 300, 30E6
        injector.operate_on_pressure = 30E6, 100 * tons_per_hour
        #injector.inject(273.15 + 30) # 30°C in K
        injector.inject(333) # FIXME
        return [producer, injector]

def select_fractures():
    fractures_centers = ComPASS.compute_face_centers()
    dz = extent[2] / shape[2]
    # select horizontal fault axis in the middle of the simulation domain
    zfrac = origin[2] + 0.5 * extent[2]
    return np.abs(fractures_centers[:, 2] - zfrac) < 0.5 * dz

# The following calls are equivalent to ComPASS.init(meshfile, logfile, outputdir)
ComPASS.init_warmup(logfile)
if comm.rank==0:
    ComPASS.build_grid(shape = shape, origin = origin, extent = extent)
    wells = make_wells(interwell_distance = extent[0]/3)
    ComPASS.set_well_geometries(wells)
    ComPASS.global_mesh_mesh_bounding_box()
    ComPASS.global_mesh_compute_all_connectivies()
    ComPASS.set_fractures(select_fractures())
    ComPASS.global_mesh_node_of_frac()
    ComPASS.global_mesh_set_dir_BC()
    ComPASS.global_mesh_frac_by_node()
    ComPASS.global_mesh_make_post_read_set_poroperm()
    ComPASS.global_mesh_make_post_read_well_connectivity_and_ip()
    ComPASS.set_well_data(wells)
    ComPASS.compute_well_indices()

ComPASS.init_phase2(outputdir)

ComPASS.main_loop(0, outputdir)
ComPASS.finalize()
