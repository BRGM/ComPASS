import ComPASS
import doublet_utils
from ComPASS.utils.units import *

ComPASS.set_output_directory_and_logfile(__file__)

def fractures_fractory(grid):
    def fractures():
        nz, Lz, Oz = grid.shape[2], grid.extent[2], grid.origin[2]
        face_centers = ComPASS.compute_face_centers()
        zfaces = face_centers[:, 2]
        face_normals = ComPASS.compute_face_normals()
        ux, uy = face_normals[:, 0], face_normals[:, 1]
        dz = Lz / nz
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = Oz + 0.5 * Lz
        return (np.abs(zfaces - zfrac) < dz) & (ux==0) & (uy==0) 
    return fractures

def face_permeability_factory(grid, channel_width=None):
    Lx, Ly, Lz = grid.extent
    if channel_width is None:
        channel_width = 0.1 * Ly
    def face_permeability():
        face_centers = ComPASS.compute_face_centers()
        nbfaces = face_centers.shape[0]
        xfc, yfc, zfc = [face_centers[:, col] for col in range(3)]
        interwell_distance = doublet_utils.interwell_distance(grid)
        xwell = grid.origin[0] + 0.5 *(Lx - interwell_distance)
        def y_channel_centerline(x):
            return 0.25 * Ly * np.sin( 2 * np.pi * (x - xwell) / interwell_distance)
        in_channel = np.abs(y_channel_centerline(xfc) - yfc) < 0.5 * channel_width
        faceperm = np.empty(nbfaces, dtype=np.double)
        faceperm[:] = 1E-13
        faceperm[in_channel] = 1E-11
        return faceperm
    return face_permeability

grid = ComPASS.Grid(
    shape = (31, 21, 3),
    extent = (3000., 2000., 100.),
    origin = (-1500., -1000., -1600.),
)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
    fracture_faces = fractures_fractory(grid),
    face_permeability = face_permeability_factory(grid),
)

@ComPASS.on_master_proc
def print_iteration_info():
    print()
    print('Time Step (iteration):', n)
    print('Current time: %.1f years' % (ComPASS.get_current_time() / year), ' -> final time:', final_time / year)
    print('Timestep: %.3f days' % (ComPASS.get_timestep() / day))
    
final_time = 30 * year
n = 0
output_frequency = 1 * year
t_output = 0
while ComPASS.get_current_time() <= final_time:
    n+= 1
    print_iteration_info()
    ComPASS.make_timestep()
    t = ComPASS.get_current_time()
    if t > t_output:
        ComPASS.output_visualization_files(n)
    # WARNING / CHECKME we may loose some outputs
    while (t_output < t):
        t_output = t_output + output_frequency
    ComPASS.timestep_summary()

ComPASS.finalize()
