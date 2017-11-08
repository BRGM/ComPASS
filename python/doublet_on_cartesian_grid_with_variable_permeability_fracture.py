import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

import numpy as np

ComPASS.load_eos('water2ph')

def fractures_factory(grid):
    def fractures():
        nz, Lz, Oz = grid.shape[2], grid.extent[2], grid.origin[2]
        face_centers = ComPASS.compute_global_face_centers()
        zfaces = face_centers[:, 2]
        face_normals = ComPASS.compute_global_face_normals()
        ux, uy = face_normals[:, 0], face_normals[:, 1]
        dz = Lz / nz
        # select horizontal fault axis in the middle of the simulation domain
        zfrac = Oz + 0.5 * Lz
        where = (np.abs(zfaces - zfrac) < dz) & (ux==0) & (uy==0)
        print('Selecting', np.sum(where), 'faces as fractures')
        return where 
    return fractures

def face_permeability_factory(grid, channel_width=None):
    Lx, Ly, Lz = grid.extent
    if channel_width is None:
        channel_width = 0.1 * Ly
    def face_permeability():
        face_centers = ComPASS.compute_global_face_centers()
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

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    grid = grid,
    wells = doublet_utils.make_wells_factory(grid),
    fracture_faces = fractures_factory(grid),
    faces_permeability = face_permeability_factory(grid),
)

#standard_loop(final_time = 30 * year, output_frequency = year)
