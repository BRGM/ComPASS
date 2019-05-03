#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
import doublet_utils
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop

import numpy as np

pres = 20. * MPa                            # initial reservoir pressure
Tres = degC2K( 70. )                        # initial reservoir temperature - convert Celsius to Kelvin degrees
Tinjection = degC2K( 30. )                  # injection temperature - convert Celsius to Kelvin degrees
Qm = 300. * ton / hour                      # production flowrate
k_matrix = 1E-13                            # matrix permeability in m^2
omega_matrix = 0.15                         # matrix porosity
K_matrix = 2                                # bulk thermal conductivity in W/m/K
omega_fracture = 0.5                        # fracture porosity
K_fracture = 2                              # bulk thermal conductivity in W/m/K
channel_width = 100.                        # fracture channel with 
channel_fracture_permeability = 1E-11       # permeability in m^2
background_fracture_permeability = 1E-12    # permeability in m^2

Lx, Ly, Lz = 3000., 2000., 100.
Ox, Oy, Oz = -1500., -1000., -1600.
nx, ny, nz = 61, 41, 3

ComPASS.load_eos('water2ph')
ComPASS.set_gravity(9.81)
ComPASS.set_fracture_thickness(1.)
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(
    shape = (nx, ny, nz),
    extent = (Lx, Ly, Lz),
    origin = (Ox, Oy, Oz),
)

def select_fractures():
    face_centers = ComPASS.compute_global_face_centers()
    zfaces = face_centers[:, 2]
    face_normals = ComPASS.compute_global_face_normals()
    ux, uy = face_normals[:, 0], face_normals[:, 1]
    dz = Lz / nz
    # select horizontal fault axis in the middle of the simulation domain
    zfrac = Oz + 0.5 * Lz
    fractures = (np.abs(zfaces - zfrac) < dz) & (ux==0) & (uy==0)
    print('Selecting', np.sum(fractures), 'faces as fractures')
    return fractures

def face_permeability():
    face_centers = ComPASS.compute_global_face_centers()
    nbfaces = face_centers.shape[0]
    xfc, yfc, zfc = [face_centers[:, col] for col in range(3)]
    interwell_distance = doublet_utils.interwell_distance(grid)
    xwell = Ox + 0.5 *(Lx - interwell_distance)
    def y_channel_centerline(x):
        return 0.25 * Ly * np.sin( 2 * np.pi * (x - xwell) / interwell_distance)
    in_channel = np.abs(y_channel_centerline(xfc) - yfc) < 0.5 * channel_width
    faceperm = np.empty(nbfaces, dtype=np.double)
    faceperm[:] = background_fracture_permeability
    faceperm[in_channel] = channel_fracture_permeability
    return faceperm[select_fractures()]

ComPASS.init(
    mesh = grid,
    wells = doublet_utils.make_wells_factory(grid),
    cell_porosity = omega_matrix,
    cell_permeability = k_matrix,
    cell_thermal_conductivity = K_matrix,
    fracture_faces = select_fractures,
    fracture_permeability = face_permeability,
    fracture_porosity = omega_fracture,
    fracture_thermal_conductivity = K_fracture,
)

doublet_utils.init_states(pres, Tres)

standard_loop(initial_timestep = day, final_time = 10 * year, output_period = 30 * day)
