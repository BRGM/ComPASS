#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

# create a vtu file with a cartesian mesh, then read it with vtk2RawMesh

import numpy as np
import ComPASS
import ComPASS.debug_utils
import ComPASS.io.mesh as io
from ComPASS.utils.vtk import vtk2RawMesh


simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

Lx, Ly, Lz = 3000.0, 2000.0, 100.0
Ox, Oy, Oz = -1500.0, -800.0, -1600.0
nx, ny, nz = 6, 7, 3

grid = ComPASS.Grid(
    shape=(nx, ny, nz),
    extent=(Lx, Ly, Lz),
    origin=(Ox, Oy, Oz),
)

omega_reservoir = 0.15  # reservoir porosity
k_reservoir = 1e-12  # reservoir permeability in m^2
K_reservoir = 2  # bulk thermal conductivity in W/m/K
simulation.init(
    mesh=grid,
    cell_porosity=omega_reservoir,
    cell_permeability=k_reservoir,
    cell_thermal_conductivity=K_reservoir,
)

io.write_mesh(simulation, "mesh_alone")
raw_mesh = vtk2RawMesh("./mesh_alone.vtu")
assert np.shape(raw_mesh.get_vertices()) == ((nx + 1) * (ny + 1) * (nz + 1), 3)
