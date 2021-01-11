#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np
import sys
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop, TimeStepManager
from ComPASS.linalg.factory import linear_solver
from ComPASS.newton import Newton


Sgres = 0.9
if len(sys.argv) > 1:
    kres = float(sys.argv[1])  # reservoir permeability in m^2
else:
    kres = 5e-17  # reservoir permeability in m^2
phires = 0.15  # reservoir porosity
K = 2.0  # bulk thermal conductivity in W/m/K
H = 2000.0  # column height
nx, ny, nz = 1, 1, 100  # discretization
Lx = (H / nz) * nx
Ly = (H / nz) * ny
Tres = degC2K(295)
Ttop = degC2K(15)
caprock_thickness = 2 * km
bottom_heat_flux = K * (Tres - Ttop) / caprock_thickness  # W/m2
gravity = 9.81
final_time = 1e5 * year

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(
    f"heatpipe-k{kres:g}.py", process_case_name=False
)

simulation.set_gravity(gravity)

grid = ComPASS.Grid(
    shape=(nx, ny, nz), extent=(Lx, Ly, H), origin=(-0.5 * Lx, -0.5 * Ly, -H),
)

simulation.init(
    mesh=grid,
    cell_permeability=kres,
    cell_porosity=phires,
    cell_thermal_conductivity=K,
)

Xres = simulation.build_state(simulation.Context.diphasic, T=Tres, Sg=Sgres)

simulation.all_states().set(Xres)

print(f"Bottom heat flux: {bottom_heat_flux} W/m^2")

face_centers = simulation.face_centers()
Neumann = ComPASS.NeumannBC()
Neumann.heat_flux = bottom_heat_flux
simulation.set_Neumann_faces(face_centers[:, 2] <= -H, Neumann)  # input flux
Neumann.heat_flux = -bottom_heat_flux
simulation.set_Neumann_faces(face_centers[:, 2] >= 0, Neumann)  # output flux

lsolver = linear_solver(simulation, direct=True)
newton = Newton(simulation, 1e-5, 8, lsolver)
tsmger = TimeStepManager(1 * day)

simulation.standard_loop(
    final_time=final_time,
    newton=newton,
    time_step_manager=tsmger,
    timeloop_statistics=True,
)

simulation.postprocess(convert_temperature=True)
