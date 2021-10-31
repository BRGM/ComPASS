# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import sys
import numpy as np

import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
from ComPASS.newton import Newton, default_Newton
from ComPASS.linalg.factory import linear_solver
from ComPASS.utils.memory import MemStatus

# column permeability in m^2 (low permeability -> bigger time steps)
k_matrix = 1e-12
# column porosity
phi_matrix = 0.15
# bulk thermal conductivity in W/m/K
K_matrix = 2.0

ComPASS.set_output_directory_and_logfile(__file__)
simulation = ComPASS.load_eos("water2ph")

grid = ComPASS.Grid(shape=(10, 10, 10), extent=(2, 1, 1), origin=(0, 0, 0),)

simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
)

mem = MemStatus("loop", skip_first=True)
for _ in range(200):
    newton = default_Newton(simulation)
    mem.update(verbose=True)
print(mem)
