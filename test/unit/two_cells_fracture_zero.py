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
from ComPASS.simulation_context import SimulationContext


# initial reservoir pressure
pL, pR = 2.0 * bar, 1.0 * bar
# initial reservoir temperature - convert Celsius degrees to Kelvin degrees
TL, TR = (
    degC2K(30),
    degC2K(20),
)
# column permeability in m^2 (low permeability -> bigger time steps)
k_matrix = 1e-12
# column porosity
phi_matrix = 0.15
# bulk thermal conductivity in W/m/K
K_matrix = 2.0

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(2, 1, 1), extent=(2, 1, 1), origin=(0, 0, 0),)


simulation.init(
    mesh=grid,
    cell_permeability=k_matrix,
    cell_porosity=phi_matrix,
    cell_thermal_conductivity=K_matrix,
    fracture_faces=np.array([0]),
    fracture_porosity=np.array([0.5]),
    fracture_permeability=np.array([1e-12]),
    fracture_thermal_conductivity=np.array([2]),
)
