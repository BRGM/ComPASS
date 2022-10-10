#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import ComPASS
from ComPASS.utils.units import *
from ComPASS.petrophysics.models.Beaude2018 import Pc, dPcdS

pres = 50.0 * bar
Tres = degC2K(350.0)

simulation = ComPASS.load_physics("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(1, 1, 1))

simulation.init(
    mesh=grid,
    cell_porosity=0.5,
    cell_permeability=1,
    cell_thermal_conductivity=1,
)

simulation.set_liquid_capillary_pressure((Pc, dPcdS))

X0 = simulation.build_state(simulation.Context.gas, p=pres, T=Tres)
simulation.all_states().set(X0)

simulation.standard_loop(initial_timestep=1, nitermax=1)

cell_states = simulation.cell_states()
assert cell_states.size() == 1
X = cell_states[0]
Context = simulation.Context
assert Context(X.context) == Context.gas
assert X.pa[0] == X.p
assert X.pa[1] == X.p - Pc(X.S[0])
