#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS

nz = 100
H = 10.0

simulation = ComPASS.load_eos("water2ph")
ComPASS.set_output_directory_and_logfile(__file__)

grid = ComPASS.Grid(shape=(1, 1, nz), extent=(10.0, 10.0, H), origin=(-5, -5, -H),)

simulation.init(
    mesh=grid, cell_permeability=1.0, cell_porosity=0.5, cell_thermal_conductivity=1.0,
)

X0 = simulation.build_state(simulation.Context.liquid, p=3.14, T=1.0)
states = simulation.all_states()
states.set(X0)
assert np.all(states.p == 3.14)
mask = np.arange(states.size()) % 2 == 0
X1 = simulation.build_state(X0)
X1.p = 42
assert X0.p == 3.14
states.set(mask, X1)
for k, pk in enumerate(states.p):
    assert (k % 2 == 0 and pk == 42) or pk == 3.14
    pass
