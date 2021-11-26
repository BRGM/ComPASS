# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


import two_cells_simulation
from ComPASS.timeloops import TimeStepManager

simulation = two_cells_simulation.make()

tsm = TimeStepManager(initial_timestep=1, increase_factor=2)

t = simulation.standard_loop(nitermax=2, time_step_manager=tsm)

# reset timestep
tsm.current_step = 1

simulation.standard_loop(
    initial_time=t, nitermax=2, time_step_manager=tsm, reset_iteration_counter=True
)
