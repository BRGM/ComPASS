# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


import ComPASS.messages as messages
import two_cells_simulation


simulation = two_cells_simulation.make()


def send_warning(tick):
    if tick.iteration == 1:
        messages.warning("Warning after one iteration")


def send_error(tick):
    if tick.iteration == 2:
        messages.error("Error after two iterations")


current_time = simulation.standard_loop(
    initial_timestep=1,
    nitermax=2,
    output_period=1,
    iteration_callbacks=[send_warning, send_error],
)
