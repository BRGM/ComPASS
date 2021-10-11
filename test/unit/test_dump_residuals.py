# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


import ComPASS.messages as messages
import two_cells_simulation


simulation = two_cells_simulation.make()


def dump_residual(newton_tick):
    cells_residuals = newton_tick.newton.convergence_scheme.residuals.cells
    print(
        f"Time loop iteration {newton_tick.timeloop_tick.iteration}",
        f"with dt {newton_tick.current_dt}",
        f"Newton iteration {newton_tick.iteration}",
    )
    print("-- Cell residuals --")
    print(cells_residuals)


current_time = simulation.standard_loop(
    initial_timestep=1, nitermax=2, newton_iteration_callbacks=[dump_residual],
)
