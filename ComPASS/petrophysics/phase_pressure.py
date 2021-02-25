#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

# p^\alpha = pref + \Delta_{pref} where \Delta_{pref} = f(S^\alpha)
def no_differences(states, rocktypes, dpadS):
    # all phase pressures set to reference pressure
    states.pa[:] = states.p[:, None]
    # assert np.all(dpadS == 0)
    # dpadS.fill(0)


def set_phase_pressure_functions(simulation, f=None):
    if f is None:
        f = no_differences
    simulation.set_fill_phase_pressure_arrays(f)
