#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

# p^\alpha = pref + \Delta_{pref} where \Delta_{pref} = f(S^\alpha)
def no_differences(states, rocktypes, pa, dpadS):
    # A lot of assertions that may slow down things...
    assert rocktypes.ndim == 1
    assert pa.ndim == 2
    assert dpadS.ndim == 2
    assert states.size() == rocktypes.shape[0] == pa.shape[0] == dpadS.shape[0]
    nb_phases = pa.shape[1]
    S = states.S
    assert S.shape[1] == nb_phases
    assert dpadS.shape[1] == nb_phases
    assert np.all((S >= 0) & (S <= 1))

    assert nb_phases == 2, "Should not be called if np==1"

    p = states.p
    pa[:, 0] = p
    pa[:, 1] = p
    assert np.all(dpadS == 0)
    # dpadS.fill(0)


def set_phase_pressure_functions(simulation, f=None):
    if f is None:
        f = no_differences
    simulation.set_fill_phase_pressure_arrays(f)
