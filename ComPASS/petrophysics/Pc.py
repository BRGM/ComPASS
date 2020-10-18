#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np


def no_capillary_pressure(states, rocktypes, Pc, dPcdS):
    # A lot of assertions that may slow down things...
    assert rocktypes.ndim == 1
    assert Pc.ndim == 2
    assert dPcdS.ndim == 3
    assert states.size() == rocktypes.shape[0] == Pc.shape[0] == dPcdS.shape[0]
    nb_phases = Pc.shape[1]
    S = states.S
    assert S.shape[1] == nb_phases
    assert dPcdS.shape[1] == dPcdS.shape[2] == nb_phases
    assert np.all((S >= 0) & (S <= 1))

    assert nb_phases == 2, "Should not be called if np==1"
    Pc.fill(0)
    dPcdS.fill(0)


def set_Pc_functions(simulation, f=None):
    if f is None:
        f = no_capillary_pressure
    simulation.set_fill_Pc_arrays(f)
