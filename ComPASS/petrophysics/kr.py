#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np


def no_interactions(states, rocktypes, kr, dkrdS):
    # A lot of assertions that may slow down things...
    # assert rocktypes.ndim == 1
    # assert kr.ndim == 2
    # assert dkrdS.ndim == 3
    # assert states.size() == rocktypes.shape[0] == kr.shape[0] == dkrdS.shape[0]
    nb_phases = kr.shape[1]
    S = states.S
    # assert S.shape[1] == nb_phases
    # assert dkrdS.shape[1] == dkrdS.shape[2] == nb_phases
    # assert np.all((S >= 0) & (S <= 1))

    assert nb_phases == 2, "Should not be called if np==1"
    kr[...] = S
    dkrdS.fill(0)
    for k in range(nb_phases):
        dkrdS[:, k, k].fill(1)


def S2(states, rocktypes, kr, dkrdS):
    nb_phases = kr.shape[1]
    S = states.S
    assert nb_phases == 2, "Should not be called if np==1"
    kr[...] = S**2
    dkrdS.fill(0)
    for k in range(nb_phases):
        dkrdS[:, k, k] = 2 * S[:, k]


def set_kr_functions(simulation, f=None):
    if f is None:
        f = S2
    simulation.set_fill_kr_arrays(f)
