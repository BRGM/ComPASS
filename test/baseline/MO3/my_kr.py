import numpy as np

# Cf. Special Panel on Geothermal Model Intercomparison Study
# Sixth Stanford Annual Workshop on Geothermal Reservoir Engineering
# December 16-18, 1980
# Problem 3: 2-D Well Test


def kr_functions(states, rocktypes, kr, dkrdS):
    assert kr.shape[1] == 2, "Should not be called if np==1"
    kr[:, 0] = 1
    kr[:, 1] = 0
    dkrdS[...] = 0
