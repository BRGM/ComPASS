import numpy as np

Pc0 = 2e5
Sg0 = 1.0 - 1.0e-2
Sl0 = 1.0 - Sg0
A = -Pc0 * np.log(Sl0) - (Pc0 / Sl0) * Sg0


@np.vectorize
def Pc(Sg):
    if Sg < Sg0:
        return -Pc0 * np.log(1.0 - Sg)
    return Pc0 * Sg / Sl0 + A


@np.vectorize
def dPcdS(Sg):
    if Sg < Sg0:
        return Pc0 / (1.0 - Sg)
    return Pc0 / Sl0
