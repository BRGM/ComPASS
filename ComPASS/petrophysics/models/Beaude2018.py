import numpy as np
import numba

Pc0 = 2e5
Sg0 = 1.0 - 1.0e-2
Sl0 = 1.0 - Sg0
A = -Pc0 * np.log(Sl0) - (Pc0 / Sl0) * Sg0


@numba.vectorize(["float64(float64)"])
def Pc(Sg):
    if Sg < Sg0:
        return -Pc0 * np.log(1.0 - Sg)
    return Pc0 * Sg / Sl0 + A


@numba.vectorize(["float64(float64)"])
def dPcdS(Sg):
    if Sg < Sg0:
        return Pc0 / (1.0 - Sg)
    return Pc0 / Sl0


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 1, n)
    pc = Pc(Sg)
    # pc = np.hstack([laws[1][0](Sg), laws[2][0](Sg)])
    # Sg = np.hstack([Sg, Sg])
    plt.plot(1.0 - Sg, pc, ".")
    # plt.ylim(1e3, 1.2e5)
    plt.yscale("log")
    plt.xlabel("Sl")
    plt.ylabel("Beaude2018")
    plt.savefig("Pc_Beaude.png")
