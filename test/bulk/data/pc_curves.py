import numpy as np

Pc0 = 2.0e5
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
    if Sg < S0:
        return Pc0 / (1.0 - Sg)
    return Pc0 / Sl0


def get():
    return Pc, dPcdS


def plot_curves():
    import matplotlib.pylab as plt

    plt.clf()
    S = np.linspace(0, 1)
    plt.plot(100 * S, Pc(S) / 1e5)
    plt.xlabel("gas saturation (%)")
    plt.ylabel("capillary pressure $p_c = p_{air}-p_{water}$ (bar)")
    plt.savefig("pc.png")


if __name__ == "__main__":
    plot_curves()
