import time
import math
import numpy as np

try:
    from scipy.special import erf, erfc
except ImportError:
    # this is not really efficient just a workaround when scipy is missing
    # good for government work
    erf = np.vectorize(math.erf)
    erfc = np.vectorize(math.erfc)

import ComPASS.utils.mpl_backends as mpl_backends

plt = mpl_backends.import_pyplot(False)

# la température
def exact_sol(x, t):
    return (
        1
        / 2
        * Cm
        * np.exp(gamma * t)
        * (
            np.exp(x * (U - P) / (2 * Dx))
            * erfc((Rd * x - P * t) / (2 * (Dx * Rd * t) ** 0.5))
            + np.exp(x * (U + P) / (2 * Dx))
            * erfc((Rd * x + P * t) / (2 * (Dx * Rd * t) ** 0.5))
        )
    )


def test():
    nx = int(Lx / dx)
    x = np.linspace(0, Lx, nx)
    # xx = np.meshgrid(x)
    for k in range(nt):
        t = (k + 1) * dt
        c = exact_sol(x, t)
        if plt and ((k + 1) % 10 == 0):
            h = plt.plot(x, c)
            plt.title("t=" + str(t))
            plt.show()


if __name__ == "__main__":
    Lx = 100.0
    dx = 2
    tfinal = 100
    dt = 0.2
    Phi = 0.2
    b = 5e03
    rhoCp = 1000 * 2000
    Cm = 60  # température pour x=0
    gamma = 0
    mu = 0
    rhow = 1000
    U = 1  # /(Phi+rhoCp*(1-Phi)/(rhow*b))*b
    Dx = 1
    nt = math.floor(tfinal / dt)
    Rd = rhow * Phi * b + rhoCp * (1 - Phi)
    P = (U**2 + 4 * Dx * Rd * (mu + gamma)) ** (1 / 2)
    print("avant test")
    test()
