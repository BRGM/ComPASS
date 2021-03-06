import time
from math import *
import numpy as np
import ComPASS
import ComPASS.utils.mpl_backends as mpl_backends

plt = mpl_backends.import_pyplot(False)

try:
    from scipy.special import erf, erfc
except ImportError:
    # this is not really efficient just a workaround when scipy is missing
    # good for government work
    erf = np.vectorize(erf)
    erfc = np.vectorize(erfc)

# Dx = 1
# P = (U**2 + 4*Dx*Rd*(mu+gamma))**(1/2)


def exact_sol_bs(x, y, t):
    exact_sol = (
        epsilon
        * C0
        / (2 * n0)
        * (
            erfc((x - mu * t) / sqrt(4 * Dl * t))
            + np.exp(mu * x / Dl) * erfc((x + mu * t) / sqrt(4 * Dl * t))
        )
    )
    for i in range(1, n):
        l1 = (
            C0
            / (i * np.pi)
            * np.sin(i * np.pi * epsilon / n0)
            * np.cos(i * np.pi * y / n0)
        )
        l2 = np.exp(
            x
            / 2
            * (mu / Dl - sqrt((mu / Dl) ** 2 + (2 * i * np.pi / n0) ** 2 * Dt / Dl))
        )
        l3 = erfc(
            (x - Dl * t * sqrt((mu / Dl) ** 2 + (2 * i * np.pi / n0) ** 2 * Dt / Dl))
            / sqrt(4 * Dl * t)
        )
        l4 = np.exp(
            x
            / 2
            * (mu / Dl + sqrt((mu / Dl) ** 2 + (2 * i * np.pi / n0) ** 2 * Dt / Dl))
        )
        l5 = erfc(
            (x + Dl * t * sqrt((mu / Dl) ** 2 + (2 * i * np.pi / n0) ** 2 * Dt / Dl))
            / sqrt(4 * Dl * t)
        )
        exact_sol = exact_sol + l1 * (l2 * l3 + l4 * l5)
    return exact_sol


def f(tau, x, y):
    return (
        tau ** (-3 / 2)
        * (
            erf((epsilon + y) / sqrt(4 * Dt * tau))
            + erf((epsilon - y) / sqrt(4 * Dt * tau))
        )
        * np.exp(-((x - tau * mu) ** 2) / (4 * Dl * tau))
    )


def exact_sol_ld(x, y, t):
    n = 100
    s = 0
    for i in range(n):
        tau = t * (np.cos((2 * (i + 1) - 1) * np.pi / (4 * n))) ** 2
        s = s + x * t * C0 / (64 * np.pi * Dl) ** (1 / 2) * np.pi / n * np.sin(
            (2 * (i + 1) - 1) * np.pi / (2 * n)
        ) * f(tau, x, y)
    return s


def test():
    nx = 100
    ny = 80
    nz = 1
    x = np.linspace(0, Lx, nx)
    y = np.linspace(-Ly / 2, Ly / 2, ny)
    xx, yy = np.meshgrid(x, y)
    for k in range(nt):
        t = (k + 1) * dt
        print(k)
        cbs = exact_sol_bs(xx, yy, t) + Tright
        cld = exact_sol_ld(xx, yy, t) + Tright
        if plt and ((k + 1) % 10 == 0):
            fig = plt.figure(1)
            plt.subplot(211)
            cs1 = plt.contourf(x, y, cbs)
            fig.colorbar(cs1)
            plt.axis("scaled")
            plt.title("Bruch Street")
            plt.subplot(212)
            cs2 = plt.contourf(x, y, cld)
            fig.colorbar(cs2)
            plt.axis("scaled")
            plt.title("Leij Dane")
            plt.suptitle("t=" + str(t))
            plt.draw()
            plt.pause(1)
            plt.savefig("exact_sol_" + str(t) + ".png")
            plt.clf()


if __name__ == "__main__":
    tfinal = 100  # 2e+5
    dt = 1  # 1e+3
    nt = floor(tfinal / dt)
    Phi = 0.2
    rhow = 1  # 1000
    b = 1  # 4.2e+3  #5e+03
    rhoCp = 1  # 800*2000
    Rd = rhow * Phi * b + rhoCp * (1 - Phi)
    print(Rd)
    # U = 1/(Phi+rhoCp*(1-Phi)/(rhow*b))*b*rhow #*densité molaire à priori constante
    U = 1  # 6.6e-4 #vitesse de darcy
    mu = b * rhow * U / Rd  # vitesse de déplacement de la thermique
    # mu = 1
    print(mu)
    Tleft, Tright = 60, 100
    C0 = Tleft - Tright  # température en x=0
    gamma = 0  # pas sur dans l'expression de P
    Dl = 10  # 8e+6/Rd  #conductivité thermique/Rd
    Dt = 10  # 8e+6/Rd
    print(Dt)
    Lx, Ly, Lz = 100.0, 80.0, 1.0
    n0 = Ly
    epsilon = 8
    n = 100
    test()
