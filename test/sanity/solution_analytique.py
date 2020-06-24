import numpy as np
import matplotlib.pyplot as plt
from scipy import special
import time
import math


def f(tau, x, y):
    return (
        tau ** (-3 / 2)
        * (
            special.erf((a + y) / math.sqrt(4 * Dt * tau))
            + special.erf((a - y) / math.sqrt(4 * Dt * tau))
        )
        * np.exp(-((x - tau * v) ** 2) / (4 * Dl * tau))
    )  # y pas encore défini


def exact_sol(x, y, t):
    n = 100
    s = 0
    for i in range(n):
        tau = t * (np.cos((2 * (i + 1) - 1) * np.pi / (4 * n))) ** 2
        s = s + x * t * c0 / (64 * np.pi * Dl) ** (1 / 2) * np.pi / n * np.sin(
            (2 * (i + 1) - 1) * np.pi / (2 * n)
        ) * f(tau, x, y)
    return s


def test():
    nx = math.floor(Lx / dx)
    ny = math.floor(Ly / dy)
    x = np.linspace(0, Lx, nx)
    y = np.linspace(-Ly / 2, Ly / 2, ny)
    xx, yy = np.meshgrid(x, y)

    for k in range(nt):
        t = (k + 1) * dt
        c2 = exact_sol(xx, yy, t)
        c1 = exact_sol(x, 0.0, t)
        if (k + 1) % 10 == 0:
            fig = plt.figure(1)
            plt.subplot(211)
            cs = plt.contourf(x, y, c2)
            fig.colorbar(cs)
            plt.subplot(212)
            plt.plot(x, c1)
            plt.title("t=" + str(t))
            plt.draw()
            plt.pause(0.1)
            plt.clf()


if __name__ == "__main__":
    Dl = 0.1
    Dt = 0.01
    c0 = 10
    a = 8
    Pe = 10
    v = 1
    Lx, Ly = 100, 40
    dx, dy = 2, 2  # dx=xy=2
    tfinal = 100
    dt = 0.2
    nt = math.floor(tfinal / dt)
    test()
