# This program computes the analutical solution from XXX
# which is available through the function solution
# Then it dumps to the file SO1-analytical.csv
# the analytical solution which is sampled on space according to
# length Lx and discretization dx
# at times given by vector t
# if matplotlib is available it will also draw a picture

import numpy as np
from scipy.special import erfc


def build_solution(
    porosity=0.2,
    permeability=1,  # D
    rhow=1000,  # kg/m3
    Cw=4200.0,  # J/K/kg
    Kw=0.6,  # W/m °K)
    rhor=2600,  # kg/m3
    Cr=800,  # J/K/kg
    Kr=2.0,  # W/m °K)
    Tres=5,  # °C
    Tinj=33,  # °C
    Qinj=1.0,
    Syz=1.0,
):
    # porosity weighted barycenter / arithmetic mean
    def poro_ar_mean(xpor, xrock):
        return xpor * porosity + xrock * (1.0 - porosity)

    # auxiliary variables
    rhoCw = rhow * Cw
    rhoCr = rhor * Cr
    rhoCeq = poro_ar_mean(rhoCw, rhoCr)  # J/K/m3
    Keq = poro_ar_mean(Kw, Kr)  # W/(m.K°)
    udarcy = Qinj / Syz
    upore = udarcy / porosity
    rhowCwupore = rhow * Cw * udarcy / rhoCeq
    KeqovrhoCeq = Keq / (porosity * rhoCeq)
    # function that is returned
    def solution(t, x):
        x = np.asarray(x)
        where = (rhowCwupore / KeqovrhoCeq) * x < 700.0
        result = np.zeros(x.shape, dtype=np.double)
        xscale = 2 * np.sqrt(KeqovrhoCeq * t)
        xw = x[where]
        DT = Tinj - Tres
        xt = rhowCwupore * t
        result[where] = Tres + DT / 2 * (
            erfc((xw - xt) / xscale)
            + np.exp(rhowCwupore * xw / KeqovrhoCeq) * erfc((xw + xt) / xscale)
        )
        where = np.logical_not(where)
        xw = x[where]
        result[where] = Tres + DT / 2 * erfc((xw - xt) / xscale)
        return result

    return solution


Lx, Ly, Lz = 870, 10.0, 10.0
dx = 10.0

# below you can speciy any parameter that is accepted by build_solution
solution = build_solution(Qinj=3e-5, Syz=Ly * Lz,)

day = 86400.0  # seconds
x = np.arange(0, Lx + 0.5 * dx, dx)
t = day * np.hstack([0.1, np.arange(20, 1020, 30), 1e8 / day])
table = np.array([solution(ti, x) for ti in t])

# dump table to csv
csv_table = np.hstack(
    [
        np.reshape(t / day, (-1, 1)),  # time in days
        np.reshape(t, (-1, 1)),  # time in seconds
        table,
    ]
)
np.savetxt(
    "SO1-analytical.csv",
    csv_table,
    delimiter=";",
    header="days;seconds;" + ";".join("%.1f" % xi for xi in x),
)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("No figure is drawn as matpplotlib is not installed.")
    raise
else:
    plt.clf()
    # Draw last time
    plt.plot(x, table[-1])
    plt.xlabel("distance (m)")
    plt.ylabel("temperature (deg C)")
    plt.grid(True)
    plt.savefig("SO1-solution.png")
