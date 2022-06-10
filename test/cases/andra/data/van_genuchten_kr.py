import numpy as np
import numba


# The relative permeability of the liquid phase is expressed by integrating
# the Mualem prediction model in the van Genuchten capillarity model
def vanGenuchten_kr(Pr, Slr, Sgr, n, Sl_reg=0.99):
    m = 1.0 - 1.0 / n
    scale = 1.0 / (1.0 - Slr - Sgr)
    Slb_reg = (Sl_reg - Slr) * scale
    krl_reg = (1.0 - (1.0 - Slb_reg ** (1.0 / m)) ** m) ** 2 * np.sqrt(Slb_reg)
    reg_slope = (1.0 - krl_reg) / (1.0 - Sl_reg)
    reg_origin = (krl_reg - Sl_reg) / (1.0 - Sl_reg)

    @numba.vectorize(["float64(float64)"])
    def fkrg(Sg):
        Slb = (1.0 - Sg - Slr) * scale
        if Slb <= 0.0:  # Sg >= 1-Slr
            return 1.0
        if Slb < 1.0:
            return (1.0 - Slb ** (1.0 / m)) ** (2.0 * m) * np.sqrt(1.0 - Slb)
        return 0.0  # Sg <= to Sgr

    @numba.vectorize(["float64(float64)"])
    def dfkrgdS(Sg):
        Slb = (1.0 - Sg - Slr) * scale
        if Slb <= 0.0:  # Sg >= 1-Slr
            return 0.0
        if Slb < 1.0:
            ss = 1.0 - Slb ** (1.0 / m)
            return scale * (
                2 * Slb ** (1.0 / m - 1.0) * ss ** (2.0 * m - 1.0) * np.sqrt(1.0 - Slb)
                + ss ** (2.0 * m) / 2.0 / np.sqrt(1.0 - Slb)
            )
        return 0.0

    @numba.vectorize(["float64(float64)"])
    def fkrl(Sl):
        Slb = (Sl - Slr) * scale
        if Slb <= 0.0:  # Sl <= Slr
            return 0.0
        if Slb < Slb_reg:  # Sl far enough from regularization
            return (1.0 - (1.0 - Slb ** (1.0 / m)) ** m) ** 2 * np.sqrt(Slb)
        else:  # Sl < 1 - Sgr
            # regularization : line such that f(Sl_reg) = krl_reg
            return reg_slope * Sl + reg_origin

    @numba.vectorize(["float64(float64)"])
    def dfkrldS(Sl):
        Slb = (Sl - Slr) * scale
        if Slb <= 0.0:  # Sl <= Slr
            return 0.0
        if Slb < Slb_reg:  # Sl far enough from regularization
            ss = 1.0 - (1.0 - Slb ** (1.0 / m)) ** m
            ds = Slb ** (1.0 / m - 1.0) * (1.0 - Slb ** (1.0 / m)) ** (m - 1.0)

            return scale * (2.0 * ds * ss * np.sqrt(Slb) + ss**2 / 2.0 / np.sqrt(Slb))
        else:  # Sl < 1 - Sgr
            # regularization
            return reg_slope

    return fkrg, fkrl, dfkrgdS, dfkrldS


kr_laws = {
    1: vanGenuchten_kr(15.0e6, 0.4, 0, 1.49),
    2: vanGenuchten_kr(2.0e6, 0.01, 0, 1.54),
}


def kr_functions(states, rocktypes, kr, dkrdS):
    nb_phases = kr.shape[1]
    S = states.S
    Sg = S[:, 0]
    Sl = S[:, 1]
    assert nb_phases == 2, "Should not be called if np==1"
    kr.fill(0)
    dkrdS.fill(0)
    # Van Genuchten
    for rt, (krg, krl, dkrgdS, dkrldS) in kr_laws.items():
        mask = rocktypes == rt
        kr[mask, 0] = krg(Sg[mask])
        dkrdS[mask, 0, 0] = dkrgdS(Sg[mask])
        kr[mask, 1] = krl(Sl[mask])
        dkrdS[mask, 1, 1] = dkrldS(Sl[mask])


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 1, n)
    krg = np.hstack([kr_laws[1][0](Sg), kr_laws[2][0](Sg)])
    Sg = np.hstack([Sg, Sg])
    plt.plot(Sg, krg, "o")
    plt.savefig("krg.png")

    from scipy.optimize import approx_fprime

    comp = lambda sat, l, eps=1e-10: l[2](sat) - approx_fprime(sat, l[0], eps)
    check = [comp(s, kr_laws[1]) for s in np.linspace(0, 0.6, 100)]
