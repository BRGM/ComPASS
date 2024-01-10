import numpy as np
import numba

# FIXME: is this a pybind11 bug?
# This is used as global variable
# If not use the reference counter of the closure
# returned by _convert_pc_to_phase_pressure_function
# goes mad at the end of the program execution
holder = None


def vanGenuchten(Pr, Slr, Sgr, n, Slb_reg=0.99):
    minus_one_over_m = n / (1.0 - n)  # usualy we define m = 1 - 1/n and slb ** (-1/m)
    # the law is regularized to allow for Sg = Sgr
    Pc1 = Pr * (Slb_reg**minus_one_over_m - 1.0) ** (1.0 / n)
    alpha = Pc1 / (Slb_reg - 1.0)
    eps = 1.0e-7  # eps is not so small because the Pc is high

    @numba.vectorize(["float64(float64)"])
    def Pc(Sg):
        Slb = (1.0 - Sg - Slr) / (1.0 - Sgr - Slr)

        assert Slb > eps  # stay at a certain distance of the asymptote

        if Slb < Slb_reg:  # Sg far enough from Sgr
            return Pr * (Slb**minus_one_over_m - 1.0) ** (1.0 / n)
        else:  # Sg close or equal to Sgr
            # the law is regularized to allow for Sg = Sgr
            return alpha * (Slb - Slb_reg) + Pc1
        # if Slb > 1.0 :  # maybe useful ?
        #     return 1.e15

    @numba.vectorize(["float64(float64)"])
    def dPcdS(Sg):
        Slb = (1.0 - Sg - Slr) / (1.0 - Sgr - Slr)
        dSlbdS = -1.0 / (1.0 - Sgr - Slr)
        assert Slb > eps  # stay at a certain distance of the asymptote

        if Slb < Slb_reg:
            return (
                Pr
                * minus_one_over_m
                / n
                * dSlbdS
                * Slb ** (minus_one_over_m - 1.0)
                * (Slb**minus_one_over_m - 1.0) ** (1.0 / n - 1.0)
            )
        else:
            return alpha * dSlbdS
        # if Slb > 1.0:
        #     return 1.e15

    return Pc, dPcdS


# dictionary where the key contains the rocktype and the value
# contains two functions : the capillary pressure Pc(Sg)
# and the derivative dPcdS(Sg)
laws = {
    # rocktype=1, Pr=15e6, Slr=0.4, Sgr=0, n=1.49, Slb_reg=0.99
    1: vanGenuchten(15.0e6, 0.4, 0, 1.49),
    # rocktype=2, Pr=2e6, Slr=0.01, Sgr=0, n=1.54, Slb_reg=0.99
    2: vanGenuchten(2.0e6, 0.01, 0, 1.54),
    4: vanGenuchten(2.0e6, 0, 0, 1.5),  # Cigeo backfill in access drifts
    5: vanGenuchten(1.6e7, 0, 0, 1.6),  # Cigeo Bentonite
    6: vanGenuchten(1.5e6, 0, 0, 1.5),  # Cigeo EDZ
    7: vanGenuchten(1.5e7, 0, 0, 1.5),  # Cigeo argilite (natural medium)
    8: vanGenuchten(
        1.5e7, 0, 0, 1.5
    ),  # Cigeo acquifer (on the top of the natural medium)
}


def phase_pressure(laws):

    # @numba.njit()
    # compute P_alpha and its partial derivatives
    # knowing the states X
    # This function defines that Pref = Palpha[0] which is the gas pressure
    def phase_pressure_function(X, rocktypes, dpadS):
        for rt, (Pc, dPcdS) in laws.items():
            mask = rocktypes == rt
            X.pa[mask, 1] = -Pc(X.S[mask, 0])
            # derivative of p_alpha wrt S_alpha, only non zero when alpha = liquid
            # d(pl)/d(Sl) = - d(pc)/d(Sl) = d(pc)/d(Sg)
            dpadS[mask, 1] = dPcdS(X.S[mask, 0])
        X.pa[:, 1] += X.p  # Pl = Pg - Pc, X.p being the reference pressure ie Pg
        X.pa[:, 0] = X.p
        # dpadS[:, 0] = 0.0

    return phase_pressure_function


def set_vanGenuchten_capillary_pressure(simulation):
    # FIXME: cf. holder definition above
    global holder

    holder = phase_pressure(laws)
    simulation.set_phase_pressure_functions(holder)


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 0.599, n)
    pc = laws[1][0](Sg)
    # pc = np.hstack([laws[1][0](Sg), laws[2][0](Sg)])
    # Sg = np.hstack([Sg, Sg])
    plt.plot(Sg, pc, ".")
    plt.ylim(1e-1, 1e15)
    plt.yscale("log")
    plt.xlabel("Sg")
    plt.ylabel("van Genuchten Pc (rt = 1)")
    plt.savefig("Pc.png")
