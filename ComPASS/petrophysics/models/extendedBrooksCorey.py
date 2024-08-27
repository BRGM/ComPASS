import numpy as np
import numba
import math

# FIXME: is this a pybind11 bug?
# This is used as global variable
# If not use the reference counter of the closure
# returned by _convert_pc_to_phase_pressure_function
# goes mad at the end of the program execution
holder = None


def no_pc():
    @numba.vectorize(["float64(float64)"])
    def Pc(Sg):
        return 0.0

    @numba.vectorize(["float64(float64)"])
    def dPcdS(Sg):
        return 0.0

    return Pc, dPcdS


# S_tilde is a value chosen to smooth the asymptote
# introduce a line to extend the curve
def extendedBrooksCorey(Slimm, Sgimm, Pe, c2, Pcmax, S_tilde):
    # S_tilde must be greater than Slimm
    assert S_tilde > Slimm
    scale_l = 1.0 - Slimm
    c2_pow = -1.0 / c2
    erf_scale = math.sqrt(math.pi) / 2.0 / Pcmax
    ddSerf_scale = Pcmax * 2.0 / math.sqrt(math.pi)

    # determine line parameters such that the derivative is continius
    def compute_a_b():
        norm_s_tilde = (S_tilde - Slimm) / scale_l
        Pc_tilde_s_tilde = Pe * norm_s_tilde**c2_pow
        Pc_s_tilde = Pcmax * math.erf(Pc_tilde_s_tilde * erf_scale)
        ddSPctilde_s_tilde = c2_pow * Pe / scale_l * norm_s_tilde ** (c2_pow - 1)
        a = (
            ddSerf_scale
            * math.exp(-((Pc_tilde_s_tilde * erf_scale) ** 2))
            * ddSPctilde_s_tilde
            * erf_scale
        )
        b = Pc_s_tilde
        return a, b

    # lines paramters
    a, b = compute_a_b()

    @numba.vectorize(["float64(float64)"])
    def Pc(Sg):
        if (1.0 - Sg) >= S_tilde:
            norm_s = (1.0 - Sg - Slimm) / scale_l
            Pc_tilde = Pe * norm_s**c2_pow
            return Pcmax * math.erf(Pc_tilde * erf_scale)
        else:
            return a * ((1.0 - Sg) - S_tilde) + b

    @numba.vectorize(["float64(float64)"])
    def dPcdS(Sg):
        if (1.0 - Sg) >= S_tilde:
            norm_s = (1.0 - Sg - Slimm) / scale_l
            Pc_tilde = Pe * norm_s**c2_pow
            ddSPctilde = -c2_pow * Pe / scale_l * norm_s ** (c2_pow - 1)
            return (
                ddSerf_scale
                * np.exp(-((Pc_tilde * erf_scale) ** 2))
                * ddSPctilde
                * erf_scale
            )
        else:
            # todo:
            # error: the derivative wrt Sg should be -a, but the convergence
            # is not obtained with -a whereas it is quite easy with a.
            # The derivative value should have an impact on the convergence
            # behaviour but not on the value after convergence.
            # todo: write a 1D test
            return a

    return Pc, dPcdS


# dictionary where the key contains the rocktype and the value
# contains two functions : the capillary pressure Pc(Sg)
# and the derivative dPcdS(Sg)
SPE11a_laws = {
    0: no_pc(),
    # rocktype=1: (Slimm=0.32, Sgimm=0.1, Pe=1500, c2=2, Pcmax=9.5e4)
    1: extendedBrooksCorey(0.32, 0.1, 1500, 2, 9.5e4, 0.325),
    2: extendedBrooksCorey(0.14, 0.1, 300, 2, 9.5e4, 0.145),
    3: extendedBrooksCorey(0.12, 0.1, 100, 2, 9.5e4, 0.125),
    4: extendedBrooksCorey(0.12, 0.1, 25, 2, 9.5e4, 0.125),
    5: extendedBrooksCorey(0.12, 0.1, 10, 2, 9.5e4, 0.125),
    6: extendedBrooksCorey(0.10, 0.1, 1, 2, 9.5e4, 0.105),
}

# SPE11b_laws = {
#     0: no_pc(),
#     # rocktype=1: (Slimm=0.32, Sgimm=0.1, Pe=1500, c2=1.5, Pcmax=3.e7, S_tilde=0.38)
#     1: extendedBrooksCorey(0.32, 0.1, 1500, 1.5, 3.0e7, 0.327),
#     2: extendedBrooksCorey(0.14, 0.1, 300, 1.5, 3.0e7, 0.147),
#     3: extendedBrooksCorey(0.12, 0.1, 100, 1.5, 3.0e7, 0.127),
#     4: extendedBrooksCorey(0.12, 0.1, 25, 1.5, 3.0e7, 0.127),
#     5: extendedBrooksCorey(0.12, 0.1, 10, 1.5, 3.0e7, 0.127),
#     6: extendedBrooksCorey(0.10, 0.1, 1, 1.5, 3.0e7, 0.107),
# }
SPE11b_laws = {
    0: no_pc(),
    # rocktype=1: (Slimm=0.32, Sgimm=0.1, Pe=1500, c2=1.5, Pcmax=3.e7, S_tilde=0.38)
    1: extendedBrooksCorey(0.32, 0.1, 1500, 1.5, 3.0e7, 0.35),
    2: extendedBrooksCorey(0.14, 0.1, 300, 1.5, 3.0e7, 0.145),
    3: extendedBrooksCorey(0.12, 0.1, 100, 1.5, 3.0e7, 0.127),
    4: extendedBrooksCorey(0.12, 0.1, 25, 1.5, 3.0e7, 0.127),
    5: extendedBrooksCorey(0.12, 0.1, 10, 1.5, 3.0e7, 0.122),
    6: extendedBrooksCorey(0.10, 0.1, 1, 1.5, 3.0e7, 0.112),
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


def set_extendedBrooksCorey_pc_SPE11a(simulation):
    # FIXME: cf. holder definition above
    global holder

    holder = phase_pressure(SPE11a_laws)
    simulation.set_phase_pressure_functions(holder)


def set_extendedBrooksCorey_pc_SPE11b(simulation):
    # FIXME: cf. holder definition above
    global holder

    holder = phase_pressure(SPE11b_laws)
    simulation.set_phase_pressure_functions(holder)


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 1, n)
    pc1 = SPE11b_laws[1][0](Sg)
    pc2 = SPE11b_laws[2][0](Sg)
    pc3 = SPE11b_laws[3][0](Sg)
    pc4 = SPE11b_laws[4][0](Sg)
    pc5 = SPE11b_laws[5][0](Sg)
    pc6 = SPE11b_laws[6][0](Sg)
    plt.plot(1.0 - Sg, pc1, ".", label="rt=1")
    plt.plot(1.0 - Sg, pc2, ".", label="rt=2")
    plt.plot(1.0 - Sg, pc3, ".", label="rt=3")
    plt.plot(1.0 - Sg, pc4, ".", label="rt=4")
    plt.plot(1.0 - Sg, pc5, ".", label="rt=5")
    plt.plot(1.0 - Sg, pc6, ".", label="rt=6")
    # plt.ylim(1e3, 1.2e5)
    # plt.xlim(0.1, 0.35)
    # plt.xticks(np.arange(0.1, 0.35, 0.02))
    plt.yscale("log")
    plt.xlabel("Sl")
    plt.ylabel("SPE11b extended Brooks Corey Pc")
    plt.legend()
    plt.savefig("Pc_1.png")
