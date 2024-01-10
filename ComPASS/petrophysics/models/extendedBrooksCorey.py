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


def extendedBrooksCorey(Slimm, Sgimm, Pe, c2, Pcmax):
    scale_l = 1.0 - Slimm
    c2_pow = -1.0 / c2
    erf_scale = math.sqrt(math.pi) / 2.0 / Pcmax
    ddSerf_scale = Pcmax * 2.0 / math.sqrt(math.pi)

    @numba.vectorize(["float64(float64)"])
    def Pc(Sg):
        norm_s = (1.0 - Sg - Slimm) / scale_l
        if norm_s < 0.0:
            return Pcmax
        Pc_tilde = Pe * norm_s**c2_pow
        return Pcmax * math.erf(Pc_tilde * erf_scale)

    @numba.vectorize(["float64(float64)"])
    def dPcdS(Sg):
        norm_s = (1.0 - Sg - Slimm) / scale_l
        if norm_s < 0.0:
            return 0.0
        Pc_tilde = Pe * norm_s**c2_pow
        ddSPctilde = -c2_pow * Pe / scale_l * norm_s ** (c2_pow - 1)
        return (
            ddSerf_scale
            * math.exp(-((Pc_tilde * erf_scale) ** 2))
            * ddSPctilde
            * erf_scale
        )

    return Pc, dPcdS


# dictionary where the key contains the rocktype and the value
# contains two functions : the capillary pressure Pc(Sg)
# and the derivative dPcdS(Sg)
laws = {
    0: no_pc(),
    # rocktype=1: (Slimm=0.32, Sgimm=0.1, Pe=1500, c2=2, Pcmax=9.5e4)
    1: extendedBrooksCorey(0.32, 0.1, 1500, 2, 9.5e4),
    2: extendedBrooksCorey(0.14, 0.1, 300, 2, 9.5e4),
    3: extendedBrooksCorey(0.12, 0.1, 100, 2, 9.5e4),
    4: extendedBrooksCorey(0.12, 0.1, 25, 2, 9.5e4),
    5: extendedBrooksCorey(0.12, 0.1, 10, 2, 9.5e4),
    6: extendedBrooksCorey(0.10, 0.1, 1, 2, 9.5e4),
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


def set_extendedBrooksCorey_capillary_pressure(simulation):
    # FIXME: cf. holder definition above
    global holder

    holder = phase_pressure(laws)
    simulation.set_phase_pressure_functions(holder)


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 1, n)
    pc1 = laws[1][0](Sg)
    pc3 = laws[3][0](Sg)
    pc4 = laws[4][0](Sg)
    pc5 = laws[5][0](Sg)
    plt.plot(1.0 - Sg, pc1, ".")
    plt.plot(1.0 - Sg, pc3, ".")
    plt.plot(1.0 - Sg, pc4, ".")
    plt.plot(1.0 - Sg, pc5, ".")
    # plt.ylim(1e3, 1.2e5)
    plt.yscale("log")
    plt.xlabel("Sl")
    plt.ylabel("extended Brooks Corey Pc (rt = 1, 3, 4, 5)")
    plt.savefig("Pc.png")
