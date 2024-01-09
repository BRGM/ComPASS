import numpy as np
import numba

# FIXME: is this a pybind11 bug?
# This is used as global variable
# If not use the reference counter of the closure
# returned by _convert_pc_to_phase_pressure_function
# goes mad at the end of the program execution
holder = None


def BrooksCorey(Slimm, Sgimm, Pe, c2):
    scale_g = 1.0 - Sgimm
    scale_l = 1.0 - Slimm
    c2_pow = -1.0 / c2

    @numba.vectorize(["float64(float64)"])
    def Pc(Sg):
        norm_s = max((1.0 - Sg - Slimm) / scale_l, 0)
        return Pe * norm_s**c2_pow

    @numba.vectorize(["float64(float64)"])
    def dPcdS(Sg):
        norm_s = (1.0 - Sg - Slimm) / scale_l
        if norm_s < 0:
            return 0.0
        return -c2_pow * Pe / scale_l * norm_s ** (c2_pow - 1)

    return Pc, dPcdS


# dictionary where the key contains the rocktype and the value
# contains two functions : the capillary pressure Pc(Sg)
# and the derivative dPcdS(Sg)
laws = {
    # rocktype=1, Slimm=0, Sgimm=0.1, Pe, c2=2
    1: BrooksCorey(0.0, 0.1, 50, 2),
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


def set_BrooksCorey_capillary_pressure(simulation):
    # FIXME: cf. holder definition above
    global holder

    holder = phase_pressure(laws)
    simulation.set_phase_pressure_functions(holder)


if __name__ == "__main__":
    import matplotlib.pylab as plt

    n = 1000
    Sg = np.linspace(0, 0.999, n)
    pc = laws[1][0](Sg)
    # pc = np.hstack([laws[1][0](Sg), laws[2][0](Sg)])
    # Sg = np.hstack([Sg, Sg])
    plt.plot(Sg, pc, ".")
    # plt.ylim(1e-1, 1e15)
    plt.yscale("log")
    plt.xlabel("Sg")
    plt.ylabel("Brooks Corey Pc (rt = 1)")
    plt.savefig("Pc.png")
