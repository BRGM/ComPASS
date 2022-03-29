import numpy as np

# Cf. Special Panel on Geothermal Model Intercomparison Study
# Sixth Stanford Annual Workshop on Geothermal Reservoir Engineering
# December 16-18, 1980
# Problem 4: Expanding 2 Phase System with Drainage (same as problem 2)


def Sstar(Sl):
    result = np.zeros_like(Sl)
    mask = Sl > 0.35
    result[mask] = (Sl[mask] - 0.35) / 0.65
    return result


def dSstardSl(Sl):
    result = np.zeros_like(Sl)
    mask = Sl > 0.35
    result[mask] = 1 / 0.65
    return result


def fkrl(Sl):
    return Sstar(Sl) ** 4


def dfkrldS(Sl):
    return 4 * dSstardSl(Sl) * (Sstar(Sl) ** 3)


def fkrg(Sl):
    S = Sstar(Sl)
    return ((1 - S) ** 2) * (1 - S**2)


def dfkrgdS(Sl):
    S = Sstar(Sl)
    dS = dSstardSl(Sl)
    return -2 * dS * (1 - S) * ((1 - S**2) + (1 - S) * S)


def kr_functions(states, rocktypes, kr, dkrdS):
    assert kr.shape[1] == 2, "Should not be called if np!=2"
    dkrdS[...] = 0
    Sl = states.S[:, 1]
    kr[:, 0] = fkrg(Sl)
    dkrdS[:, 0, 0] = -dfkrgdS(Sl)  # derivative is w.r.t. Sg
    kr[:, 1] = fkrl(Sl)
    dkrdS[:, 1, 1] = dfkrldS(Sl)


if __name__ == "__main__":

    S = np.linspace(0, 1, 100)

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("No figure is drawn as matpplotlib is not installed.")
    else:
        plt.clf()
        plt.plot(S, fkrl(S), label="krl")
        plt.plot(S, fkrg(S), label="krg")
        plt.xlabel("Sl")
        plt.legend()
        plt.savefig("kr.png")

    ds = 1e-7

    def check(f, df):
        dfapprox = (f(S + ds) - f(S - ds)) / (2 * ds)
        return np.allclose(df(S), dfapprox)

    print(check(Sstar, dSstardSl))
    print(check(fkrl, dfkrldS))
    print(check(fkrg, dfkrgdS))
