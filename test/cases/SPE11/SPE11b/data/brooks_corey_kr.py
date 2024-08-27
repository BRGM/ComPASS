import numpy as np

# Brooks Corey type rel perm
def fkr(S, Simm, alpha):
    Sred = (S - Simm) / (1.0 - Simm)
    mask = Sred > 0.0
    Sn = np.zeros_like(S)
    Sn[mask] = Sred[mask]
    # derivative of Sn**alpha
    dfkrdS = np.zeros_like(S)
    if np.isscalar(Simm):
        dfkrdS[mask] = alpha / (1.0 - Simm) * Sn[mask] ** (alpha - 1)
    else:
        dfkrdS[mask] = alpha / (1.0 - Simm[mask]) * Sn[mask] ** (alpha - 1)
    return Sn**alpha, dfkrdS


def make_kr_functions(Sgimm, Slimm, alpha):
    def kr_functions(states, rocktypes, kr, dkrdS):
        nb_phases = kr.shape[1]
        S = states.S
        Sg = S[:, 0]
        Sl = S[:, 1]
        assert nb_phases == 2, "Should not be called if np==1"
        # kr.fill(0)
        dkrdS.fill(0)
        kr[:, 0], dkrdS[:, 0, 0] = fkr(Sg, Sgimm, alpha)
        kr[:, 1], dkrdS[:, 1, 1] = fkr(Sl, Slimm[rocktypes], alpha)

    return kr_functions


# SPE11a: Sgimm = 0.1, Slimm depends on the rocktype, alpha = 2
Slimm = np.array([0, 0.32, 0.14, 0.12, 0.12, 0.12, 0.10])
kr_functions_SPE11a = make_kr_functions(0.1, Slimm, 2)
# SPE11b:
# Sgimm = 0.1, Slimm depends on the rocktype (same as test a), alpha = 1.5
kr_functions_SPE11b = make_kr_functions(0.1, Slimm, 1.5)


if __name__ == "__main__":
    import matplotlib.pylab as plt

    def check(f, S, Simm, alpha, ds=1e-7):
        fplus, _ = f(S + ds, Simm, alpha)
        fminus, _ = f(S - ds, Simm, alpha)
        _, df = f(S, Simm, alpha)
        dfapprox = (fplus - fminus) / (2 * ds)
        return np.allclose(df, dfapprox)

    n = 1000
    S = np.linspace(0, 1, n)
    krg, dkrg = fkr(1 - S, 0.1, 1.5)
    assert check(fkr, S, 0.1, 1.5)
    krl1, _ = fkr(S, Slimm[1], 1.5)
    krl3, _ = fkr(S, Slimm[3], 1.5)
    plt.plot(S, krg, ".", label="krg")
    plt.plot(S, krl1, ".", label="krl, rt=1")
    plt.plot(S, krl3, ".", label="krl, rt=3")
    plt.xlabel("Sl")
    plt.ylabel("SPE11b Brooks Corey kr")
    plt.legend()
    plt.savefig("kr_brooks_corey.png")
