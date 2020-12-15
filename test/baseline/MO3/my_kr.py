import numpy as np


def fkrl(S):
    result = np.zeros_like(S)
    S0 = 0.3
    dSref = 1.0 - S0
    mask = S > S0
    result[mask] = ((S[mask] - S0) / dSref) ** 2
    return result


def dfkrldS(S):
    result = np.zeros_like(S)
    S0 = 0.3
    dSref = 1.0 - S0
    mask = S > S0
    result[mask] = (2.0 / dSref) * ((S[mask] - S0) / dSref)
    return result


def fkrg(S):
    result = np.ones_like(S)
    S0 = 0.7
    mask = S < S0
    result[mask] = (S[mask] / S0) ** 2
    return result


def dfkrgdS(s):
    result = np.zeros_like(s)
    S0 = 0.7
    mask = S < S0
    result[mask] = (2 / S0) * (S[mask] / S0)
    return result


def kr_functions(states, rocktypes, kr, dkrdS):
    nb_phases = kr.shape[1]
    S = states.S
    Sg = S[:, 0]
    Sl = S[:, 1]
    assert nb_phases == 2, "Should not be called if np==1"
    kr[:, 0] = fkrg(Sg)
    dkrdS[:, 0] = dfkrgdS(Sg)
    kr[:, 1] = fkrldS(Sl)
    dkrdS[:, 1] = dfkrldS(Sl)


# simulation.set_kr_functions(kr_functions)

# import matplotlib.pylab as plt

# S = np.linspace(0, 1)
# plt.plot(S, fkrl(S), label="krl")
# plt.plot(S, fkrg(1-S), label="krg")
# plt.xlabel("Sl")
# plt.legend()

# plt.show()

if __name__ == "__main__":
    S = np.linspace(0, 1)
    ds = 1e-7

    def check(f, df):
        dfapprox = (f(S + ds) - f(S - ds)) / (2 * ds)
        return np.allclose(df(S), dfapprox)

    print(check(fkrl, dfkrldS))
    print(check(fkrg, dfkrgdS))
