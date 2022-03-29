import numpy as np


def fkr(S):
    return S**2


def dfkrdS(S):
    return 2 * S


def kr_functions(states, rocktypes, kr, dkrdS):
    nb_phases = kr.shape[1]
    S = states.S
    Sg = S[:, 0]
    Sl = S[:, 1]
    assert nb_phases == 2, "Should not be called if np==1"
    kr.fill(0)
    dkrdS.fill(0)
    kr[:, 0] = fkr(Sg)
    dkrdS[:, 0, 0] = dfkrdS(Sg)
    kr[:, 1] = fkr(Sl)
    dkrdS[:, 1, 1] = dfkrdS(Sl)


# import matplotlib.pylab as plt

# S = np.linspace(0, 1)
# plt.plot(S, fkrl(S), label="krl")
# plt.plot(S, fkrg(1-S), label="krg")
# plt.xlabel("Sl")
# plt.legend()

# plt.show()

# if __name__ == "__main__":
#     S = np.linspace(0, 1)
#     ds = 1e-7

#     def check(f, df):
#         dfapprox = (f(S + ds) - f(S - ds)) / (2 * ds)
#         return np.allclose(df(S), dfapprox)

#     print(check(fkrl, dfkrldS))
#     print(check(fkrg, dfkrgdS))
