import numpy as np

from ComPASS.properties.water2ph import liquid, gas, diphasic

Tsat = diphasic.Tsat
dTsatdp = diphasic.dTsatdp
mul = lambda p: liquid.mu(p, Tsat(p))
rhol = lambda p: liquid.rho(p, Tsat(p))
hl = lambda p: liquid.h(p, Tsat(p))
mug = lambda p: gas.mu(p, Tsat(p))
rhog = lambda p: gas.rho(p, Tsat(p))
hg = lambda p: gas.h(p, Tsat(p))


# from iapws import IAPWS97
# from iapws.iapws97 import _TSat_P

# dTsatdp = lambda p: 0.5*(_TSat_P(1e-6*(p+1))-_TSat_P(1e-6*(p-1)))
# mul = lambda p: IAPWS97(P=1e-6*p, x=0).mu
# rhol = lambda p: IAPWS97(P=1e-6*p, x=0).rho
# hl = lambda p: 1e3 * IAPWS97(P=1e-6*p, x=0).h
# mug = lambda p: IAPWS97(P=1e-6*p, x=1).mu
# rhog = lambda p: IAPWS97(P=1e-6*p, x=1).rho
# hg = lambda p: 1e3 * IAPWS97(P=1e-6*p, x=1).h


krl = krg = lambda S: S
# krl = krg = lambda S: (S-0.3)/0.7 if S>0.3 else 0

Pc = lambda S: S / (1.0001 - S)
dPcdS = lambda S: 1.0001 / (1.0001 - S) ** 2


def build(
    k=1e-16,
    K=2,
    g=9.81,
    phi=2.0,
    # dSdz=0,
):
    def dpdz_mass(p, Sg):
        rholp = rhol(p)
        rhogp = rhog(p)
        Ml = rholp * k * krl(1 - Sg) / mul(p)
        Mg = rhogp * k * krg(Sg) / mug(p)
        # return (Ml*dSdz*dPcdS(Sg)-(Ml*rholp+Mg*rhogp)*g)/(Ml+Mg)
        return -(Ml * rholp + Mg * rhogp) * g / (Ml + Mg)

    def dpdz_energy(p, Sg):
        rholp = rhol(p)
        rhogp = rhog(p)
        Ml = rholp * k * krl(1 - Sg) / mul(p)
        Mg = rhogp * k * krg(Sg) / mug(p)
        hlp = hl(p)
        hgp = hg(p)
        # return (phi+Ml*hlp*dSdz*dPcdS(Sg)+(Ml*rholp*hlp+Mg*rhogp*hgp)*g)/(Ml*hlp+Mg*hgp-K*dTsatdp(p))
        return -(phi + (Ml * rholp * hlp + Mg * rhogp * hgp) * g) / (
            Ml * hlp + Mg * hgp + K * dTsatdp(p)
        )

    return dpdz_mass, dpdz_energy


import matplotlib.pyplot as plt

# Sg = np.linspace(0.9, 1)
# plt.clf()
# for p in [1e6, 2e6, 7e6, 8e6, 1e7]:
#     k=5e-17
#     K=2.
#     phi = K*(Tsat(p)-293.15)/2000.
#     print("Bottom heat flux:", phi, f"Tsat({p*1.e-6:.1f}MPa)={Tsat(p)-273.15}")
#     dpm, dpe = build(k=k, K=K, phi=phi)
#     # plt.plot(Sg, [dpdz_mass(p, Sgk) for Sgk in Sg])
#     # plt.plot(Sg, [dpdz_energy(p, Sgk) for Sgk in Sg])
#     plt.plot(Sg, [dpm(p, Sgk)-dpe(p, Sgk) for Sgk in Sg], label=f"p={p*1.e-6:.1f}MPa")
# # plt.ylim(-10, 10)
# plt.legend()
# plt.plot((Sg.min(),Sg.max()),(0,0), '-k')


from scipy.optimize import root_scalar, minimize_scalar

plt.clf()
p = 8e6
K = 2.0
phi = K * (Tsat(p) - 293.15) / 2000.0
print("Bottom heat flux:", phi, f"Tsat({p*1.e-6:.1f}MPa)={Tsat(p)-273.15}")

k = np.linspace(-18, -10)
mtol = 1e-6
xtol = 1e-8
sol = []
for ki in k:
    dpm, dpe = build(k=10**ki, K=K, phi=phi)
    f = lambda S: abs(dpm(p, S) - dpe(p, S))
    res = minimize_scalar(f, method="brent", bracket=(0.0, 1.0), tol=mtol)
    message = f"k={10**ki:.2g} : S={res.x:.3g} |f(S)|={f(res.x):.5g}"
    x = res.x
    f = lambda S: dpm(p, S) - dpe(p, S)
    if f(0.5 * (1 + x)) * f(1.0) < 0:
        res = root_scalar(f, method="toms748", bracket=(0.5 * (1 + x), 1.0), xtol=xtol)
        message += f" root at S={res.root:.8e}"
        sol.append((ki, res.root))
    print(message)

sol = np.array(sol)
plt.plot(sol[:, 0], sol[:, 1])
plt.yscale("log")

plt.savefig("heat-pipe.png")
