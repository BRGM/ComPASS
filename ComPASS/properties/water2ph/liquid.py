import numpy as np
from .diphasic import psat, dpsatdT


def h(p, T):
    a = -14.4319e3
    b = 4.70915e3
    c = -4.87534
    d = 1.45008e-2
    T0 = 273.0
    dT = T - T0
    return a + b * dT + c * dT**2 + d * dT**3


def dhdp(p, T):
    return 0


def dhdT(p, T):
    a = -14.4319e3
    b = 4.70915e3
    c = -4.87534
    d = 1.45008e-2
    T0 = 273.0
    dT = T - T0
    return b + 2 * c * dT + 3 * d * dT**2


def mu(p, T):
    ss = (
        0.021482 * (T - 273.0 - 8.435 + np.sqrt(8078.4 + (T - 273.0 - 8.435) ** 2))
        - 1.2
    )
    return 1.0e-3 / ss


def dmudp(p, T):
    return 0


def dmudT(p, T):
    ss = (
        0.021482 * (T - 273.0 - 8.435 + np.sqrt(8078.4 + (T - 273.0 - 8.435) ** 2))
        - 1.2
    )
    dssdT = 0.021482 * (
        1.0 + (T - 273.0 - 8.435) / np.sqrt(8078.4 + (T - 273.0 - 8.435) ** 2)
    )
    return -1.0e-3 * dssdT / ss**2


def rho(p, T):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    ss = rho0 + a * T + b * T**2
    cw = a1 + a2 * p + T * (b1 + b2 * p) + T**2 * (c1 + c2 * p)
    return ss * (1.0 + cw * (p - psat(T)))


def drhodp(p, T):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    ss = rho0 + a * T + b * T**2
    cw = a1 + a2 * p + T * (b1 + b2 * p) + T**2 * (c1 + c2 * p)
    dcwdp = a2 + b2 * T + c2 * T**2
    return ss * (dcwdp * (p - psat(T)) + cw)


def drhodT(p, T):
    rho0 = 780.83795
    a = 1.6269192
    b = -3.0635410e-3
    a1 = 2.4638e-9
    a2 = 1.1343e-17
    b1 = -1.2171e-11
    b2 = 4.8695e-20
    c1 = 1.8452e-14
    c2 = -5.9978e-23
    ss = rho0 + a * T + b * T**2
    dssdT = a + 2 * b * T
    cw = a1 + a2 * p + T * (b1 + b2 * p) + T**2 * (c1 + c2 * p)
    dcwdT = b1 + b2 * p + 2 * T * (c1 + c2 * p)
    return dssdT * (1.0 + cw * (p - psat(T))) + ss * (
        dcwdT * (p - psat(T)) - cw * dpsatdT(T)
    )
