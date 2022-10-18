import numpy as np


def psat(T):
    return 1.0e-3 * (T - 273.0) ** 4.0


def dpsatdT(T):
    return 4.0e-3 * (T - 273.0) ** 3.0


def Tsat(p):
    return 273.0 + 100.0 * (1.0e-5 * p) ** 0.25


def dTsatdp(p):
    return 25.0e-5 * (1.0e-5 * p) ** -0.75
