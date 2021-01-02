def h(p, T):
    a = 1990.89e3
    b = 1901.6
    return a + T * b


def dhdp(p, T):
    return 0


def dhdT(p, T):
    b = 1901.6
    return b


def mu(p, T):
    return (0.361 * T - 10.2) * 1.0e-7


def dmudp(p, T):
    return 0


def dmudT(p, T):
    return 3.61e-8


def rho(p, T):
    u = 0.018016
    R = 8.3145
    return p * u / (R * T)


def drhodp(p, T):
    u = 0.018016
    R = 8.3145
    return u / (R * T)


def drhodT(p, T):
    u = 0.018016
    R = 8.3145
    return -p * u / (R * T ** 2)
