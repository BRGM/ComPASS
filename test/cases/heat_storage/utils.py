import numpy as np


def hydrostatic_pressure(zbottom, ztop, ptop, rho, nz, gravity=9.81, nbsteps=100):
    """
    :param zbottom:
    :param ztop: 
    :param ptop: pressure at ptop 
    :param rho: a function of p and z that gives the fluid density
    :param nz: the number of interpolation points
    :param gravity: acceleration of gravity (defaults to 9.81)
    :param nbsteps: the number of integration points between interpolation points 
    """
    assert zbottom < ztop
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    p = ptop
    pressures = [p]
    for zbot, ztop in zip(z[1:], z[:-1]):
        zeta = ztop
        dz = (ztop - zbot) / nbsteps
        for _ in range(nbsteps):
            p += gravity * rho(p, zeta) * dz
            zeta -= dz
        pressures.append(p)
    pressures = np.asarray(pressures)
    return lambda zeta: np.interp(zeta, z[::-1], pressures[::-1])
