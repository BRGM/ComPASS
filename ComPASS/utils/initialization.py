import types
import numpy as np


def hydrostatic_pressure_profile(
    simulation, zbottom, ztop, nz, ptop, geotherm, nbsteps=100
):
    assert zbottom < ztop
    rho = simulation.liquid_molar_density
    gravity = simulation.get_gravity()
    if not callable(geotherm):
        geotherm = lambda _: geotherm
    if isinstance(geotherm, types.FunctionType):
        geotherm = np.vectorize(geotherm)
    z = np.linspace(zbottom, ztop, nz)[::-1]  # from top to bottom
    p = ptop
    pressures = [p]
    for zbot, ztop in zip(z[1:], z[:-1]):
        zeta = ztop
        dz = (ztop - zbot) / nbsteps
        for _ in range(nbsteps):
            p += gravity * rho(p, geotherm(zeta)) * dz
            zeta -= dz
        pressures.append(p)
    pressures = np.asarray(pressures)
    return lambda zeta: np.interp(zeta, z[::-1], pressures[::-1])
