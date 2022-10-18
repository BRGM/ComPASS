import numpy as np
from ComPASS.utils.dbgtools import check_derivatives, check_partial_derivatives
from ComPASS.utils.units import *

p = np.logspace(np.log10(bar), np.log10(50 * MPa), 10)
T = np.linspace(273.15, 573.15, 10)


def test_diphasic():
    from ComPASS.properties.water2ph.diphasic import psat, dpsatdT, Tsat, dTsatdp

    assert check_derivatives(psat, dpsatdT, T)
    assert check_derivatives(Tsat, dTsatdp, p)


def test_phases():
    from ComPASS.properties.water2ph import gas, liquid

    def check(phase):
        assert check_partial_derivatives(
            phase.rho, [phase.drhodp, phase.drhodT], [p, T], output_error=True
        ), "error in rho derivatives"
        assert check_partial_derivatives(
            phase.mu, [phase.dmudp, phase.dmudT], [p, T], output_error=True
        ), "error in mu derivatives"
        assert check_partial_derivatives(
            phase.h, [phase.dhdp, phase.dhdT], [p, T], output_error=True
        ), "error in h derivatives"
        pass

    check(gas)
    check(liquid)
