from itertools import product
import numpy as np


def check_derivatives(f, dfdh, h, dh=1e-7, atol=1e-6, rtol=1e-5, output_maxerror=False):
    a = (f(h + dh) - f(h - dh)) / (2 * dh)
    b = dfdh(h)
    if output_maxerror:
        print(np.fabs(a - b).max())
    return np.allclose(a, b, atol=atol, rtol=rtol)


def check_partial_derivatives(
    f, dfdX, h, dh=1e-7, atol=1e-6, rtol=1e-5, output_error=False
):
    n = len(h)
    assert len(dfdX) == n
    assert dh > 0
    for i, dfdXi in enumerate(dfdX):
        dx = np.zeros(n, dtype=np.double)
        dx[i] = dh
        for X in product(*h):
            x = np.array(X)
            a = (f(*(x + dx)) - f(*(x - dx))) / (2 * dh)
            b = dfdXi(*x)
            if abs(a - b) > atol + rtol * abs(b):
                if output_error:
                    print(f"h={X} i={i}: {a} vs. {b}")
                return False
    return True
