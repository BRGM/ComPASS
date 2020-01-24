from collections import namedtuple
import numpy as np


class Solution(namedtuple("params", ["Pw", "qw", "mu", "km", "rho", "rw"])):
    def __call__(self, x, y):
        Pw, qw, mu, km, rho, rw = self
        r = np.sqrt(x ** 2 + y ** 2)
        res = Pw * np.ones(r.shape)
        mask = r > rw
        res[mask] = Pw + ((qw * mu) / (2 * np.pi * km * rho)) * np.log(r[mask] / rw)
        return res


if __name__ == "__main__":
    S = Solution(0, 1, 1, 1, 1, 0.1)
    x = np.linspace(1e-3, 10)
    assert np.allclose(S(x, 0), S(0, x))
