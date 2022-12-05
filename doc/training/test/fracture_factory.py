import numpy as np


class AASFracture:
    """
    Axis Aligned Square Fracture
    """

    def __init__(self, axis, center, width):
        center = np.array(center, copy=True)
        center.shape = -1
        dim = len(center)
        assert dim > 0
        assert axis in list(range(dim)), "wrong axis: 0:Ox 1:Oy 2:Oz..."
        assert width > 0
        self.dim = dim
        self.axis = axis
        self.hw = 0.5 * float(width)
        self.center = center
        self.projected_center = self._remove_axis(center)

    def _remove_axis(self, a):
        axis = self.axis
        return np.hstack([a[..., :axis], a[..., axis + 1 :]])

    def __call__(self, xyz, threshold=1.0e-10):
        assert threshold >= 0
        xyz = np.asarray(xyz)
        axis = self.axis
        xi = self.center[axis]
        if self.dim == 1:
            return np.abs(xyz - xi) < self.hw + threshold
        xyz.shape = -1, self.dim
        on_axis = np.abs(xyz[..., axis] - xi) < threshold
        on_fracture = (
            np.linalg.norm(
                self._remove_axis(xyz[on_axis]) - self.projected_center,
                ord=np.inf,
                axis=-1,
            )
            < self.hw + threshold
        )
        on_axis[on_axis] = on_fracture
        return on_axis


class AASFractureNetwork:
    def __init__(
        self, nb_fractures, origin, shape, steps, max_width_in_steps, seed=12345
    ):
        assert nb_fractures >= 1
        assert max_width_in_steps >= 1
        assert all([nx > 0 for nx in shape])
        assert all([dx > 0 for dx in steps])
        rng = np.random.default_rng(seed)
        rint = rng.integers
        assert len(origin) == len(shape) == len(steps)
        dim = len(origin)
        fractures = []
        for _ in range(nb_fractures):
            axis = rint(dim)
            center = [Ox + rint(nx) * dx for Ox, nx, dx in zip(origin, shape, steps)]
            width = rint(1, max_width_in_steps + 1) * steps[axis]
            fractures.append(AASFracture(axis, center, width))
        self.fractures = fractures

    def __call__(self, points, threshold=1.0e-10):
        fractures = self.fractures
        assert len(fractures) >= 1
        result = fractures[0](points, threshold=threshold)
        for F in fractures[1:]:
            result = result | F(points, threshold=threshold)
        return result


if __name__ == "__main__":

    F = AASFracture(0, 1.0, 1.0)
    print(F(0.0))
    print(F(1.0))
    print(F(2.0))
    print(F(np.linspace(0, 2)))

    nx, ny = 11, 11
    points = np.stack(
        np.meshgrid(np.linspace(0, 1, nx), np.linspace(0, 1, ny)), axis=-1
    )
    points.shape = nx * ny, 2

    for axis in [0, 1]:
        F = AASFracture(axis, (0.5, 0.5), 0.6)
        where = F(points)
        where.shape = nx, ny
        print(where.astype(np.int8))

    FN = AASFractureNetwork(
        10, (0, 0), (nx, ny), (1 / (nx - 1), 1 / (ny - 1)), 6, seed=1
    )
    where = FN(points)
    where.shape = nx, ny
    print(where.astype(np.int8))
