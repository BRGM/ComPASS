import numpy as np

from MeshTools import Hexahedron, HybridMesh
from MeshTools.utils import axis_extrusion


def geometric_progression(n, u0, S, tol=1e-8):
    assert 0 < u0 < S
    assert n > 1
    R = lambda a: u0 * (a ** (n + 1) - 1) / (a - 1) - S
    dRda = lambda a: u0 * ((n * (a - 1) + 1) * a**n / (a - 1) ** 2)
    a = np.exp(np.log(S / u0) / n)
    assert a > 1
    Ra = R(a)
    while abs(Ra) > tol:
        a -= Ra / dRda(a)
        Ra = R(a)
    return a


def make_mesh(r, theta, layer_thickness):

    r = np.asarray(r, dtype="d")

    assert r.ndim == 1
    assert r.size > 1 and np.all(r[:-1] < r[1:])

    try:
        h = float(layer_thickness)
        layer_thickness = [
            h,
        ]
    except TypeError:
        pass

    def radius_vertices(alpha, z):
        res = np.vstack(
            [
                r * np.cos(alpha),
                r * np.sin(alpha),
                np.tile(z, r.shape[0]),
            ]
        )
        return res.T

    vertices = np.vstack(
        [
            radius_vertices(0, 0),
            radius_vertices(theta, 0),
        ]
    )

    nr = r.shape[0]
    ni = np.arange(0, nr - 1)
    ni.shape = -1, 1

    bottom_faces = np.hstack([ni, ni + 1, ni + nr + 1, ni + nr])

    vertices, hexs = axis_extrusion(vertices, bottom_faces, layer_thickness)
    hexs = [Hexahedron(nodes) for nodes in hexs]

    return HybridMesh.create(vertices, hexs)


def make_grid(r, dy, layer_thickness):

    r = np.asarray(r, dtype="d")

    assert r.ndim == 1
    assert r.size > 1 and np.all(r[:-1] < r[1:])

    try:
        h = float(layer_thickness)
        layer_thickness = [
            h,
        ]
    except TypeError:
        pass

    def base_vertices(y, z):
        res = np.vstack(
            [
                r,
                np.tile(y, r.shape[0]),
                np.tile(z, r.shape[0]),
            ]
        )
        return res.T

    vertices = np.vstack(
        [
            base_vertices(0, 0),
            base_vertices(dy, 0),
        ]
    )

    nr = r.shape[0]
    ni = np.arange(0, nr - 1)
    ni.shape = -1, 1

    bottom_faces = np.hstack([ni, ni + 1, ni + nr + 1, ni + nr])

    vertices, hexs = axis_extrusion(vertices, bottom_faces, layer_thickness)
    hexs = [Hexahedron(nodes) for nodes in hexs]

    return HybridMesh.create(vertices, hexs)


def make_logmesh(nr, rmin, rmax, theta, layer_thickness):

    assert nr > 1
    assert 0 < rmin < rmax

    return make_mesh(
        np.logspace(np.log10(rmin), np.log10(rmax), nr), theta, layer_thickness
    )


def make_geometric_mesh(n, u0, S, theta, layer_thickness, tol=1e-8):
    a = geometric_progression(n, u0, S, tol)
    r = [u0]
    for _ in range(n):
        r.append(r[-1] * a)
    return make_mesh(np.cumsum(r), theta, layer_thickness)


if __name__ == "__main__":
    from MeshTools import to_vtu

    to_vtu(make_logmesh(20, 0.1, 1e3, np.pi / 4.0, 10.0), "hextest")
    n, u0, S = 10, 0.1, 1000
    a = geometric_progression(n, u0, S)
    r = [u0]
    for _ in range(n):
        r.append(r[-1] * a)
    print(r)
    print(np.cumsum(r))
