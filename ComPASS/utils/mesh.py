import numpy as np
from ComPASS import RawMesh


def _concatenate_packs(packs):
    result = []
    for pack in packs:
        result.extend(pack)
    return result


def extrude(points, cells, thickness=1.0, axis=2):

    points = np.asarray(points)
    assert points.ndim == 2
    cells = np.asarray(cells)
    assert cells.ndim == 2

    if points.shape[1] == 2:
        if axis == 0:
            points = np.hstack([np.zeros(points.shape[0])[:, None], points])
        elif axis == 2:
            points = np.hstack([points, np.zeros(points.shape[0])[:, None]])
        else:
            assert axis == 1
            points = np.hstack(
                [
                    points[:, 0][:, None],
                    np.zeros(points.shape[0])[:, None],
                    points[:, 1][:, None],
                ]
            )
    assert points.shape[1] == 3

    vertices = np.vstack([points, points])
    npts = points.shape[0]
    vertices[npts:, axis] = vertices[:npts, axis] + thickness

    nc = cells.shape[0]  # number of (2D) cells
    col = lambda j: cells[:, j][:, None]

    face_nodes = [
        cells,  # front
        cells + npts,  # back
    ]
    for j in range(cells.shape[1]):
        face_nodes.append(
            np.hstack([col(j - 1), col(j - 1) + npts, col(j) + npts, col(j)])
        )
    if cells.shape[1] == 4:  # all quads
        face_nodes = np.vstack(face_nodes)
    else:
        face_nodes = _concatenate_packs(face_nodes)

    cell_faces = np.hstack(
        [np.arange(k * nc, (k + 1) * nc)[:, None] for k in range(2 + cells.shape[1])]
    )
    cell_nodes = np.hstack([cells, cells + npts])

    return RawMesh(
        vertices=vertices,
        face_nodes=face_nodes,
        cell_nodes=cell_nodes,
        cell_faces=cell_faces,
    )


def extrude_triangle_strips(Lx, nx, dz, nbstrips=1):
    """
    :param Lx: length of the strips along Ox
    :param nx: number of points along Ox
    :param dz: height of the strips
    :param nbstrips: number of strips
    """
    assert nbstrips > 0
    nz = nbstrips + 1
    assert nx >= 2 and nz >= 2  # 2 points for at least one cell
    x = np.linspace(0, Lx, nx)
    dx = x[1] - x[0]
    assert dx > 0
    # x staggered
    xs = np.hstack([x[0], x[1:] - 0.5 * dx, x[-1]])
    pts = []
    for k in range(nz):
        if k % 2 == 0:
            pts.append(np.vstack([x, np.tile(k * dz, nx)]).T)
        else:
            pts.append(np.vstack([xs, np.tile(k * dz, nx + 1)]).T)
    pts = np.vstack(pts)
    triangles = []
    for k in range(nz - 1):
        O = k * nx + k // 2
        if k % 2 == 0:
            triangles.append(
                np.vstack(
                    [
                        np.array([[O, O + nx, O + nx + 1], [O, O + nx + 1, O + 1]]) + i
                        for i in range(nx - 1)
                    ]
                    + [[O + nx - 1, O + 2 * nx - 1, O + 2 * nx]]
                )
            )
        else:
            triangles.append(
                np.vstack(
                    [
                        np.array(
                            [[O, O + nx + 1, O + 1], [O + nx + 1, O + nx + 2, O + 1]]
                        )
                        + i
                        for i in range(nx - 1)
                    ]
                    + [[O + nx - 1, O + 2 * nx, O + nx]]
                )
            )
    triangles = np.vstack(triangles)
    return extrude(pts, triangles, axis=1, thickness=dz)


def extrude_quad_strips(Lx, nx, dz, nbstrips=1):
    """
    :param Lx: length of the strips along Ox
    :param nx: number of points along Ox
    :param dz: height of the strips
    :param nbstrips: number of strips
    """
    assert nbstrips > 0
    nz = nbstrips + 1
    assert nx >= 2 and nz >= 2  # 2 points for at least one cell
    x = np.linspace(0, Lx, nx)
    dx = x[1] - x[0]
    assert dx > 0
    pts = []
    for k in range(nz):
        pts.append(np.vstack([x, np.tile(k * dz, nx)]).T)
    pts = np.vstack(pts)
    row = np.vstack([np.array([0, 1, nx + 1, nx]) + k for k in range(nx - 1)])
    quads = np.vstack([row + k * nx for k in range(nz - 1)])
    return extrude(pts, quads, axis=1, thickness=dz)


if __name__ == "__main__":
    import vtkwriters as vtkw

    mesh = extrude([(0, 0), (1, 0), (0, 1), (1, 1)], [(0, 1, 3, 2)])
    vtkw.write_vtu(vtkw.vtu_doc(mesh.vertices, mesh.cell_nodes), "test")
    vtkw.write_vtp(vtkw.vtp_doc(mesh.vertices, mesh.face_nodes), "test-faces")
