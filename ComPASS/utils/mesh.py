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


if __name__ == "__main__":
    import vtkwriters as vtkw

    mesh = extrude([(0, 0), (1, 0), (0, 1), (1, 1)], [(0, 1, 3, 2)])
    vtkw.write_vtu(vtkw.vtu_doc(mesh.vertices, mesh.cell_nodes), "test")
    vtkw.write_vtp(vtkw.vtp_doc(mesh.vertices, mesh.face_nodes), "test-faces")
