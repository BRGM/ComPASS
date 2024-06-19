import numpy as np
import meshio
from ComPASS.utils.mesh import extrude


def convert_mesh(filename, thickness=None):
    mesh = meshio.read(filename)
    points = np.copy(mesh.points)
    # points are already 3D but coordinate along 3rd dimension is 0
    # y <-> z # use advance slicing: https://stackoverflow.com/a/4857981
    points[:, [1, 2]] = points[:, [2, 1]]
    if thickness is None:
        z = np.unique(points[:, 2])
        dz = z[1] - z[0]
        thickness = dz
    assert thickness > 0
    raw_mesh = extrude(points, mesh.cells[0].data, thickness=thickness, axis=1)
    # specify the vtk format of the cell (define the cell nodes order)
    nb_cells = len(raw_mesh.cell_nodes)
    raw_mesh.cell_types = np.full((nb_cells,), 12, dtype=int)
    rocktype = np.array(mesh.cell_data["gmsh:physical"])[0]
    # rocktype[rocktype == 7] = 0
    # meshio.write_points_cells(
    # "cells_rkt_RawMesh.vtu",
    # points=mesh.points,
    # cells=mesh.cells,
    # cell_data={"rocktype": [rocktype]},
    # )
    # print(f"VTU file with rock types written to cells rtk")
    return raw_mesh, rocktype


if __name__ == "__main__":
    import sys
    import vtkwriters as vtkw

    if len(sys.argv) > 1:

        filename = sys.argv[1]
        thickness = 15.0

        mesh, rocktype = convert_mesh(filename, thickness)

        vtkw.write_vtu(
            vtkw.vtu_doc(
                mesh.vertices, mesh.cell_nodes, celldata={"rocktype": rocktype}
            ),
            "test",
        )
