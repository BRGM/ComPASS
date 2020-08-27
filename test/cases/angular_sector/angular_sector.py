import numpy as np

import pyCGAL.kernels.Epick.Mesh_2 as Mesh_2
from pyCGAL.kernels.Epick import Point_2 as Point
from pyCGAL.kernels.Epick import distance

from MeshTools import Wedge, HybridMesh
from MeshTools.utils import axis_extrusion


def mesh_2(R=1000, theta=np.pi / 6, ds=None, clean_mesh=True):
    """
    Mesh an horizontal angular sector and extrude it in 3D.
    :param R: radius of the angular sector
    :param theta: angle of the angular sector (in radians)
    :ds: may be a float with the target edge size or a pair with
         the target edge size around the origin and at the boundary
    :returns: a MeshTools mesh object
    """

    if ds is None:
        ds = 0.01 * R, 0.1 * R
    try:
        ds_O, ds_R = ds
    except TypeError:
        ds_O = ds_R = float(ds)

    cdt = Mesh_2.Constrained_Delaunay_triangulation_2()
    O = Point(0, 0)
    vO = cdt.insert(O)
    n = int((R * theta) / ds_R) + 1
    vAi = cdt.insert(Point(R, 0))
    cdt.insert_constraint(vO, vAi)
    for th in np.linspace(0, theta, n)[1:]:
        vAj = cdt.insert(Point(R * np.cos(th), R * np.sin(th)))
        cdt.insert_constraint(vAi, vAj)
        vAi = vAj
    cdt.insert_constraint(vAi, vO)

    def sizing_field(P):
        r = distance(O, P)
        return ((R - r) / R) * ds_O + (r / R) * ds_R

    if ds_O == ds_R:
        sizing_field = lambda P: ds_O

    criteria = Mesh_2.Delaunay_mesh_adaptative_size_criteria_2(
        sizing_field=sizing_field
    )
    Mesh_2.refine_Delaunay_mesh_2(cdt, criteria)

    Mesh_2.lloyd_optimize_mesh_2(cdt)

    vertices, triangles = Mesh_2.as_arrays(cdt)

    if clean_mesh:  # FIXME: we should use CGAL Epeck kernel to generate the mesh
        threshold = (0.1 * ds_O) ** 2
        x, y = [vertices[:, j] for j in range(2)]
        A, B, C = [triangles[:, j] for j in range(3)]
        area = 0.5 * np.fabs(
            (x[B] - x[A]) * (y[C] - y[A]) - (y[B] - y[A]) * (x[C] - x[A])
        )
        kept_triangles = np.nonzero(area > threshold)[0]
        kept_vertices, new_triangles = np.unique(
            triangles[kept_triangles], return_inverse=True
        )
        new_triangles.shape = -1, 3
        vertices, triangles = vertices[kept_vertices], new_triangles

    return vertices, triangles


def extruded_sector(R=1000, theta=np.pi / 6, ds=None, layer_thicknesses=None):
    """
    Mesh an horizontal angular sector and extrude it in 3D.
    :param R: radius of the angular sector
    :param theta: angle of the angular sector (in radians)
    :ds: may be a float with the target edge size or a pair with
         the target edge size around the origin and at the boundary
    :layer_thicknesses: a sequence of layer thicknesses
    :returns: a MeshTools mesh object
    """

    vertices, triangles = mesh_2(R, theta, ds)

    if layer_thicknesses is None:
        layer_thicknesses = [1.0]

    v3D = np.zeros((vertices.shape[0], 3), dtype="d")
    v3D[:, :2] = vertices

    vertices, wedges = axis_extrusion(v3D, triangles, layer_thicknesses)
    wedges = [Wedge(nodes) for nodes in wedges]
    return HybridMesh.create(vertices, wedges)
