import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np

import mCGAL.kernels.Epick.Mesh_2 as Mesh_2
from mCGAL.kernels.Epick import Point_2 as Point
from mCGAL.kernels.Epick import Segment_2 as Segment
from mCGAL.kernels.Epick import distance

# C ------------- (L) ------------- D
# |                                 |
# |                                 |
# (h)                               (h)
# |                                 |
# |                                 |
# A ------ W1 --- (d) --- W2 ------ B

d = 1200  # inter well distance
L = 10000  # field length
h = 0.5 * (L - d)  # half field width

sw = 25.0  # cell size around wells
siw = 100.0  # cell size inter wells
sff = 1000.0  # cell size far field
lff = 2 * d  # far field limit

cdt = Mesh_2.Constrained_Delaunay_triangulation_2()
vA = cdt.insert(Point(-0.5 * L, 0))
vB = cdt.insert(Point(0.5 * L, 0))
vC = cdt.insert(Point(-0.5 * L, h))
vD = cdt.insert(Point(0.5 * L, h))
# cdt.insert(Point(0, 0.5*h))
W1 = Point(-0.5 * d, 0)
vW1 = cdt.insert(W1)
W2 = Point(0.5 * d, 0)
vW2 = cdt.insert(W2)
cdt.insert_constraint(vA, vW1)
cdt.insert_constraint(vW1, vW2)
cdt.insert_constraint(vW2, vB)
cdt.insert_constraint(vB, vD)
cdt.insert_constraint(vC, vD)
cdt.insert_constraint(vA, vC)

W1W2 = Segment(W1, W2)


def sizing_field(P):
    W1P = distance(W1, P)
    W2P = distance(W2, P)
    W1W2P = distance(W1W2, P)

    def f(r, s):
        if r > lff:
            return sff
        return ((lff - r) / lff) * s + (r / lff) * sff

    return min(f(W1P, sw), f(W2P, sw), f(W1W2P, siw))


criteria = Mesh_2.Delaunay_mesh_adaptative_size_criteria_2(sizing_field=sizing_field)
Mesh_2.refine_Delaunay_mesh_2(cdt, criteria)

Mesh_2.lloyd_optimize_mesh_2(cdt, max_iteration_number=100)

print(f"Number of vertices: {cdt.number_of_vertices()}")
print(f"Number of triangles: {cdt.number_of_faces()}")

vertices, triangles = Mesh_2.as_arrays(cdt)

# We reset border coordinates just in case....
epsilon = 1e-3 * sw
assert epsilon > 0
x, y = vertices[:, 0], vertices[:, 1]
vertices[x < x.min() + epsilon, 0] = -0.5 * L
vertices[x > x.max() - epsilon, 0] = 0.5 * L
vertices[y < epsilon, 1] = 0
vertices[y > y.max() - epsilon, 1] = h

# plt.gca().set_aspect("equal")
# plt.triplot(vertices[:, 0], vertices[:, 1], triangles)
# plt.plot([-0.5*d, 0.5*d], [0, 0], 'rx')
# plt.show()

# Go to 3D
tmp = vertices
vertices = np.empty((tmp.shape[0], 3), dtype=np.double)
vertices[:, :2] = tmp

import MeshTools as MT
from MeshTools.utils import axis_extrusion

thickness = 1.0
nb_layers = 3 * 3
vertices, cells = axis_extrusion(vertices, triangles, np.tile(thickness, nb_layers))
nt = triangles.shape[0]
assert cells.shape == (nb_layers * nt, 6)
reservoir = np.zeros(cells.shape[0], dtype=np.bool)
reservoir[nt : 2 * nt] = True
np.savez("mesh", vertices=vertices, cells=cells, reservoir=reservoir)

elements = [MT.Wedge(MT.idarray(cell)) for cell in cells]
mesh = MT.TetWedgeMesh.create(vertices, elements)
MT.to_vtu(mesh, "mesh.vtu", celldata={"reservoir": reservoir.astype(np.int)})
