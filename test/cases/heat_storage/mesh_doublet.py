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

np.savetxt("doublet.nodes", np.around(vertices, decimals=1), fmt="%g")
np.savetxt("doublet.triangles", triangles, fmt="%d")
