import numpy as np
import MeshTools as MT
import MeshTools.URG as urg

Point = urg.Point
Vector = urg.Vector
Plane = urg.Plane
Hull = urg.Hull

hull = Hull([
    Plane(Point(-1, 0, 0), Vector(-1, 0, 0)),
    Plane(Point(1, 0, 0), Vector(1, 0, 0)),
    Plane(Point(0, -1, 0), Vector(0, -1, 0)),
    Plane(Point(0, 1, 0), Vector(0, 1, 0)),
    Plane(Point(0, 0, -1), Vector(0, 0, -1)),
    Plane(Point(0, 0, 1), Vector(0, 0, 1)),   
])

for bi, boundary in enumerate(hull.boundaries()):
    mesh = MT.TSurf.make(*boundary.as_arrays())
    MT.to_vtu(mesh, "hull_boundary_%02d.vtu" % bi)


P = Plane(Point(0, 0, 0), Vector(1, 1.5, 1))
S = hull.intersect(P)
mesh = MT.TSurf.make(*S.as_arrays())

MT.to_vtu(mesh, "intersection.vtu")

