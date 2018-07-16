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


P1 = Plane(Point(0, 0, 0), Vector(1, 1.5, 1))
P2 = Plane(Point(0, 0, -0.2), Vector(0.15, 0.1, 1))
S1 = hull.intersect(P1)
S2 = hull.intersect(P2)

urg.corefine(S1, S2)

on_side = urg.On_side(P1, Point(0,0, 0.5))
S2.keep(on_side)

for si, surface in enumerate([S1, S2]):
    mesh = MT.TSurf.make(*surface.as_arrays())
    MT.to_vtu(mesh, "surface%02d.vtu" % si)

