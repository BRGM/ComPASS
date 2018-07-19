from collections import namedtuple
import numpy as np
import MeshTools as MT
import MeshTools.vtkwriters as vtkw
import MeshTools.URG as urg

Point = urg.Point
Vector = urg.Vector
Plane = urg.Plane
Hull = urg.Hull

def extract_data(filename):
    data = np.loadtxt(filename, skiprows=1, dtype=np.double)
    index = np.asarray(data[:, 0], dtype=np.int)
    return data[:, 1:], index

normals, ni = extract_data("URG_nvecs.out")
origins, oi = extract_data("URG_origins.out")

assert np.all(ni==oi)

xmin, ymin, zmin = origins.min(axis=0)
xmax, ymax, zmax = origins.max(axis=0)

Lh, Lv = 1E4, 1E3

# the hull is defined by the intersection of the negative side of an arbitrary
# set of planes
hull = Hull([
    Plane(Point(xmin - Lh, ymin, zmin), Vector(-1, 0, 0)),
    Plane(Point(xmax + Lh, ymin, zmin), Vector( 1, 0, 0)),
    Plane(Point(xmin, ymin - Lh, zmin), Vector( 0,-1, 0)),
    Plane(Point(xmin, ymax + Lh, zmin), Vector( 0, 1, 0)),
    Plane(Point(xmin, ymin, zmin - Lv), Vector( 0, 0,-1)),
    Plane(Point(xmin, ymin, zmax + Lv), Vector( 0, 0, 1)),
])

for bi, boundary in enumerate(hull.boundaries()):
    urg.simplify(boundary)
    #mesh = MT.TSurf.make(*boundary.as_arrays())
    #MT.to_vtu(mesh, "boundary%02d.vtu" % bi)

PlaneDefinition = namedtuple('PlaneDefinition', ['origin', 'normal'])   
plane_definitions = [PlaneDefinition(Point(*Oi), Vector(*ni))
                     for Oi, ni in zip(origins, normals)]
planes = [Plane(*definition) for definition in plane_definitions]

# select subset of planes
n = 5
planes = planes[:n]

# output normals as vtu
vtkw.write_vtu(
        vtkw.points_as_vtu_doc(
                origins[:n],
                pointdata = {'normals': normals[:n]},     
    ),
    'normals.vtu'
)

meshes = [hull.intersect(plane) for plane in planes]
for mesh in meshes:
    mesh.simplify_connected_components()

n = len(planes)
for i in range(n):
    print(i)
    Pi = planes[i]
    Mi = meshes[i]
    for j in range(i+1, n):
        Oj = plane_definitions[j].origin
        Mj = meshes[j]
        urg.corefine(Mi, Mj)
        Mi.simplify_connected_components()
        Mj.simplify_connected_components()
        Mj.keep(urg.On_side(Pi, Oj))

for mi, surface in enumerate(meshes):
    mesh = MT.TSurf.make(*surface.as_arrays())
    MT.to_vtu(mesh, "surface%03d.vtu" % mi,
        celldata = {'component': surface.connected_components()}
    )


