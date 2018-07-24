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
n = 20
planes = planes[:n]

# output normals as vtu
vtkw.write_vtu(
        vtkw.points_as_vtu_doc(
                origins[:n],
                pointdata = {'normals': normals[:n]},     
    ),
    'normals.vtu'
)

surfaces = [hull.intersect(plane) for plane in planes]
for S in surfaces:
    S.mark_all_vertices_as_corners()
    S.simplify_connected_components()

for k, S in enumerate(surfaces):
    mesh = MT.TSurf.make(*S.as_arrays())
    MT.to_vtu(mesh, "original_surface%03d.vtu" % k,
        celldata = {'component': S.connected_components()}
    )

def dump(k):
    print('dumping', k)
    Sk = surfaces[k]
    print('   constraints:', Sk.constraints())
    np.savetxt('constrained_midpoints%03d.txt' % k, Sk.constrained_midpoints())
    mesh = MT.TSurf.make(*Sk.as_arrays())
    MT.to_vtu(mesh, "surface%03d.vtu" % k,
        celldata = {'component': Sk.connected_components()}
    )
    print('dumped', k)

n = len(planes)
for j in range(1, n):
    print(j)
    Oj = plane_definitions[j].origin
    Sj = surfaces[j]
    for i in range(j):
        #if j>=9:
        print('#'*10, i, '<->', j)
        Pi = planes[i]
        Si = surfaces[i]
        #if i == 0 and j>=17 :
        #    for k in (j, i):
        #        dump(k)
        #if (i, j) == (0, 14) :
        #    for k in (6, i, j):
        #        dump(k)
        #    1/0
        nb_intersection_edges = urg.corefine(Si, Sj)
        if nb_intersection_edges>0:
            #print('>>>>>>>>> (', Si.index, Sj.index, ')', Si.number_of_constraints(Sj), 'vs', Sj.number_of_constraints(Si))
            print('Constraints', Si.index, ':', Si.constraints())
            print('           ', Sj.index, ':', Sj.constraints())
            #if j>=9:
            #    meshj = MT.TSurf.make(*Sj.as_arrays())
            #    MT.to_vtu(meshj, "surface%03d.vtu" % j,
            #        celldata = {'component': Sj.connected_components()}
            #    )
            #if i==4 and j==9:
            #    relaxed = []
            #else:
            #    relaxed = Sj.keep_connected(urg.On_side(Pi, Oj))
            relaxed = Sj.keep_connected(urg.On_side(Pi, Oj))
            for sk in relaxed:
                #print('relaxing', surfaces[sk].index, 'from', Sj.index)
                #surfaces[sk].remove_constraint(Sj.index)
                #if(j>=16 and surfaces[sk].index==0):
                #    dump(6)
                #    dump(0)
                print('simplifying relaxed surface', surfaces[sk].index, 'from', Sj.index)
                surfaces[sk].simplify_connected_components()
                if(surfaces[sk].index==0):
                    dump(sk)
            #if True or j<6:
            print('simplifying', Si.index)
            Si.simplify_connected_components()
            if(i==0):
                dump(i)
            #if j>=9:
            print('simplifying', Sj.index)
            Sj.simplify_connected_components()
        #if (i, j) == (6, 14) :
        #    for k in (0, i, j):
        #        dump(k)
        #    1/0
    if j>=16:
        dump(0)
for k, S in enumerate(surfaces):
    np.savetxt('constrained_midpoints%03d.txt' % k, S.constrained_midpoints())
    mesh = MT.TSurf.make(*S.as_arrays())
    MT.to_vtu(mesh, "surface%03d.vtu" % k,
        celldata = {'component': S.connected_components()}
    )


