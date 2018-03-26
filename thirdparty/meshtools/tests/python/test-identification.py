import numpy as np
import MeshTools as MT

Point = MT.Point
Triangle = MT.Triangle
Tet = MT.Tetrahedron
Wedge = MT.Wedge

pts = [
	Point(( 0., 0., -1.)),
	Point(( 1., -1., 0.)),
	Point(( 1., 1., 0.)),
	Point(( -1., 0., 0.)),
	Point(( 0., 0., 1.))
]

pts2 = [
	Point(( 2., 0., -1.)),
	Point(( 2., -1., 0.)),
	Point(( 2., 1., 0.))
]

tets = [
	Tet((1, 2, 3, 0)),
	Tet((1, 2, 3, 4))
]

wedges = [
    Wedge((0, 1, 2, 5, 6, 7))
]

mesh =MT.TetMesh.Mesh()
vertices = mesh.vertices
for P in pts:
    vertices.append(P)

cellnodes = mesh.connectivity.cells.nodes
for elt in tets:
    cellnodes.append(elt)

mesh.connectivity.update_from_cellnodes()

print("face id:", mesh.face_id(Triangle((1, 2, 3))))
