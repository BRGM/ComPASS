import numpy as np
import MeshTools as MT


def test_meshtool():

    Point = MT.Point
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

    mesh = MT.HybridMesh.Mesh()
    vertices = mesh.vertices
    for P in pts + pts2:
        vertices.append(P)

    cellnodes = mesh.connectivity.cells.nodes
    for elt in tets + wedges:
        cellnodes.append(elt)

    mesh.connectivity.update_from_cellnodes()

    MT.to_vtu(mesh, 'mesh-test3.vtu')

if __name__=='__main__':
    test_meshtools()
