import numpy as np
import MeshTools as MT

def test_meshtools():

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

    mesh = MT.TetMesh.Mesh()
    vertices = mesh.vertices
    for P in pts:
        vertices.append(P)

    cellnodes = mesh.connectivity.cells.nodes
    for elt in tets:
        cellnodes.append(elt)

    mesh.connectivity.update_from_cellnodes()

    print('nb vertices:', mesh.nb_vertices)
    print(mesh.vertices)

    print('connectivity:', mesh.connectivity)
    print('cells:', mesh.connectivity.cells)
    cellnodes = mesh.connectivity.cells.nodes
    print('nodes:', cellnodes)

    print('mesh faces:', mesh.connectivity.faces)
    for face in mesh.connectivity.faces.nodes:
        print(face)

    v = MT.as_coordinate_array(mesh.vertices)
    print(mesh.vertices[0])
    print(v)

    print("face centers:")
    fc = mesh.face_centers()
    print(fc)
    print(MT.as_coordinate_array(fc))

    print("cells nodes as COC data")
    pointers, nodes = mesh.cells_nodes_as_COC()
    print(pointers)
    print(nodes)
    print("cells faces as COC data")
    pointers, faces = mesh.cells_faces_as_COC()
    print(pointers)
    print(faces)

    mesh = MT.HybridMesh.Mesh()
    vertices = mesh.vertices
    for P in pts + pts2:
        vertices.append(P)

    cellnodes = mesh.connectivity.cells.nodes
    for elt in tets + wedges:
        cellnodes.append(elt)

    mesh.connectivity.update_from_cellnodes()
    print('hybrid mesh cells:')
    for cell in mesh.connectivity.cells.nodes:
        print(cell)
    print('hybrid mesh faces:')
    for fi, face in enumerate(mesh.connectivity.faces.nodes):
        print('face', fi, ':' , face,
              'check', mesh.connectivity.faces.id(face),
              'type', type(face))
    print('face centers:')
    fc = mesh.face_centers()
    print(MT.as_coordinate_array(fc))

    boundaries = mesh.connectivity.boundary_faces()
    print('boundary faces:', boundaries)

    cells_types = set(type(elt) for elt in tets + wedges)
    print('Mesh has', len(cells_types), 'different types.')
    print(set(cells_types))

    print("cells nodes as COC data")
    pointers, nodes = mesh.cells_nodes_as_COC()
    print(pointers)
    print(nodes)
    print("faces nodes as COC data")
    pointers, nodes = mesh.faces_nodes_as_COC()
    print(pointers)
    print(nodes)

if __name__=='__main__':
    test_meshtools()
