import numpy as np
import MeshTools as MT
import MeshTools.vtkwriters as vtkw

def test_hexaedron():

    Point = MT.Point
    Hexa = MT.Hexahedron

    pts = [
	    Point(( 0., 0., 0.)),
	    Point(( 1., 0., 0.)),
	    Point(( 1., 1., 0.)),
	    Point(( 0., 1., 0.)),
	    Point(( 0., 0., 1.)),
	    Point(( 1., 0., 1.)),
	    Point(( 1., 1., 1.)),
	    Point(( 0., 1., 1.)),
    ]

    hexas = [
	    Hexa((0, 1, 2, 3, 4, 5, 6, 7)),
    ]

    mesh = MT.HybridMesh.Mesh()
    vertices = mesh.vertices
    for P in pts:
        vertices.append(P)

    cellnodes = mesh.connectivity.cells.nodes
    for elt in hexas:
        cellnodes.append(elt)

    mesh.connectivity.update_from_cellnodes()

    offsets, cellsnodes = mesh.cells_nodes_as_COC()
    vtkw.write_vtu(
        vtkw.vtu_doc_from_COC(
                mesh.vertices_array(), 
                np.array(offsets[1:], copy=False), # no first zero offset for vtk 
                np.array(cellsnodes, copy=False),
                mesh.cells_vtk_ids()
        ),
        'hexa.vtu'
    ) 

if __name__=='__main__':
    test_hexaedron()
