import numpy as np
import MeshTools as MT
import MeshTools.vtkwriters as vtkw

def print_collection(iterable):
    for item in iterable:
        print(item)

def explore_mesh(mesh, name):
    print('Nb vertices:', mesh.nb_vertices)
    print('Nb cells:', mesh.nb_cells)
    print('Nb faces:', mesh.nb_faces)
    print('Vertices')
    print_collection(mesh.vertices)
    print('Cell nodes')
    print_collection(mesh.connectivity.cells.nodes)
    print('Cell face')
    print_collection(mesh.connectivity.cells.faces)
    print('Face nodes')
    print_collection(mesh.connectivity.faces.nodes)
    print('Boundary faces:', mesh.boundary_faces())
    fcenters = MT.as_coordinate_array(mesh.face_centers())
    print('All face centers: (with shape', fcenters.shape, ')')
    print(fcenters)
    ccenters = MT.as_coordinate_array(mesh.cell_centers())
    print('All cell centers: (with shape', ccenters.shape, ')')
    print(ccenters)
    vertices = MT.as_coordinate_array(mesh.vertices)
    cellnodes = np.array([np.array(cell) for cell in mesh.connectivity.cells.nodes])
    print(cellnodes)
    vtkw.write_vtu(
        vtkw.vtu_doc(vertices, cellnodes),
        name + '.vtu'
    )
    vtkw.write_vtu(
        vtkw.vtu_doc(fcenters, np.reshape(np.arange(fcenters.shape[0]), (-1,1))),
        name + '_face_centers.vtu'
    )

def test_meshtools():

    vertices = np.array([(0,0,-1), 
        (0.5*np.sqrt(3), 0.5*np.sqrt(3), 0),
        (-1, 0, 0),
        (0.5*np.sqrt(3), -0.5*np.sqrt(3), 0),
        (0,0,1)], dtype=np.double)
    cells = MT.idarray([(0, 1, 2, 3), (1, 2, 3, 4)])

    mesh = MT.TetMesh.make(vertices, MT.idarray(cells))

    explore_mesh(mesh, 'tets')

    horizontal_facet = MT.Triangle((1, 2, 3))
    print('Horizontal facet:', mesh.connectivity.faces.id(horizontal_facet))

    vertices = np.array([
        (0,0,0), (1,0,0), (1,1,0), (0,1,0), 
        (0,0,1), (1,0,1), (1,1,1), (0,1,1), 
        (0,0,2), (1,0,2), (1,1,2), (0,1,2), 
    ], dtype=np.double)
    # Pisa tower...
    vertices[:, 0]+= 0.1 * vertices[:,2] 
    cells = MT.idarray([(0, 1, 2, 3, 4, 5, 6, 7), (4, 5, 6, 7, 8, 9, 10, 11)])

    mesh =  MT.HexMesh.make(vertices, MT.idarray(cells))

    explore_mesh(mesh, 'hexs')

if __name__=='__main__':
    test_meshtools()
