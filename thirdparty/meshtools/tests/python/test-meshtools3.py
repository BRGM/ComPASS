import numpy as np
import MeshTools

Point = MeshTools.Point
Tet = MeshTools.Tetrahedron
Wedge = MeshTools.Wedge

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

mesh = MeshTools.HybridMesh.Mesh()
vertices = mesh.vertices
for P in pts + pts2:
    vertices.append(P)

cellnodes = mesh.connectivity.cells.nodes
for elt in tets + wedges:
    cellnodes.append(elt)

mesh.connectivity.update_from_cellnodes()

MeshTools.to_vtu(mesh, 'mesh-test3.vtu')

#data = mesh.connectivity.faces.nodes
#
#from mpi4py import MPI
#
#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#
#if rank == 0:
#    req = comm.Isend(data, dest=1, tag=11)
#    req.wait()
#elif rank == 1:
#    req = comm.Irecv(source=0, tag=11)
#    data = req.wait()
#    print('from proc 1', type(data))
#



#print('hybrid mesh cells:')
#for cell in mesh.connectivity.cells.nodes:
#    print(cell)

#print('hybrid mesh faces:')
#for fi, face in enumerate(mesh.connectivity.faces.nodes):
#    print('face', fi, ':' , face,
#          'check', mesh.connectivity.faces.id(face),
#          'type', type(face))
#print('face centers:')
#fc = mesh.face_centers()
#print(MeshTools.as_coordinate_array(fc))

#boundaries = mesh.connectivity.boundary_faces()
#print('boundary faces:', boundaries)

#cells_types = set(type(elt) for elt in tets + wedges)
#print('Mesh has', len(cells_types), 'different types.')
#print(set(cells_types))

#print("cells nodes as COC data")
#pointers, nodes = mesh.cells_nodes_as_COC()
#print(pointers)
#print(nodes)
#print("faces nodes as COC data")
#pointers, nodes = mesh.faces_nodes_as_COC()
#print(pointers)
#print(nodes)
