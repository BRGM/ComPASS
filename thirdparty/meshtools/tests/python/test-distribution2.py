import numpy as np
from mpi4py import MPI

import GridTools as GT
import MeshTools as MT

nx, ny, nz = 4, 4, 1

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    nbparts = comm.size
    mesh = MT.HexMesh.make(*GT.grid2hexs((nx, ny, nz), idtype=MT.idtype()))
    cell_color = np.zeros(mesh.nb_cells, dtype='i')
    for i in range(1, nbparts):
        cell_color[i*(mesh.nb_cells//nbparts):]+= 1
    distribution = mesh.distribute(nbparts, cell_color)
    node_color = np.empty(mesh.nb_vertices, dtype='i')
    face_color = np.empty(mesh.nb_faces, dtype='i')
    face_location = []
    for proc, dist in enumerate(distribution):
        cells, nodes, faces = [v.array_view() for v in dist[0]]
        ghost_cells_index, ghost_nodes_index, ghost_faces_index = dist[1]
        node_color[nodes[:ghost_nodes_index]] = proc 
        face_color[faces[:ghost_faces_index]] = proc 
        face_location.append(mesh.locate_faces_with_cell(cells, faces))
    cellnodes = mesh.connectivity.cells.nodes.raw_array()
    for proc in range(1, nbparts):
        cells, nodes, faces = [v.array_view() for v in distribution[proc][0]]
        comm.send((cells.shape[0], nodes.shape[0], faces.shape[0], cellnodes.shape[1]), dest=proc, tag=11)
        comm.Send([cellnodes[cells], MPI.BYTE], dest=proc, tag=12)
else:
    nbcells, nbnodes, nbfaces, raw_cell_size = comm.recv(source=0, tag=11)
    cellnodes = np.empty((nbcells, raw_cell_size), dtype=np.byte)
    comm.Recv([cellnodes, MPI.BYTE], source=0, tag=12)
    cellnodes = MT.HexMesh.CellsNodes.from_raw_array(cellnodes)
    for cell in cellnodes:
        print('proc', rank, 'received:', cell)
