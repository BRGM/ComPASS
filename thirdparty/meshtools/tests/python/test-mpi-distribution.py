import numpy as np
from mpi4py import MPI

import GridTools as GT
import MeshTools as MT

import sys
import datetime
import vtkwriters as vtkw

nx, ny, nz = 5, 5, 6

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def np2mpi(nptype):
    return MPI._typedict[np.dtype(nptype).char]
mpi_id_type = np2mpi(MT.idtype())

def geometric_coloring(mesh):
    nbparts = comm.size
    cell_color = np.zeros(mesh.nb_cells, dtype='i')
    for i in range(1, nbparts):
        cell_color[i*(mesh.nb_cells//nbparts):]+= 1
    return cell_color

def create_and_distribute_mesh(mesh_constructor, constructor_arguments, coloring_algorithm):
    # retrieve mesh class from underlying MeshTools module
    mesh_module_name = mesh_constructor.__module__.split('.')[-1]
    assert hasattr(MT, mesh_module_name)
    start = datetime.datetime.now()
    if rank == 0:
        global_mesh = mesh_constructor(*constructor_arguments)
        print('coloring after', (datetime.datetime.now() - start).total_seconds())
        cell_color = coloring_algorithm(global_mesh)
        print('distributing after', (datetime.datetime.now() - start).total_seconds())
        distribution = global_mesh.distribute(cell_color)
        node_color = np.empty(global_mesh.nb_vertices, dtype='i')
        face_color = np.empty(global_mesh.nb_faces, dtype='i')
        face_location = []
        # we make a first loop to locate faces with relatively to cell
        for proc, dist in enumerate(distribution):
            cells, nodes, faces = [v.array_view() for v in dist[0]]
            #ghost_cells_index, ghost_nodes_index, ghost_faces_index = dist[1]
            #node_color[nodes[:ghost_nodes_index]] = proc 
            #face_color[faces[:ghost_faces_index]] = proc 
            face_location.append(global_mesh.locate_faces_with_cell(cells, faces))
        vertices = global_mesh.vertices_array()
        cellnodes = global_mesh.connectivity.cells.nodes.raw_array()
        print('starting distribution after', (datetime.datetime.now() - start).total_seconds())
        nbparts = len(distribution)
        assert nbparts==comm.size
        for proc in range(1, nbparts):
            cells, nodes, faces = [v.array_view() for v in distribution[proc][0]]
            # ?? send owns or deduce owns from color ? 
            comm.send((cells.shape[0], nodes.shape[0], faces.shape[0], cellnodes.shape[1]), dest=proc, tag=11)
            comm.send(distribution[proc][1], dest=proc, tag=111)
            # send global ids ->
            assert cells.dtype==MT.idtype()
            comm.Send([cells, mpi_id_type], dest=proc, tag=12)
            assert nodes.dtype==MT.idtype()
            comm.Send([nodes, mpi_id_type], dest=proc, tag=13) 
            assert faces.dtype==MT.idtype()
            comm.Send([faces, mpi_id_type], dest=proc, tag=14) 
            # send part elements ->
            assert vertices.dtype==np.double
            comm.Send([vertices[nodes], np2mpi(np.double)], dest=proc, tag=15) 
            comm.Send([cellnodes[cells], MPI.BYTE], dest=proc, tag=16)
            face_cell, face_position = face_location[proc]
            assert face_cell.dtype==MT.idtype()
            comm.Send([face_cell, mpi_id_type], dest=proc, tag=17) 
            assert face_position.dtype==np.uint8
            comm.Send([face_position, np2mpi(np.uint8)], dest=proc, tag=18) 
        # part that stays on master proc
        cells_gid, nodes_gid, unsorted_faces_gid = [v.array_view() for v in distribution[0][0]]
        ghost_indexes = distribution[0][1]
        vertices = np.copy(vertices[nodes_gid])
        cellnodes = cellnodes[cells_gid]
        face_cell, face_position = face_location[0]
    else:
        nbcells, nbnodes, nbfaces, raw_cell_size = comm.recv(source=0, tag=11)
        ghost_indexes = comm.recv(source=0, tag=111)
        # -> receive global ids
        cells_gid = np.empty((nbcells,), dtype=MT.idtype())
        comm.Recv([cells_gid, mpi_id_type], source=0, tag=12) 
        nodes_gid = np.empty((nbnodes,), dtype=MT.idtype())
        comm.Recv([nodes_gid, mpi_id_type], source=0, tag=13) 
        unsorted_faces_gid = np.empty((nbfaces,), dtype=MT.idtype())
        comm.Recv([unsorted_faces_gid, mpi_id_type], source=0, tag=14) 
        # -> receive part elements
        vertices = np.empty((nbnodes, 3), dtype=np.double)
        comm.Recv([vertices, np2mpi(np.double)], source=0, tag=15)
        #print('proc', comm.rank, ':', vertices.min(axis=0), vertices.max(axis=0))
        cellnodes = np.empty((nbcells, raw_cell_size), dtype=np.byte)
        comm.Recv([cellnodes, MPI.BYTE], source=0, tag=16)
        face_cell = np.empty((nbfaces,), dtype=MT.idtype())
        comm.Recv([face_cell, mpi_id_type], source=0, tag=17) 
        face_position = np.empty((nbfaces,), dtype=np.uint8)
        comm.Recv([face_position, np2mpi(np.uint8)], source=0, tag=18) 
    # Rebuild local meshes
    print('rebuilding on', rank, 'after', (datetime.datetime.now() - start).total_seconds())
    Mesh = getattr(MT, mesh_module_name)
    mesh = Mesh.create_from_remap(vertices, cellnodes, nodes_gid)
    # CHECKME: invert cell_gid, would this be faster in C++?
    gid2lid = {cells_gid[k]: k for k in range(cells_gid.shape[0])}
    local_face_cell = np.array([gid2lid[gci] for gci in face_cell])
    faces_gid_order = mesh.identify_faces_from_positions(local_face_cell, face_position)
    faces_gid = unsorted_faces_gid[faces_gid_order]
    sys.stdout.flush()
    comm.Barrier()
    print('total processing time on', rank, ':', (datetime.datetime.now() - start).total_seconds())
    offsets, cellsnodes = mesh.cells_nodes_as_COC()
    ghost_tag = np.zeros(mesh.nb_cells, dtype='i')
    ghost_tag[ghost_indexes[0]:] = 1
    vtkw.write_vtu(
        vtkw.vtu_doc_from_COC(
            mesh.vertices_array(), 
            np.array(offsets[1:], copy=False), # no first zero offset for vtk 
            np.array(cellsnodes, copy=False),
            mesh.cells_vtk_ids(),
            celldata = {'ghost': ghost_tag}
        ),
        'localmesh-%03d.vtu' % rank
    )
    return mesh

mesh = create_and_distribute_mesh(MT.HexMesh.make, GT.grid2hexs((nx, ny, nz), idtype=MT.idtype()), geometric_coloring)

