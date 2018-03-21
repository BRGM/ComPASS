import numpy as np
import GridTools as GT
import MeshTools as MT
import vtkwriters as vtkw

nx, ny, nz = 4, 4, 1
nbparts = 4

mesh = MT.HexMesh.make(*GT.grid2hexs((nx, ny, nz), idtype=MT.idtype()))
cell_color = np.zeros(mesh.nb_cells, dtype='i')
for i in range(1, nbparts):
    cell_color[i*(mesh.nb_cells//nbparts):]+= 1
distribution = mesh.distribute(nbparts, cell_color)

node_color = np.empty(mesh.nb_vertices, dtype='i')
face_color = np.empty(mesh.nb_faces, dtype='i')
for proc, dist in enumerate(distribution):
    cells, nodes, faces = [v.array_view() for v in dist[0]]
    ghost_cells_index, ghost_nodes_index, ghost_faces_index = dist[1]
    node_color[nodes[:ghost_nodes_index]] = proc 
    face_color[faces[:ghost_faces_index]] = proc

offsets, cellsnodes = mesh.cells_nodes_as_COC()
vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        mesh.vertices_array(), 
        np.array(offsets[1:], copy=False), # no first zero offset for wtk 
        np.array(cellsnodes, copy=False),
        mesh.cells_vtk_ids(),
        pointdata={'color': node_color},
        celldata={'color': cell_color}
    ),
    'distribution.vtu'
)
    
offsets, facesnodes = mesh.faces_nodes_as_COC()
vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        mesh.vertices_array(), 
        np.array(offsets[1:], copy=False), # no first zero offset for wtk 
        np.array(facesnodes, copy=False),
        mesh.faces_vtk_ids(),
        pointdata={'color': node_color},
        celldata={'color': face_color}
    ),
    'distribution_faces.vtu'
)
