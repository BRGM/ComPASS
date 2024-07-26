import numpy as np
import meshio
from ComPASS.messages import error
from ComPASS import RawMesh


def vtk2RawMesh(
    vtk_file,
):
    vtk_mesh = meshio.read(vtk_file)
    # get mesh connectivity
    faces_nodes = []
    cells_faces = []
    cells_nodes = []
    cells_types = []
    fi = 0
    for vtktype, c_nodes in vtk_mesh.cells_dict.items():
        nb_cells = len(c_nodes)
        # tetrahedral cells
        if vtktype == "tet" or vtktype == "tetra":
            c_types = np.full((nb_cells,), 10, dtype=int)
            cells_types.extend(c_types)
            cells_nodes.extend(c_nodes)
            for nodes in c_nodes:
                faces_nodes.append([nodes[0], nodes[1], nodes[2]])
                faces_nodes.append([nodes[0], nodes[1], nodes[3]])
                faces_nodes.append([nodes[1], nodes[2], nodes[3]])
                faces_nodes.append([nodes[2], nodes[0], nodes[3]])
                cells_faces.append(np.arange(fi, fi + 4))
                fi += 4

        elif vtktype == "voxel":
            c_types = np.full((nb_cells,), 11, dtype=int)
            cells_types.extend(c_types)
            cells_nodes.extend(c_nodes)
            for nodes in c_nodes:
                faces_nodes.append([nodes[0], nodes[1], nodes[3], nodes[2]])
                faces_nodes.append([nodes[4], nodes[5], nodes[7], nodes[6]])
                faces_nodes.append([nodes[0], nodes[1], nodes[5], nodes[4]])
                faces_nodes.append([nodes[1], nodes[3], nodes[7], nodes[5]])
                faces_nodes.append([nodes[3], nodes[2], nodes[6], nodes[7]])
                faces_nodes.append([nodes[2], nodes[0], nodes[4], nodes[6]])
                cells_faces.append(np.arange(fi, fi + 6))
                fi += 6

        # wedge cells
        elif vtktype == "wedge":
            c_types = np.full((nb_cells,), 13, dtype=int)
            cells_types.extend(c_types)
            cells_nodes.extend(c_nodes)
            for nodes in c_nodes:
                faces_nodes.append(nodes[:3])
                faces_nodes.append(nodes[3:])
                faces_nodes.append([nodes[0], nodes[1], nodes[4], nodes[3]])
                faces_nodes.append([nodes[0], nodes[3], nodes[5], nodes[2]])
                faces_nodes.append([nodes[1], nodes[4], nodes[5], nodes[2]])
                cells_faces.append(np.arange(fi, fi + 5))
                fi += 5

        # check the pyramid keyword !
        elif vtktype == "pyramid":
            c_types = np.full((nb_cells,), 14, dtype=int)
            cells_types.extend(c_types)
            cells_nodes.extend(c_nodes)
            for nodes in c_nodes:
                faces_nodes.append(nodes[:4])
                faces_nodes.append([nodes[0], nodes[1], nodes[4]])
                faces_nodes.append([nodes[1], nodes[2], nodes[4]])
                faces_nodes.append([nodes[2], nodes[3], nodes[4]])
                faces_nodes.append([nodes[3], nodes[0], nodes[4]])
                cells_faces.append(np.arange(fi, fi + 5))
                fi += 5

        else:
            error("cell type not recognized in vtk2RawMesh")

    # get mesh data
    # cell data
    cells_data_dict = {}
    for data_name, data_dict in vtk_mesh.cell_data_dict.items():
        cells_data = []
        for data_celltype, data in data_dict.items():
            # "if" are necessary and in same order as previously to keep the cells order
            if data_celltype == "tet" or data_celltype == "tetra":
                cells_data.extend(data)
            elif data_celltype == "voxel":
                cells_data.extend(data)
            elif data_celltype == "wedge":
                cells_data.extend(data)
            elif data_celltype == "pyramid":
                cells_data.extend(data)

        cells_data_dict[data_name] = np.array(cells_data)

    return RawMesh(
        vertices=np.asarray(vtk_mesh.points),
        face_nodes=faces_nodes,
        # same nb of nodes in each cell
        cell_nodes=np.asarray(cells_nodes),
        # same nb of faces in each cell
        cell_faces=np.asarray(cells_faces),
        cell_types=np.asarray(cells_types),
        **cells_data_dict,
    )
