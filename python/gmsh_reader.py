

# cf. format description at:
# http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format

import numpy as np
import thirdparties
import MeshTools as MT
import vtkwriters as vtkw

# conversion GMesh element code -> MeshTools object
# cf codelist at: http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
element_factory = {
        2: MT.Triangle,
        3: MT.Quad,
        4: MT.Tetrahedron,
        5: MT.Hexahedron,
        6: MT.Wedge,
        }

filename = 'case-tests/andra/andra_gallery.msh'

def retrieve_nodes(stream):
    result = []
    line = stream.readline().strip()
    nbnodes = int(line)
    for i in range(nbnodes):
        line = stream.readline().strip()    
        assert line
        line = line.split()
        assert int(line[0]) == i+1
        result.append(tuple(float(s) for s in line[1:]))
    line = stream.readline().strip()
    assert line.startswith('$EndNodes')
    return np.array(result, dtype=np.double)

def retrieve_elements(stream):
    result = []
    line = stream.readline().strip()
    nbelts = int(line)
    for i in range(nbelts):
        line = stream.readline().strip()    
        assert line
        line = tuple(int(s) for s in line.split())
        assert line[0] == i+1
        classid = line[1]
        nbtags = line[2]
        tags = line[3:3+nbtags]
        nodes = tuple(i-1 for i in line[3+nbtags:]) # node indexing starts at 1
        result.append((element_factory[classid](nodes), tags))
    line = stream.readline().strip()
    assert line.startswith('$EndElements')
    return result

with open(filename) as f:
    line = f.readline().strip()
    while line:
        if line.startswith('$Nodes'):
            nodes = retrieve_nodes(f)
        if line.startswith('$Elements'):
            elements = retrieve_elements(f)
        line = f.readline().strip()

# select heaxedra only
hexaedra = [elt for elt, tags in elements if type(elt) is MT.Hexahedron]
mesh = MT.HexMesh.create(nodes, hexaedra)

vertices = MT.as_coordinate_array(mesh.vertices)
cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, cellnodes),
    filename.replace('.msh', '_hexaedra_only.vtu')
)

# filter out 2D elements
elements_3D = [(elt, tags) for elt, tags in elements
               if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)]

mesh = MT.HybridMesh.create(nodes, [elt for elt, tags in elements_3D])
physical = np.array([tags[0] for elt, tags in elements_3D])

offsets, cellnodes = mesh.cells_nodes_as_COC()

vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        nodes, 
        np.array(offsets[1:], copy=False), # no first zero offset for wtk 
        np.array(cellnodes, copy=False),
        mesh.vtk_ids(),
        celldata={'physical': physical}
    ),
    filename.replace('.msh', '_3D.vtu')
)
    
# Indentify rear faces and export them as 2D mesh with the physical tag
fc = MT.as_coordinate_array(mesh.face_centers())
bf = MT.as_id_array(mesh.boundary_faces())
where = fc[bf, 0]==0 # x coordinate == 0 ?
rear_faces = bf[where]
# retrieve rear_faces tag using their cells
bc = MT.as_id_array(mesh.boundary_cells())
rear_faces_tag = physical[bc[where]]

face_nodes = mesh.connectivity.faces.nodes
rear_elements = [face_nodes[fk] for fk in rear_faces]
rear_celltypes = np.array([elt.vtk_id() for elt in rear_elements])
rear_cellnodes = [np.array(elt) for elt in rear_elements]
rear_offsets = np.cumsum([a.shape[0] for a in rear_cellnodes])
rear_cellnodes = np.hstack(rear_cellnodes)

vtkw.write_vtu(
    vtkw.vtu_doc_from_COC(
        nodes, 
        rear_offsets,
        rear_cellnodes,
        rear_celltypes,
        celldata={'physical': rear_faces_tag}
    ),
    filename.replace('.msh', '_rear_face.vtu')
)
