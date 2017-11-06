

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


#vertices = MT.as_coordinate_array(nodes)
#cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

hexaedra = np.array([np.array(elt) 
            for elt, tags in elements
            if type(elt) is MT.Hexahedron])
mesh = MT.HexMesh.make(nodes, hexaedra)

vertices = MT.as_coordinate_array(mesh.vertices)
cellnodes = np.array([np.array(nodes) for nodes in mesh.connectivity.cells.nodes])

vtkw.write_vtu(
    vtkw.vtu_doc(vertices, cellnodes),
    filename.replace('.msh', '.vtu')
)

#vtk_ids = mesh.vtk_ids()

