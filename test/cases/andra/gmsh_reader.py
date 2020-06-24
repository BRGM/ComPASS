#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


# cf. format description at:
# http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format

import numpy as np
import MeshTools as MT
import MeshTools.vtkwriters as vtkw

# conversion GMesh element code -> MeshTools object
# cf codelist at: http://gmsh.info/doc/texinfo/gmsh.html#MSH-ASCII-file-format
element_factory = {
    2: MT.Triangle,
    3: MT.Quad,
    4: MT.Tetrahedron,
    5: MT.Hexahedron,
    6: MT.Wedge,
}


def retrieve_nodes(stream):
    result = []
    line = stream.readline().strip()
    nbnodes = int(line)
    for i in range(nbnodes):
        line = stream.readline().strip()
        assert line
        line = line.split()
        assert int(line[0]) == i + 1
        result.append(tuple(float(s) for s in line[1:]))
    line = stream.readline().strip()
    assert line.startswith("$EndNodes")
    return np.array(result, dtype=np.double)


def retrieve_elements(stream):
    result = []
    line = stream.readline().strip()
    nbelts = int(line)
    for i in range(nbelts):
        line = stream.readline().strip()
        assert line
        line = tuple(int(s) for s in line.split())
        assert line[0] == i + 1
        classid = line[1]
        nbtags = line[2]
        tags = line[3 : 3 + nbtags]
        nodes = tuple(i - 1 for i in line[3 + nbtags :])  # node indexing starts at 1
        result.append((element_factory[classid](nodes), tags))
    line = stream.readline().strip()
    assert line.startswith("$EndElements")
    return result


def retrieve_mesh_elements(filename):
    with open(filename) as f:
        line = f.readline().strip()
        while line:
            if line.startswith("$Nodes"):
                nodes = retrieve_nodes(f)
            if line.startswith("$Elements"):
                elements = retrieve_elements(f)
            line = f.readline().strip()
    return nodes, elements
