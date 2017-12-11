# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import os
import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import thirdparties
import MeshTools as MT
import vtkwriters as vtkw
import gmsh_reader

ComPASS.load_eos('water2ph')

filename = 'case-tests/andra/andra_gallery.msh'

nodes, elements = gmsh_reader.retrieve_mesh_elements(filename)

# filter out 2D elements
elements = [(elt, tags) for elt, tags in elements
            if type(elt) in (MT.Tetrahedron, MT.Wedge, MT.Hexahedron)]
tags = [tags for elt, tags in elements]
physical = np.array([tag[0] for tag in tags])
elements = [elt for elt, tags in elements]

mesh = MT.HybridMesh.create(nodes, elements)

def set_physical_flags():
    flags = ComPASS.global_cellflags()
    flags[:] = physical

ComPASS.set_output_directory_and_logfile(__file__)

ComPASS.init(
    mesh = mesh,
    set_global_flags = set_physical_flags
)
#set_initial_values()
#set_boundary_conditions()
#
#standard_loop(final_time = 30 * year, output_period = year)

ComPASS.dumps.dump_own_element_numbers()
ComPASS.dumps.dump_mesh()
ComPASS.dumps.dump_states()
