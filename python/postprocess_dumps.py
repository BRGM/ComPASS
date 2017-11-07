# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import glob, os
import numpy as np
import thirdparties
import vtkwriters as vtkw

directory = 'output-test-simulation_on_gmsh_mesh'
assert os.path.isdir(directory)
todir = lambda s : os.path.join(directory, s)

for meshfile in glob.glob(todir('mesh_proc_*.npz')):
    head, tail = os.path.split(meshfile)
    proc = tail[5:]
    datafile = todir('state_' + proc)
    assert os.path.exists(datafile)
    print('processing *_' + proc)
    mesh = np.load(meshfile)
    data = np.load(datafile)
    extract_data = lambda s: { name.replace(s, ''): data[name] for name in data.files if name.startswith(s) }
    vtkw.write_vtu(
        vtkw.vtu_doc_from_COC(
            mesh['vertices'], 
            mesh['cellnodes_offsets'], 
            mesh['cellnodes_values'], 
            mesh['celltypes'], 
            pointdata = extract_data('node_'),
            celldata = extract_data('cell_'),
        ),
        todir(proc.replace('.npz', '.vtu'))
    )

