# -*- coding: utf-8 -*-
#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#


import numpy as np
import ComPASS
from ComPASS.utils.units import *
from ComPASS.timeloops import standard_loop
import thirdparties
import MeshTools as MT

H = 3 * km

ComPASS.load_eos('water2ph')
ComPASS.set_output_directory_and_logfile(__file__)

basename = '47K/input.2'

# extract information from TetGen file
nodefile = basename + '.node'
elementfile = basename + '.ele'
facefile = basename + '.face'
vertices = np.loadtxt(nodefile, skiprows=1, usecols=(1, 2, 3))
tets = np.loadtxt(elementfile, skiprows=1, usecols=(1, 2, 3, 4), dtype=int)
tets-= 1 # start indexing at 0
tets = MT.idarray(tets)
faces = np.loadtxt(facefile, skiprows=1, usecols=(1, 2, 3), dtype=int)
faces-= 1 # start indexing at 0
faces = MT.idarray(faces)

# filter out box boundaries (unit box)
not_boundary = np.ones(faces.shape[0], dtype=bool)
for iaxis in range(3):
    coordi = vertices[:, iaxis][faces]
    for value in (0, 1):
        on_boundary = np.all(coordi==value, axis=1)
        not_boundary[on_boundary] = False
faces = faces[not_boundary]

# scale domain
vertices*= H
mesh = MT.TetMesh.make(vertices, tets)
mesh_faces = mesh.connectivity.faces
# should be available through C++
fracture_faces = [mesh_faces.id(MT.Triangle(triangle)) for triangle in faces]
fracture_faces = np.array(fracture_faces)


ComPASS.init(
    mesh = mesh,
    fracture_faces = lambda: fracture_faces,
)

#print('Old computation')
#print(ComPASS.old_compute_face_centers())
#print('New computation')
#print(ComPASS.compute_face_centers())
assert np.all((len(ComPASS.old_compute_face_centers())==0 and len(ComPASS.compute_face_centers())==0) or
               ComPASS.old_compute_face_centers()==ComPASS.compute_face_centers())
#print('Fractures old computation')
#print(ComPASS.old_compute_fracture_centers())
#print('Fractures new computation')
#print(ComPASS.compute_fracture_centers())
assert np.all((len(ComPASS.old_compute_fracture_centers())==0 and len(ComPASS.compute_fracture_centers())==0) or
               ComPASS.old_compute_face_centers()==ComPASS.compute_face_centers())

