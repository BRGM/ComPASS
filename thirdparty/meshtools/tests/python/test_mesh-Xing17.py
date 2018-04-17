# -*- coding: utf-8 -*-

import numpy as np
import MeshTools as MT

def test_mesh_Xing17():

    H = 3E3

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

    MT.to_vtu(mesh, 'Xing17.vtu')

    tsurf = MT.TSurf.create(vertices, [mesh.connectivity.faces.nodes[fk] for fk in fracture_faces])

    MT.to_vtu(tsurf, 'Xing17-fractures.vtu', celldata={'face': fracture_faces})

if __name__=='__main__':
    test_mesh_Xing17()
