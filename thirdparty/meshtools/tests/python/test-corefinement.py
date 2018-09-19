import numpy as np
import MeshTools as MT
import MeshTools.Corefinement as CR

meshes = CR.test()

as_tsurf = lambda t: MT.TSurf.make(t[0], MT.idarray(t[1]))

for mi, mesh in enumerate(meshes):
    print("*** Mesh ", mi)
    print(mesh)
    MT.to_vtu(as_tsurf(mesh), 'tsurf%02d.vtu' % mi, ofmt='ascii')
