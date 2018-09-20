import numpy as np
from skimage.measure import marching_cubes_lewiner as marching_cubes

import MeshTools as MT
import MeshTools.Corefinement as CR


shape = nx, ny, nz = (10,)*3
steps = tuple(np.linspace(-1, 1, n) for n in shape)
coordinates = np.meshgrid(*steps, indexing='ij')
#points = np.array(point_arrays)
#points.shape = (-1, 3) 
points = np.stack(coordinates, axis=-1)
points.shape = (-1, 3)

image1 = CR.f1(points)
image1.shape = shape
tsurf1 = marching_cubes(image1, level=[rank_values])


image2 = CR.f2(points)
image2.shape = shape
tsurf2 = marching_cubes(image2)

#meshes = CR.test()

as_tsurf = lambda t: MT.TSurf.make(t[0], MT.idarray(t[1]))

for mi, triangles in enumerate((tsurf1, tsurf2)):
    print("*** MC ", mi)
    verts, faces, normals, values = triangles
    #print(mesh)
    MT.to_vtu(as_tsurf((verts, faces)), 'mc_tsurf%02d.vtu' % mi, ofmt='ascii')

#for mi, mesh in enumerate(meshes):
#    print("*** Mesh ", mi)
#    #print(mesh)
#    MT.to_vtu(as_tsurf(mesh), 'tsurf%02d.vtu' % mi, ofmt='ascii')
