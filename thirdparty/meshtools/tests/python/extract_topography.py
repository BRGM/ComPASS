import numpy as np
import MeshTools as MT
import MeshTools.CGALWrappers as CGAL

c3t3 = CGAL.C3t3('J/J.c3t3.binary.cgal')
topo_facets = np.where(c3t3.facet_tags[:, 0]==-1)[0]

topo = MT.TSurf.make(c3t3.vertices, MT.idarray(c3t3.facets[topo_facets]))
MT.to_vtu(topo, 'topo.vtu')

dtm = CGAL.dtm_from_triangles(c3t3.vertices, c3t3.facets[topo_facets])
depth = dtm.depths(c3t3.vertices)
mesh = MT.TetMesh.make(c3t3.vertices, MT.idarray(c3t3.cells))
MT.to_vtu(mesh, 'mesh.vtu', pointdata={'depth': depth})



