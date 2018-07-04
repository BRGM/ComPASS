import numpy as np
import MeshTools as MT
import MeshTools.Corefinement as CR

tm1, tm2 = CR.test()
print(tm1)
print(tm2)

as_tsurf = lambda t: MT.TSurf.make(t[0], MT.idarray(t[1]))
m1 = as_tsurf(tm1)
MT.to_vtu(m1, 'tsurf1.vtu', ofmt='ascii')
MT.to_vtu(as_tsurf(tm2), 'tsurf2.vtu')
