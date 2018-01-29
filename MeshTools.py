import numpy as np
from _MeshTools import *

idarray = lambda a: np.asarray(a, dtype=idtype())

## Tet Volumes
#A = mesh.vertices[mesh.cellnodes[:, 0]]
#AB, AC, AD = (mesh.vertices[mesh.cellnodes[:, i+1]] - A for i in range(3))
#det = AB[:, 0] * AC[:, 1] * AD[:, 2]
#det+= AB[:, 1] * AC[:, 2] * AD[:, 0]
#det+= AB[:, 2] * AC[:, 0] * AD[:, 1]
#det-= AB[:, 2] * AC[:, 1] * AD[:, 0]
#det-= AB[:, 0] * AC[:, 2] * AD[:, 1]
#det-= AB[:, 1] * AC[:, 0] * AD[:, 2]
#vol = np.abs(det) / 6.

#print(vol.min(), vol.max())
