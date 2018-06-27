import numpy as np
import MeshTools.PetrelMesh as PM

import matplotlib.pyplot as plt
import matplotlib.tri as tri

vertices = np.array([
    (0, 0),
    (0, 0.2),
    (0, 0.4),
    (0, 0.6),
    (0, 1),
    (1, 0),
    (1, 0.6),
    (1, 0.2),
    (1, 1.),
], dtype=np.double)

S1 = np.array([
    (0, 5),
    (1, 6),
    (3, 8),
], dtype=np.int32)

S2 = np.array([
    (0, 5),
    (2, 7),
    (4, 8),
], dtype=np.int32)

tvertices, triangles, component = PM.split_and_mesh(vertices, S1, S2)

plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(tvertices[:,0], tvertices[:,1], triangles)
for S in S1:
    plt.plot(vertices[S, 0], vertices[S, 1], 'r') 
for S in S2:
    plt.plot(vertices[S, 0], vertices[S, 1], 'k') 
tcenters = np.vstack([tvertices[triangle].mean(axis=0) for triangle in triangles])
for i, xy in enumerate(tcenters):
    x, y = xy
    plt.text(x, y, '%d' % component[i])
plt.show()



