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

segments = np.array([
    (0, 5),
    (1, 6),
    (3, 8),
    (0, 5),
    (2, 7),
    (4, 8),
], dtype=np.int32)

tvertices, triangles, component, faces = PM.mesh(vertices, segments)

assert np.all(tvertices[:vertices.shape[0]]==vertices)

plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(tvertices[:,0], tvertices[:,1], triangles, color='grey')
for S in segments:
    plt.plot(vertices[S, 0], vertices[S, 1], 'r') 
tcenters = np.vstack([tvertices[triangle].mean(axis=0) for triangle in triangles])
for i, xy in enumerate(tcenters):
    x, y = xy
    plt.text(x, y, '%d' % component[i])
for i, xy in enumerate(tvertices):
    x, y = xy
    plt.text(x, y, '%d' % i, color='b')
for fi, face in enumerate(faces):
    print("face", fi, ":", face)
plt.show()
