#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

# WARNING: Set your python path adequately
# e.g.: export PYTHONPATH=/home/simon/ComPASS/python
import MeshTools as MT
import vtkwriters as vtkw

# radius
inner_radius = 1.0
external_radius = 100.0
nr = 10
nsectors = 12
thickness = 10.0
nslices = 5

ntheta = nsectors
ny = nslices + 1
nnodes = nr * ntheta * ny

# you may use np.logspace or any other distribution
r = np.linspace(inner_radius, external_radius, nr)
y = np.linspace(0, thickness, ny)

vertices = np.zeros((nnodes, 3), dtype=np.double)
for j, thetaj in enumerate(np.linspace(0, 2 * np.pi, nsectors + 1)[:-1]):
    vertices[j * nr : (j + 1) * nr, 0] = r * np.cos(thetaj)
    vertices[j * nr : (j + 1) * nr, 2] = r * np.sin(thetaj)
    vertices[j * nr : (j + 1) * nr, 1] = y[0]
for k, yk in enumerate(y[1:]):
    vertices[(k + 1) * (nr * ntheta) : (k + 2) * (nr * ntheta)] = vertices[
        k * (nr * ntheta) : (k + 1) * (nr * ntheta)
    ]
    vertices[(k + 1) * (nr * ntheta) : (k + 2) * (nr * ntheta), 1] = yk

ncr, ncth, ncy = nr - 1, nsectors, nslices
ncells = ncr * ncth * ncy
cells = np.zeros((ncells, 8), dtype=MT.idtype())

tmp = np.arange(ncr)
cells[:ncr, 0] = tmp
cells[:ncr, 1] = tmp + 1
cells[:ncr, 2] = tmp + 1 + nr
cells[:ncr, 3] = tmp + nr
for j in range(1, ncth - 1):
    cells[j * ncr : (j + 1) * ncr, :4] = cells[(j - 1) * ncr : j * ncr, :4] + nr
cells[(ncth - 1) * ncr : ncth * ncr, :2] = cells[
    (ncth - 2) * ncr : (ncth - 1) * ncr, 3:1:-1
]
cells[(ncth - 1) * ncr : ncth * ncr, 2:4] = cells[:ncr, 1::-1]
cells[: (ncr * ncth), 4:] = cells[: (ncr * ncth), :4] + nr * ntheta
for k in range(1, ncy):
    cells[k * (ncr * ncth) : (k + 1) * (ncr * ncth), :] = (
        cells[(k - 1) * (ncr * ncth) : k * (ncr * ncth), :] + nr * ntheta
    )

mesh = MT.hexmesh(vertices, cells)

vtkw.write_vtu(vtkw.vtu_doc(mesh.vertices, mesh.cellnodes), "andra.vtu")
