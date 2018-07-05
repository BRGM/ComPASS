import numpy as np
import MeshTools as MT

x = y = np.linspace(0, 1, 10)
z = (np.logspace(0,1,10)-1)/9
mesh = MT.grid3D(steps=(x, y, z))

MT.to_vtu(mesh, 'grid_with_variable_steps.vtu')

H = 10
L = 50
nx = 30
nz = 6
x = np.hstack([
        np.linspace(0, 0.3*L, nx//2)[:-1],
        np.logspace(np.log10(0.3*L), np.log10(L), nx-nx//2)
        ])
y = [-1, 1]
zpos = H * (np.cumprod(np.tile(2, nz))) / 2**nz
z = np.hstack([-zpos[::-1], 0, zpos])
mesh = MT.grid3D(steps=(x, y, z))

MT.to_vtu(mesh, 'grid_with_fracture.vtu')
