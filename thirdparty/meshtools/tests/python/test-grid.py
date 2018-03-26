import MeshTools as MT

shape = nx, ny = 5, 8

grid = MT.grid3D(shape=shape)

MT.to_vtu(grid, 'grid.vtu')
