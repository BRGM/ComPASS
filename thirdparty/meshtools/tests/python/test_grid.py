import MeshTools as MT

def test_grid():
    shape = nx, ny = 5, 8
    grid = MT.grid3D(shape=shape)
    MT.to_vtu(grid, 'grid.vtu')

if __name__=='__main__':
    test_grid()
