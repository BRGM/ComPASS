from itertools import product
import numpy as np


def chessboard(grid_shape, macro_grid, nbcols, loops_order):
    nx, ny, nz = grid_shape
    colors = np.zeros((nx, ny, nz), dtype=np.int32)
    npx, npy, npz = macro_grid
    sx, sy, sz = nx // npx, ny // npy, nz // npz
    p = 0
    assert nbcols > 0
    loops_argsort = np.argsort(loops_order)
    assert np.all(np.array(loops_order)[loops_argsort] == np.arange(0, 3))
    corners = [
        list(range(0, nx, sx)) + [nx,],
        list(range(0, ny, sy)) + [ny,],
        list(range(0, nz, sz)) + [nz,],
    ]
    slices = [[slice(l[i - 1], l[i]) for i in range(1, len(l))] for l in corners]
    for subblock in product(*[slices[k] for k in loops_order]):
        slx, sly, slz = (subblock[k] for k in loops_argsort)
        colors[slx, sly, slz] = p % nbcols
        p += 1
    return colors


if __name__ == "__main__":
    from MeshTools import vtkwriters as vtkw

    cb = chessboard(
        (30, 20, 6), (3, 2, 3), 2, loops_order=(1, 2, 0,)  # x first then z then y
    )
#    vtkw.write_vti(
#        vtkw.block_as_vti_doc(
#            cb, location='cell', name='color'
#        ), 'chessboard'
#    )
