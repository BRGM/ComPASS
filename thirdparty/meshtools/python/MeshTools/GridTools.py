import numpy as np

class GridInfo:
    def __init__(self, shape=None, extent=None, origin=None):
        assert shape is not None or extent is not None
        if shape is None:
            shape = (1,) * len(tuple(extent))
        shape = tuple(shape)
        dim = len(shape)
        assert dim<=3
        if origin is None:
            origin = (0.,) * dim
        origin = tuple(origin)
        assert len(origin)==dim
        if extent is None:
            extent = (1.,) * dim
        extent = tuple(extent)
        assert len(extent)==dim
        if dim<3:
            shape+= (1,) * (3 - dim) # grid is always in 3D
            origin+= (0.,) * (3 - dim) # grid is always in 3D
            extent+= (1.,) * (3 - dim) # grid is always in 3D
        self.shape = shape
        self.extent = extent
        self.origin = origin

def grid_elements(node_shape=None, extent=None, cell_shape=None, steps=None,
                  quad=False):
    if not cell_shape is None:
        assert node_shape is None
        node_shape = tuple(n+1 for n in cell_shape)
    if not extent is None:
        assert steps is None
        steps = (L/(n-1) for L, n in zip(extent, node_shape))
    nx, ny, nz = node_shape    
    dx, dy, dz = steps
    x = np.arange(0, nx * dx - 0.5*dx, dx)
    y = np.arange(0, ny * dy - 0.5*dy, dy)
    z = np.arange(0, nz * dz - 0.5*dz, dz)    
    vertices = np.empty((nx*ny*nz, 3), dtype='d')
    vertices[:, 0] = np.tile(x, ny*nz)
    vertices[:, 1] = np.tile(np.hstack([np.tile(yi, nx) for yi in y]), nz)
    vertices[:, 2] = np.hstack([np.tile(zi, nx*ny) for zi in z])    
    ncx, ncy, ncz = nx-1, ny-1, nz-1
    cells = np.empty((ncx*ncy*ncz, 8), dtype='i')
    if quad:
        cells[0] = (0, 1, nx, nx+1, nx*ny, nx*ny+1, nx*ny+nx, nx*ny+nx+1)
    else: # hex
        cells[0] = (0, 1, nx+1, nx, nx*ny, nx*ny+1, nx*ny+nx+1, nx*ny+nx)
    for i in range(1, ncx):
        cells[i] = cells[i-1] + 1
    for j in range(1, ncy):
        cells[j*ncx:(j+1)*ncx] = cells[(j-1)*ncx:j*ncx] + nx
    for k in range(1, ncz):
        cells[k*(ncx*ncy):(k+1)*(ncx*ncy)] = cells[(k-1)*(ncx*ncy):k*(ncx*ncy)] + nx*ny
    return vertices, cells

def grid2tets(shape, extent=(1., 1., 1.)):
    # number of cells
    ncx, ncy, ncz = shape
    # number of nodes
    nx, ny, nz = ncx+1, ncy+1, ncz+1
    dxyz = np.array([L/d for L, d in zip(extent, shape)], dtype=np.double)
    assert np.all(dxyz>0)
    dx, dy, dz = dxyz
    ncubes = (nx-1)*(ny-1)*(nz-1)
    nnodes = nx * ny * nz
    ncenters = (nx-1) * (ny-1) * (nz-1)
    nfacesx = nx * (ny-1) * (nz-1) 
    nfacesy = ny * (nx-1) * (nz-1) 
    nfacesz = nz * (nx-1) * (ny-1) 
    nfaces = nfacesx + nfacesy + nfacesz

    vertices = np.zeros((nnodes + ncenters + nfaces, 3), dtype=np.double)

    offc = nnodes
    offfx = offc + ncenters
    offfy = offfx + nfacesx
    offfz = offfy + nfacesy
    corner = lambda i, j, k: i*ny*nz + j*nz + k
    center = lambda i, j, k: offc + i*(ny-1)*(nz-1) + j*(nz-1) + k
    fcx = lambda i, j, k: offfx + i*(ny-1)*(nz-1) + j*(nz-1) + k
    fcy = lambda i, j, k: offfy + i*(ny)*(nz-1) + j*(nz-1) + k
    fcz = lambda i, j, k: offfz + i*(ny-1)*(nz) + j*(nz) + k

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
               vertices[corner(i,j,k)] = (i*dx, j*dy, k*dz)
    for i in range(nx-1):
        for j in range(ny-1):
            for k in range(nz-1):
               vertices[center(i,j,k)] = ((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz)
    for i in range(nx):
        for j in range(ny-1):
            for k in range(nz-1):
               vertices[fcx(i,j,k)] = (i*dx, (j+0.5)*dy, (k+0.5)*dz)
    for i in range(nx-1):
        for j in range(ny):
            for k in range(nz-1):
               vertices[fcy(i,j,k)] = ((i+0.5)*dx, j*dy, (k+0.5)*dz)
    for i in range(nx-1):
        for j in range(ny-1):
            for k in range(nz):
               vertices[fcz(i,j,k)] = ((i+0.5)*dx, (j+0.5)*dy, k*dz)

    tets = np.zeros((24*ncubes, 4), dtype=np.int)

    tet = 0
    for i in range(nx-1):
        for j in range(ny-1):
            for k in range(nz-1):
                # faces z- bottom
                tets[tet] = (corner(i,j,k), corner(i+1,j,k), fcz(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j,k), corner(i+1,j+1,k), fcz(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k), corner(i,j+1,k), fcz(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j+1,k), corner(i,j,k), fcz(i,j,k), center(i,j,k))
                tet+=1
                # faces z+ top
                tets[tet] = (corner(i,j,k+1), corner(i+1,j,k+1), fcz(i,j,k+1), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j,k+1), corner(i+1,j+1,k+1), fcz(i,j,k+1), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k+1), corner(i,j+1,k+1), fcz(i,j,k+1), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j+1,k+1), corner(i,j,k+1), fcz(i,j,k+1), center(i,j,k))
                tet+=1
                # faces x-
                tets[tet] = (corner(i,j,k), corner(i,j+1,k), fcx(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j+1,k), corner(i,j+1,k+1), fcx(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j+1,k+1), corner(i,j,k+1), fcx(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j,k+1), corner(i,j,k), fcx(i,j,k), center(i,j,k))
                tet+=1
                # faces x+
                tets[tet] = (corner(i+1,j,k), corner(i+1,j+1,k), fcx(i+1,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k), corner(i+1,j+1,k+1), fcx(i+1,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k+1), corner(i+1,j,k+1), fcx(i+1,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j,k+1), corner(i+1,j,k), fcx(i+1,j,k), center(i,j,k))
                tet+=1
                # faces y-
                tets[tet] = (corner(i,j,k), corner(i+1,j,k), fcy(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j,k), corner(i+1,j,k+1), fcy(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j,k+1), corner(i,j,k+1), fcy(i,j,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j,k+1), corner(i,j,k), fcy(i,j,k), center(i,j,k))
                tet+=1
                # faces y+
                tets[tet] = (corner(i,j+1,k), corner(i+1,j+1,k), fcy(i,j+1,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k), corner(i+1,j+1,k+1), fcy(i,j+1,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i+1,j+1,k+1), corner(i,j+1,k+1), fcy(i,j+1,k), center(i,j,k))
                tet+=1
                tets[tet] = (corner(i,j+1,k+1), corner(i,j+1,k), fcy(i,j+1,k), center(i,j,k))
                tet+=1
    
    return vertices, tets

def steps2hex(steps_along_axes, idtype=np.int):
    x, y, z = (np.asarray(steps) for steps in steps_along_axes)
    assert all(len(a.shape)==1 for a in (x, y, z))
    # number of nodes
    shape = nx, ny, nz = (steps.shape[0] for steps in (x, y, z))
    assert all(n>1 for n in shape)
    nnodes = nx * ny * nz
    # number of cells
    ncx, ncy, ncz = nx-1, ny-1, nz-1
    origin = np.array((x[0], y[0], z[0]))
    assert all(np.all(a[:-1]<a[1:]) for a in (x, y, z))
    vertices = np.zeros((nnodes, 3), dtype=np.double)
    vertices[:, 0] = np.tile(x, ny*nz)
    vertices[:, 1] = np.tile(np.hstack([np.tile(yi, nx) for yi in y]), nz)
    vertices[:, 2] = np.hstack([np.tile(zi, nx*ny) for zi in z])
    nhexs = ncx * ncy * ncz
    hexs = np.zeros((nhexs, 8), dtype=idtype)
    tmp = np.arange(ncx)
    hexs[:ncx, 0] = tmp
    hexs[:ncx, 1] = tmp + 1  
    hexs[:ncx, 2] = tmp + 1 + nx  
    hexs[:ncx, 3] = tmp + nx  
    for j in range(1, ncy):
        hexs[j*ncx:(j+1)*ncx, :4] = hexs[(j-1)*ncx:j*ncx, :4] + nx
    hexs[:(ncx*ncy), 4:] = hexs[:(ncx*ncy), :4] + nx*ny
    for k in range(1, ncz):
        hexs[k*(ncx*ncy):(k+1)*(ncx*ncy), :] = hexs[(k-1)*(ncx*ncy):k*(ncx*ncy), :] + nx*ny
    return vertices, hexs

def grid2hexs(idtype=np.int, **kwargs):
    if 'gridinfo' in kwargs:
        gridinfo = kwargs['gridinfo']
    else:
        gridinfo = GridInfo(**kwargs)
    return steps2hex([
        np.linspace(O, O + L, n + 1) # nb cells -> nb nodes
        for O, L, n in zip(gridinfo.origin, gridinfo.extent, gridinfo.shape)
    ], idtype=idtype)
