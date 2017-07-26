import numpy as np

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

def grid2hexs(shape, extent=(1., 1., 1.)):
    # number of cells
    ncx, ncy, ncz = shape
    # number of nodes
    nx, ny, nz = ncx+1, ncy+1, ncz+1
    x, y, z = [np.linspace(0, L, n) for L, n in zip(extent, (nx, ny, nz))]
    nnodes = nx * ny * nz
    vertices = np.zeros((nnodes, 3), dtype=np.double)
    vertices[:, 0] = np.tile(x, ny*nz)
    vertices[:, 1] = np.tile(np.hstack([np.tile(yi, nx) for yi in y]), nz)
    vertices[:, 2] = np.hstack([np.tile(zi, nx*ny) for zi in z])
    nhexs = ncx * ncy * ncz
    hexs = np.zeros((nhexs, 8), dtype=np.int)
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