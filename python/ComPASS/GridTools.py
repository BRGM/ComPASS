import numpy as np

def grid2tets(shape, extent=(1., 1., 1.)):
    # number of cells
    ncx, ncy, ncz = shape
    # number of nodes
    nx, ny, nz = ncx+1, ncy+1, ncz+1
    dx, dy, dz = extent
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
