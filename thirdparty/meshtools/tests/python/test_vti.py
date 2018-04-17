# -*- coding: utf-8 -*-

import numpy as np
import MeshTools.vtkwriters as vtkw


def test_vti():

    nx, ny, nz = 10, 5, 8

    xy = np.meshgrid(np.arange(nx),
                     np.arange(ny),
                     indexing='ij')
    x, y = xy

    vtkw.write_vti(
        vtkw.vti_doc(
                    (nx, ny),
                    celldata={'x': x, 'y': y}
                ),
        'test-vti2D'        
    )

    # In 2D must had a fake n+1 dimension
    #vtkw.write_vti(
    #    vtkw.vti_doc(
    #                (nx-1, ny-1),
    #                pointdata={'x': x, 'y': y}
    #            ),
    #    'test-vti2D_pointdata'        
    #)

    xyz = np.meshgrid(np.arange(nx),
                      np.arange(ny),
                      np.arange(nz),
                      indexing='ij')
    x, y, z = xyz

    vtkw.write_vti(
        vtkw.vti_doc(
                    (nx, ny, nz),
                    origin = (x[0,0,0], y[0,0,0], z[0,0,0]),
                    delta = (x[1,0,0]-x[0,0,0],
                             y[0,1,0]-y[0,0,0],
                             z[0,0,1]-z[0,0,0]),
                    celldata={'x': x, 'y': y, 'z': z}
                ),
        'test-vti3D'        
    )

    vtkw.write_vti(
        vtkw.vti_doc(
                    (nx-1, ny-1, nz-1),
                    extent = [(0, n-1) for n in (nx, ny, nz)],
                    pointdata={'x': x, 'y': y, 'z': z}
                ),
        'test-vti3D_pointdata'        
    )

if __name__=='__main__':
    test_vti()
