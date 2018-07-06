import glob, re, os
import numpy as np
import MeshTools.vtkwriters as vtkw
import ComPASS.mpi as mpi

def extract_mesh(basename):
    info = []
    with open('nodeinfo' + basename) as f:
        for line in f:
            info.append(line.strip().split(' '))
    info = np.array(info, dtype='c')
    pointdata = {
        'own': np.array(info[:,0]==b'o', dtype=np.int),
        'frac': np.array(info[:,1]==b'y', dtype=np.int),
        'dirichlet_P': np.array(info[:,2]==b'd', dtype=np.int),
        'dirichlet_T': np.array(info[:,3]==b'd', dtype=np.int),
    }
    facedata = {}
    celldata = {}
    def collect(property, dictionnary, dtype):
        filename = property + basename
        if os.path.exists(filename):
            dictionnary[property] = np.loadtxt(filename, dtype=dtype)
            #print('loaded', property, dictionnary[property].shape, 'for basename', basename) 
    collect('nodeflags', pointdata, np.int)
    collect('cellflags', celldata, np.int)
    collect('facecenters', facedata, np.double)
    collect('faceflags', facedata, np.int)
    vertices = np.loadtxt('nodes' + basename)
    #print('loaded', vertices.shape, 'vertices for basename', basename) 
    cells = np.loadtxt('cells' + basename, dtype=np.uint64)
    #print('loaded', cells.shape, 'cells for basename', basename)     
    fractures = None
    fracturefile = 'fractures' + basename
    if os.path.exists(fracturefile):
        has_fractures = False
        with open(fracturefile) as f:
            if f.readline().strip():
                has_fractures = True
        if has_fractures:
            fractures = np.loadtxt(fracturefile, dtype=np.uint64)
    return vertices, cells, pointdata, facedata, celldata, fractures

def dump_mesh(basename, vertices, cells, pointdata, facedata, celldata, fractures=None):
    assert all(np.asarray(a).shape == (vertices.shape[0],) for a in pointdata.values())
    assert all(np.asarray(a).shape == (cells.shape[0],) for a in celldata.values())
    vtkw.write_vtu(
        vtkw.vtu_doc(
            vertices, cells, 
            pointdata = pointdata, celldata=celldata, ofmt='ascii'
        ),
        'mesh' + basename,
    ) 
    if facedata:
        assert 'facecenters' in facedata
        face_centers = facedata['facecenters']
        assert all(np.asarray(a).shape == (face_centers.shape[0],) for key, a in facedata.items() if key is not 'facecenters')
        vtkw.write_vtu(
            vtkw.points_as_vtu_doc(
                face_centers, 
                pointdata = {name: value for name, value in facedata.items() if name!='facecenters'}
            ),
            'facedata' + basename + '.vtu',
        ) 
    if fractures is not None:
        vtkw.write_vtu(
            vtkw.vtu_doc(vertices, fractures),
            'frac_mesh' + basename + '.vtu',
        ) 

@mpi.on_master_proc
def extract_all_meshes():
    files = glob.glob('nodes[0-9]*')
    extensions = [s[-5:] for s in files]
    for basename in extensions:
        data = extract_mesh(basename)
        dump_mesh(basename, *data)

if __name__=='__main__':
    extract_all_meshes()

