import glob, re
import numpy as np
import MeshTools as MT
import MeshTools.vtkwriters as vtkw

def extract_meshes():
    files = glob.glob('nodes[0-9]*')
    extensions = [s[-5:] for s in files]
    for s in extensions:
        info = []
        with open('nodeinfo'+s) as f:
            for line in f:
                info.append(line.strip().split(' '))
        info = np.array(info, dtype='c')
        vertices = np.loadtxt('nodes'+s)
        tets = np.loadtxt('cells'+s, dtype=np.uint64)
        mesh = MT.TetMesh.make(vertices, tets)
        MT.to_vtu(mesh, 'mesh'+s,
                  pointdata = {
                      'own': np.array(info[:,0]==b'o', dtype=np.int),
                      'frac': np.array(info[:,1]==b'y', dtype=np.int),
                  })
        fracs = np.loadtxt('fractures'+s, dtype=np.uint64)
        vtkw.write_vtu(
            vtkw.vtu_doc(vertices, fracs),
            'frac_mesh'+s+'.vtu',
        ) 

if __name__=='__main__':
    extract_meshes()

