import numpy as np
import MeshTools.vtkwriters as vtkw

def test_pvtu():
    vtkw.write_pvtu(
        vtkw.pvtu_doc(
            np.double,
            ['toto1_.vtu', 'toto2.vtu'],
            pointdata_types = {'toto': np.double, 'tutu': np.int64},
            celldata_types = {'toto': np.float32, 'tutu': np.int8},
        ),
        'mypvtu'
    )

if __name__=='__main__':
    test_pvtu()
