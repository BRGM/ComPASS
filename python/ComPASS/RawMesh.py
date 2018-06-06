import os
import numpy as np


class RawMesh:

    def __init__(self, **kwds):
        assert 'vertices' in kwds
        assert 'cell_nodes' in kwds
        assert 'cell_faces' in kwds
        assert 'face_nodes' in kwds
        for name, value in kwds.items():
            setattr(self, name, value)
    
    @classmethod
    def convert(cls, other):
        info = { name: getattr(other, name)
                 for name in ['vertices', 'cell_nodes', 'cell_faces', 'face_nodes'] }
        return cls(**info)

