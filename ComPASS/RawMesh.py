import os
import numpy as np


class RawMesh:
    def __init__(self, **kwds):
        if len(kwds) > 0:
            self.build(**kwds)

    def build(self, **kwds):
        for param in ["vertices", "cell_nodes", "cell_faces", "face_nodes"]:
            if param not in kwds:
                raise RuntimeError(
                    f"Cannot build a RawMesh object without {param} argument."
                )
            if hasattr(self, param):
                raise RuntimeError(f"RawMesh object already has a {param} attribute.")
        for name, value in kwds.items():
            setattr(self, name, value)

    @classmethod
    def convert(cls, other):
        info = {
            name: getattr(other, name)
            for name in ["vertices", "cell_nodes", "cell_faces", "face_nodes"]
        }
        return cls(**info)
