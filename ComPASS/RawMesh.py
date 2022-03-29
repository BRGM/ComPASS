import os
import numpy as np
from vtkwriters import vtk_celltype


class CoC:
    def __init__(self, pointers, contents, subsizes):
        self.pointers = pointers
        self.contents = contents
        self.subsizes = subsizes

    @classmethod
    def from_array_of_arrays(self, arrayofarrays):
        subsizes = np.array([len(a) for a in arrayofarrays])
        return CoC(
            CoC.int_array(CoC.make_pointers(subsizes)),
            CoC.int_array(np.hstack(arrayofarrays)),
            subsizes,
        )

    @staticmethod
    def make_pointers(a):
        return np.cumsum(np.hstack([[0], a]))

    @staticmethod
    def int_array(a):
        return np.asarray(a, dtype=np.int32)

    @staticmethod
    def double_array(a):
        return np.asarray(a, dtype=np.double)

    def get_subsizes(self):
        if self.subsizes is None:
            self.subsizes = np.array(
                [
                    self.pointers[i + 1] - self.pointers[i]
                    for i in range(len(self.pointers) - 1)
                ]
            )
        return self.subsizes


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

    def convert_array_of_arrays(self, field):
        res = getattr(self, field)
        if isinstance(res, CoC):
            return res
        else:
            res = CoC.from_array_of_arrays(res)
            setattr(self, field, res)
            return res

    def get_vertices(self):
        if not isinstance(self.vertices, np.ndarray):
            self.vertices = CoC.double_array(self.vertices)
        return self.vertices

    def get_cell_nodes(self):
        return self.convert_array_of_arrays("cell_nodes")

    def get_cell_faces(self):
        return self.convert_array_of_arrays("cell_faces")

    def get_face_nodes(self):
        return self.convert_array_of_arrays("face_nodes")

    # cell and face types default to -1
    # try hint with simple geometries
    @staticmethod
    def fill_types(res, subsizes, known):
        res[:] = -1
        for p in known:
            res[subsizes == p[0]] = vtk_celltype[p[1]]

    def fill_cell_types(self, res):
        try:  # already provided
            res[:] = self.cell_types
        except AttributeError:
            RawMesh.fill_types(
                res,
                self.cell_nodes.get_subsizes(),
                [[4, "tet"], [8, "voxel"]],
            )

    def fill_face_types(self, res):
        try:  # already provided
            res[:] = self.face_types
        except AttributeError:
            RawMesh.fill_types(
                res,
                self.face_nodes.get_subsizes(),
                [[3, "triangle"], [4, "pixel"]],
            )

    @classmethod
    def convert(cls, other):
        info = {
            name: getattr(other, name)
            for name in ["vertices", "cell_nodes", "cell_faces", "face_nodes"]
        }
        return cls(**info)
