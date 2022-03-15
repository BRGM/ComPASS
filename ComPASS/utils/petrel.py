from pathlib import Path
import numpy as np
from petrelgridio.petrel_grid import PetrelGrid as _PetrelGrid

from ..RawMesh import RawMesh
from .. import mpi


class PetrelGrid:
    def __init__(self, path):
        if mpi.is_on_master_proc:
            path = Path(path)
            if not path.exists():
                raise FileNotFoundError(path)
            self.load(path)
            return
        self.mesh = None
        self.fracture_faces = None

    def load(self, path):
        grid = _PetrelGrid.build_from_files(path)
        hexa = grid.process()[0]
        vertices, cell_faces, face_nodes, fracture_faces = grid.process_faults(
            hexa, return_splitted_faces=True
        )
        cell_nodes = [
            np.unique(np.hstack([face_nodes[f] for f in faces])) for faces in cell_faces
        ]
        self.grid = grid
        self.mesh = RawMesh(
            vertices=vertices,
            cell_nodes=cell_nodes,
            face_nodes=face_nodes,
            cell_faces=cell_faces,
        )
        self.fracture_faces = fracture_faces

    @property
    def permeability(self):
        k = self.grid.permx, self.grid.permy, self.grid.permz
        assert k[0].shape == k[1].shape
        assert k[0].shape == k[2].shape
        shape = k[0].shape
        result = np.zeros((np.prod(shape), 3, 3), dtype=np.double)
        for i in range(3):
            result[:, i, i] = np.asarray(np.ravel(k[i]), dtype=np.double)
        result *= 1e-15  # to mdarcy
        return result
