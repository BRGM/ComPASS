#
# This file is part of kernel.
#
# kernel.is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

import ComPASS
# access to underlying mpi4py
from ComPASS.mpi import MPI as mpi

class Legacy:

    def __init__(self):
        # FIXME: should be put elsewhere
        assert ComPASS.kernel
        self.kernel = ComPASS.kernel
        # FIXME: this rely on the residuals being allocated elsewhere
        #        and associated with fortran code
        assert ComPASS.mesh_is_local
        self.residuals = ComPASS.Residuals()
        self.reference_conservation = np.zeros(
            ComPASS.Residuals.npv(), dtype=np.double
        )
        self.reference_closure = 0

    def norms(self):
        residuals = self.residuals
        res = np.linalg.norm(residuals.own_nodes, 1, axis=0)
        res+= np.linalg.norm(residuals.own_fractures, 1, axis=0)
        res+= np.linalg.norm(residuals.own_cells, 1, axis=0)
        res[0]+= np.linalg.norm(residuals.own_injectors, 1)
        res[0]+= np.linalg.norm(residuals.own_producers, 1)
        return res
    
    def reset_conservation_reference(self, dt):
        ref = np.zeros(ComPASS.Residuals.npv(), dtype=np.double)
        for states in [
            ComPASS.own_node_states(),
            ComPASS.own_fracture_states(),
            ComPASS.own_cell_states()
        ]:
            ref+= np.linalg.norm(states.accumulation, 1, axis=0)
        ref/= 1000. * dt
        ref+= 1. # FIXME: ?
        conservation = np.maximum(self.norms(), ref)
        mpi.COMM_WORLD.Allreduce(conservation, self.reference_conservation, mpi.MAX)
    
    def reset_closure_reference(self):
        kernel = self.kernel
        assert kernel
        closure = kernel.Residu_RelativeNorm_local_closure()
        closure = max(1., closure) # FIXME: ?
        self.reference_closure =  mpi.COMM_WORLD.allreduce(closure, mpi.MAX)
    
    def reset_references(self, dt):
        self.reset_conservation_reference(dt)
        self.reset_closure_reference()
    
    def relative_norm(self):
        assert self.reference_closure>0
        assert np.all(self.reference_conservation>0)
        kernel = self.kernel
        assert kernel
        local = max(
            kernel.Residu_RelativeNorm_local_closure() / self.reference_closure,
            np.max(self.norms()/self.reference_conservation)
        )
        return mpi.COMM_WORLD.allreduce(local, mpi.MAX)
