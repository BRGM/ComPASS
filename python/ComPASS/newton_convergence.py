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
        self.reference_pv = np.zeros(
            ComPASS.Residuals.npv(), dtype=np.double
        )
        self.reference_closure = 0

    def pv_norms(self):
        residuals = self.residuals
        local = np.linalg.norm(residuals.own_nodes, 1, axis=0)
        local+= np.linalg.norm(residuals.own_fractures, 1, axis=0)
        local+= np.linalg.norm(residuals.own_cells, 1, axis=0)
        local[0]+= np.linalg.norm(residuals.own_injectors, 1)
        local[0]+= np.linalg.norm(residuals.own_producers, 1)
        global_norms = np.zeros(ComPASS.Residuals.npv(), dtype=np.double)
        mpi.COMM_WORLD.Allreduce(local, global_norms, mpi.SUM)
        return global_norms
    
    def closure_norm(self):
        kernel = self.kernel
        assert kernel
        closure = kernel.Residu_RelativeNorm_local_closure()
        return mpi.COMM_WORLD.allreduce(closure, mpi.SUM)
    
    def reset_conservation_reference(self, dt):
        local = np.zeros(ComPASS.Residuals.npv(), dtype=np.double)
        for states in [
            ComPASS.own_node_states(),
            ComPASS.own_fracture_states(),
            ComPASS.own_cell_states()
        ]:
            local+= np.linalg.norm(states.accumulation, 1, axis=0)
        local/= 1000. * dt
        global_reference = np.zeros(ComPASS.Residuals.npv(), dtype=np.double)
        mpi.COMM_WORLD.Allreduce(local, global_reference, mpi.SUM)
        global_reference+= 1.
        self.reference_pv = np.maximum(self.pv_norms(), global_reference)
    
    def reset_closure_reference(self):
        self.reference_closure  = max(1., self.closure_norm()) # FIXME: ?
    
    def reset_references(self, dt):
        self.reset_conservation_reference(dt)
        self.reset_closure_reference()
    
    def relative_norm(self):
        assert self.reference_closure>0
        assert np.all(self.reference_pv>0), (
            'mininmum of reference value is %g' % np.min(self.reference_pv)
        )
        return max(
            self.closure_norm() / self.reference_closure,
            np.max(self.pv_norms() / self.reference_pv),
        )
