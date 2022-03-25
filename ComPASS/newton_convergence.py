#
# This file is part of kernel.
#
# kernel.is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

import numpy as np

# access to underlying MPI.py
from .mpi import MPI, is_on_master_proc
from ._kernel import get_kernel


class Legacy:
    def __init__(self, simulation):
        # FIXME: should be put elsewhere
        self.kernel = get_kernel()
        self.simulation = simulation
        # FIXME: this rely on the residuals being allocated elsewhere
        #        and associated with fortran code
        assert self.simulation.mesh_is_local
        self.residuals = self.simulation.Residuals()
        self.reference_pv = np.zeros(self.simulation.Residuals.npv(), dtype=np.double)
        self.reference_closure = 0

    def pv_norms(self):
        residuals = self.residuals
        local = np.linalg.norm(residuals.own_nodes, 1, axis=0)
        local += np.linalg.norm(residuals.own_fractures, 1, axis=0)
        local += np.linalg.norm(residuals.own_cells, 1, axis=0)
        local[0] += np.linalg.norm(residuals.own_injectors, 1)
        local[0] += np.linalg.norm(residuals.own_producers, 1)
        global_norms = np.zeros(self.simulation.Residuals.npv(), dtype=np.double)
        MPI.COMM_WORLD.Allreduce(local, global_norms, MPI.SUM)
        self._global_norms = global_norms
        return global_norms

    @property
    def latest_pv_norms(self):
        return self._global_norms

    def pv_norms_contributions(self):
        residuals = self.residuals
        npv = self.simulation.Residuals.npv()

        def reduce(v, n):
            V = np.zeros(n, dtype=np.double)
            MPI.COMM_WORLD.Allreduce(v, V, MPI.SUM)
            return V

        contributions = {}
        for name in ["own_nodes", "own_fractures", "own_cells"]:
            contributions[name] = reduce(
                np.linalg.norm(getattr(residuals, name), 1, axis=0), npv
            )
        for name in ["own_injectors", "own_producers"]:
            contributions[name] = reduce(np.linalg.norm(getattr(residuals, name), 1), 1)
        return contributions

    def closure_norm(self):
        kernel = self.kernel
        assert kernel
        closure = kernel.Residu_RelativeNorm_local_closure()
        closure = MPI.COMM_WORLD.allreduce(closure, MPI.SUM)
        self._closure = closure
        return closure

    @property
    def latest_closure_norm(self):
        return self._closure

    def reset_conservation_reference(self, dt):
        # global accumulation over the whole simulation domain (all procs involved)
        global_accumulation = self.simulation.total_accumulation(reset_states=False)
        global_reference = global_accumulation / (1000.0 * dt) + 1.0
        self.reference_pv = np.maximum(self.pv_norms(), global_reference)

    def reset_closure_reference(self):
        self.reference_closure = max(1.0, self.closure_norm())  # FIXME: ?

    def reset_references(self, dt):
        self.reset_conservation_reference(dt)
        self.reset_closure_reference()

    def relative_norm(self, display_contributions=False):
        assert self.reference_closure > 0
        assert np.all(
            self.reference_pv > 0
        ), "mininmum of reference value is %g" % np.min(self.reference_pv)
        closure_norm = self.closure_norm()
        if display_contributions:
            contributions = self.pv_norms_contributions()
            if is_on_master_proc:
                for name, value in contributions.items():
                    print(
                        f"{name} residual contribution: "
                        f"{value} ({value/self.reference_pv[:value.shape[0]]})"
                    )
                print(
                    "closure residual: "
                    f"{closure_norm} ({closure_norm/self.reference_closure})"
                )
        return max(
            closure_norm / self.reference_closure,
            np.max(self.pv_norms() / self.reference_pv),
        )
