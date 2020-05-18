from petsc4py import PETSc
import numpy as np


class LinearSystem:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self, system):
        # CHECKME: is there a way to check PETSc has been initialized here
        assert system.kernel is not None
        self.system = system
        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()

    def setUp(self, simulation):
        self.createAMPI(simulation)
        self.setAMPI(simulation)
        self.setVecs(simulation)

    def solve(self):
        return self.system.kernel.SolvePetsc_ksp_solve(self.x)

    def check_solution(self):
        self.system.kernel.SolvePetsc_check_solution(self.x)

    def createAMPI(self, simulation):

        JacA = simulation.retrieve_jacobian()
        part = simulation.retrieve_partitioning()
        nb_node_own = part.nb_node_own
        nb_frac_own = part.nb_frac_own
        nb_node_local = part.nb_node_local
        nb_frac_local = part.nb_frac_local
        nb_comp_thermique = part.nb_comp_thermique
        nb_well_inj_own = part.nb_well_inj_own
        nb_well_prod_own = part.nb_well_prod_own
        nb_well_inj_local = part.nb_well_inj_local
        nb_well_prod_local = part.nb_well_prod_local
        nb_node_own_ncpus = simulation.nb_nodes_own()
        nb_frac_own_ncpus = simulation.nb_fractures_own()
        nb_well_inj_cpus = simulation.nb_wellinj_own()
        nb_well_prod_cpus = simulation.nb_wellprod_own()
        # print(nb_node_own_ncpus, nb_frac_own_ncpus, nb_well_inj_cpus, nb_well_prod_cpus)
        columns = JacA.columns() - 1
        row_offset = JacA.row_offset()

        # local row/col size:  node own and frac own
        n_rowl = (
            (nb_node_own + nb_frac_own) * nb_comp_thermique
            + nb_well_inj_own
            + nb_well_prod_own
        )
        n_coll = n_rowl

        # global row/col size: sum of all procs
        n_rowg = 0
        n_rowg = np.sum(
            (nb_node_own_ncpus + nb_frac_own_ncpus) * nb_comp_thermique
            + nb_well_inj_cpus
            + nb_well_prod_cpus
        )
        n_colg = n_rowg

        # number of nonzeros per row in diag or off-diag portion
        d_nnz = np.zeros(n_rowl, dtype="int32")
        o_nnz = np.zeros(n_rowl, dtype="int32")
        # print('proc %i has %i nodes own %i fractures own %i components + thermal'%(mpi.proc_rank, nb_node_own, nb_frac_own, nb_comp_thermique))
        # print('create sparse matrix on proc %i with structure %i (out of %i) x %i (out of %i )'%(mpi.proc_rank, n_rowl, n_rowg, n_coll, n_colg))

        # N: Node, F: Frac, WI: well inj, WP: well prod
        # l: local, o: own
        nl_fo = nb_node_local + nb_frac_own
        nl_fl = nb_node_local + nb_frac_local
        nl_fl_wio = nb_node_local + nb_frac_local + nb_well_inj_own
        nl_fl_wil = nb_node_local + nb_frac_local + nb_well_inj_local
        nl_fl_wil_wpo = (
            nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_own
        )
        nl_fl_wil_wpl = (
            nb_node_local + nb_frac_local + nb_well_inj_local + nb_well_prod_local
        )

        # rows associated with node own and frac own
        for i in range(nb_node_own + nb_frac_own):
            for j in range(row_offset[i], row_offset[i + 1]):

                lj = columns[j]

                if lj < nb_node_own:  # node own
                    d_nnz[
                        i * nb_comp_thermique : (i + 1) * nb_comp_thermique
                    ] += nb_comp_thermique

                elif lj >= nb_node_own and lj < nb_node_local:  # ghost node
                    o_nnz[
                        i * nb_comp_thermique : (i + 1) * nb_comp_thermique
                    ] += nb_comp_thermique

                elif lj >= nb_node_local and lj < nl_fo:  # frac own
                    d_nnz[
                        i * nb_comp_thermique : (i + 1) * nb_comp_thermique
                    ] += nb_comp_thermique

                elif lj >= nl_fo and lj < nl_fl:  # frac ghost
                    o_nnz[
                        i * nb_comp_thermique : (i + 1) * nb_comp_thermique
                    ] += nb_comp_thermique

                elif lj >= nl_fl and lj < nl_fl_wio:  # well inj own
                    d_nnz[i * nb_comp_thermique : (i + 1) * nb_comp_thermique] += 1

                elif lj >= nl_fl_wio and lj < nl_fl_wil:  # well inj ghost
                    o_nnz[i * nb_comp_thermique : (i + 1) * nb_comp_thermique] += 1

                elif lj >= nl_fl_wil and lj < nl_fl_wil_wpo:  # well prod own
                    d_nnz[i * nb_comp_thermique : (i + 1) * nb_comp_thermique] += 1

                else:  # well prod ghost
                    o_nnz[i * nb_comp_thermique : (i + 1) * nb_comp_thermique] += 1

        # rows associated with well inj own and well prod own
        start = (nb_node_own + nb_frac_own) * nb_comp_thermique

        for i in range(nb_well_inj_own + nb_well_prod_own):

            li = i + nb_node_own + nb_frac_own

            for j in range(row_offset[li], row_offset[li + 1]):

                lj = columns[j]

                if lj < nb_node_own:  # node own
                    d_nnz[start + i] += nb_comp_thermique

                elif nb_node_own <= lj and lj < nb_node_local:  # node ghost
                    o_nnz[start + i] += nb_comp_thermique

                elif nb_node_local <= lj and lj < nl_fo:  # frac own
                    d_nnz[start + i] += nb_comp_thermique

                elif nl_fo <= lj and lj < nl_fl:  # frac ghost
                    o_nnz[start + i] += nb_comp_thermique

                elif nl_fl <= lj and lj < nl_fl_wio:  # well inj own
                    d_nnz[start + i] += 1

                elif nl_fl_wio <= lj and lj < nl_fl_wil:  # well inj ghost
                    o_nnz[start + i] += 1

                elif nl_fl_wil <= lj and lj < nl_fl_wil_wpo:  # well prod own
                    d_nnz[start + i] += 1

                elif nl_fl_wil_wpo <= lj and lj < nl_fl_wil_wpl:  # well prod ghost
                    o_nnz[start + i] += 1
                else:
                    print("error in create ampi")

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))

    def setAMPI(self, simulation):

        JacA = simulation.retrieve_jacobian()
        part = simulation.retrieve_partitioning()
        rlg = part.rowl_to_rowg()
        clg = part.coll_to_colg()
        blocks = JacA.blocks()
        columns = JacA.columns() - 1
        row_offset = JacA.row_offset()
        block_size = JacA.block_size
        nb_rows = JacA.nb_rows
        nb_node_own = part.nb_node_own
        nb_frac_own = part.nb_frac_own
        nb_node_local = part.nb_node_local
        nb_frac_local = part.nb_frac_local
        nb_comp_thermique = part.nb_comp_thermique
        nb_well_inj_own = part.nb_well_inj_own
        nb_well_prod_own = part.nb_well_prod_own

        m, n = nb_comp_thermique, nb_comp_thermique
        idxm, idxn = (
            np.zeros(nb_comp_thermique, dtype="int32"),
            np.zeros(nb_comp_thermique, dtype="int32"),
        )
        arange = np.arange(nb_comp_thermique, dtype="int32")

        for i in range(nb_node_own + nb_frac_own):
            for j in range(row_offset[i], row_offset[i + 1]):

                row = rlg[i] - 1  # 0-based in petsc
                col = clg[columns[j]] - 1  # 0-based in petsc
                if columns[j] < nb_node_local + nb_frac_local:
                    idxm[:] = row + arange
                    idxn[:] = col + arange
                    self.A.setValues(idxm, idxn, blocks[j, :, :], addv=None)

                else:  # col is wellinj or wellprod, insert JacA.blocks[j,:,0]
                    idxm[:] = row + arange
                    idxn[0] = col

                    self.A.setValues(idxm, idxn[0], blocks[j, :, 0], addv=None)

        for i in range(
            nb_node_own + nb_frac_own,
            nb_node_own + nb_frac_own + nb_well_inj_own + nb_well_prod_own,
        ):
            for j in range(row_offset[i], row_offset[i + 1]):

                row = rlg[i] - 1  # 0-based in petsc
                col = clg[columns[j]] - 1  # 0-based in petsc

                # col is node or frac, insert JacA.blocks[j,0,:]
                if columns[j] < nb_node_local + nb_frac_local:
                    idxn[:] = col + arange
                    idxm[0] = row

                    self.A.setValues(idxm[0], idxn, blocks[j, 0, :], addv=None)

                else:  # col is wellinj or wellprod, insert JacA.blocks[j,0,0]
                    idxm[0] = row
                    idxn[0] = col

                    self.A.setValues(idxm[0], idxn[0], blocks[j, 0, 0], addv=None)

        self.A.assemblyBegin()
        self.A.assemblyEnd()

    def setVecs(self, simulation):

        JacRHS = simulation.retrieve_right_hand_side()
        part = simulation.retrieve_partitioning()
        nb_node_own = part.nb_node_own
        nb_frac_own = part.nb_frac_own
        nb_comp_thermique = part.nb_comp_thermique
        nb_well_inj_own = part.nb_well_inj_own
        nb_well_prod_own = part.nb_well_prod_own

        assert self.A.isAssembled(), "A must be assembled before calling setVecs"

        self.RHS = self.A.createVecs(side="right")
        self.x = self.A.createVecs(side="left")

        for i in range(nb_node_own + nb_frac_own):
            for j in range(nb_comp_thermique):

                self.RHS[i * nb_comp_thermique + j] = JacRHS[i, j]

        start = (nb_node_own + nb_frac_own) * nb_comp_thermique
        blockstart = nb_node_own + nb_frac_own

        for i in range(nb_well_inj_own + nb_well_prod_own):
            self.RHS[i + start] = JacRHS[i + blockstart, 0]

        self.RHS.assemblyBegin()
        self.RHS.assemblyEnd()
