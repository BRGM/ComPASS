from .__init__ import PETSc, mpi, np


class NullPC(object):
    """ A PC-Python Context which returns zero. Applied on the
    T/s field in the CPRAMG procedure """

    def setUp(self, pc):
        pass

    def apply(self, pc, x, y):
        y.set(0.0)


class BlockJacobi(PETSc.PC):
    def __init__(self, linear_system, comm=PETSc.COMM_WORLD):
        self.create(comm)
        self.setOperators(linear_system.A, linear_system.A)
        self.setType(PETSc.PC.Type.BJACOBI)
        self.setFactorLevels(1)

    def __repr__(self):
        return "Block Jacobi performing ILU(1) on each block"


class CPRAMG(PETSc.PC):
    """ PETSc implementation of the CPR-AMG preconditioning procedure"""

    def __init__(self, linear_system, amg_type="hypre", comm=PETSc.COMM_WORLD):
        # CPR-AMG is a multiplicative composite PC :
        # CPR-AMG = M2*(I - A*M1) + M1
        # Where M1 is an AMG procedure on the pressure unknowns,
        # set using the Petsc fieldsplit type.
        # and M2 = ILU(A)

        self.create(comm)
        self.setOperators(linear_system.A, linear_system.A)
        self.setType(PETSc.PC.Type.COMPOSITE)
        self.setCompositeType(PETSc.PC.CompositeType.MULTIPLICATIVE)
        self.addCompositePC(PETSc.PC.Type.FIELDSPLIT)
        self.addCompositePC(PETSc.PC.Type.BJACOBI)
        self.setUp()

        # We now define M1 as an additive fieldsplit pc
        # M1 = M11 + M12
        # with M11 an AMG v-cycle procedure on the pressure field
        # and M12 = 0 a Null PC which returns zero on the T/s field
        lsbuilder = linear_system.lsbuilder
        block_size = lsbuilder.get_block_size()
        (sizes, d_nnz, o_nnz) = lsbuilder.get_non_zeros()
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes
        n_wells = lsbuilder.get_n_wells()
        my_rowstart = lsbuilder.get_rowstart(mpi.proc_rank)

        non_well_nrowl = n_rowl - n_wells
        fs_pc = self.getCompositePC(0)
        p_indices = my_rowstart + np.arange(
            non_well_nrowl, step=block_size, dtype="int32"
        )
        p_indices = np.concatenate(
            (
                p_indices,
                my_rowstart
                + np.arange(start=non_well_nrowl, stop=n_rowl, dtype="int32"),
            )
        )

        # The Index Set for the pressure field
        p_IS = PETSc.IS().createGeneral(p_indices, comm=comm)
        # The Index Set for the temperature and saturation field
        # is simply the complement of the pressure IndexSet
        rest_IS = p_IS.complement(0, n_rowg)
        fs_pc.setFieldSplitIS(("pressure", p_IS), ("rest", rest_IS))
        fs_pc.setFieldSplitType(PETSc.PC.CompositeType.ADDITIVE)
        sub_ksp_list = fs_pc.getFieldSplitSubKSP()

        # Algebraic multigrid procedure on the pressure field
        pressure_ksp = sub_ksp_list[0]
        pressure_ksp.setType(PETSc.KSP.Type.PREONLY)
        pressure_pc = pressure_ksp.getPC()

        if amg_type not in ("gamg", "hypre"):
            mpi.master_print(
                f"Unknown AMG type in CPR-AMG preconditioner : {amg_type}; supported types : hypre, gamg\nUsing default hypre"
            )
            amg_type = "hypre"
        self.amg_type = amg_type

        if self.amg_type == "hypre":
            try:
                pressure_pc.setType(PETSc.PC.Type.HYPRE)
                PETSc.Options().setValue(
                    "-sub_0_fieldsplit_pressure_pc_hypre_boomeramg_strong_threshold",
                    0.5,
                )
            except PETSc.Error:  # If HYPRE is not available use default AMG from PETSc
                mpi.master_print(
                    "Hypre BoomerAMG is not available, using PETSc's GAMG procedure instead"
                )
        elif self.amg_type == "gamg":
            pressure_pc.setType(PETSc.PC.Type.GAMG)
            pressure_pc.setGAMGType("agg")

        rest_ksp = sub_ksp_list[1]
        rest_ksp.setType(PETSc.KSP.Type.PREONLY)
        null_pc = rest_ksp.getPC()
        null_pc.setType(PETSc.PC.Type.PYTHON)
        null_pc.setPythonContext(NullPC())

    def __repr__(self):
        return f"CPR-AMG preconditioner using {self.amg_type} AMG procedure on the pressure unknowns"
