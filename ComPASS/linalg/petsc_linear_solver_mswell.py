#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from .__init__ import mpi, PETSc
from .solver import *
from .exceptions import IterativeSolverFailure, DirectSolverFailure
from .preconditioners import CPRAMG


class PetscLinearSystemMSWell:
    """
    A structure that holds and manages the linear system as
    Petsc objects
    """

    def __init__(self, simulation):
        """
        :param simulation: an initialised simulation object to get the data from
        """
        self.x = PETSc.Vec()
        self.A = PETSc.Mat()
        self.RHS = PETSc.Vec()

        self.lsbuilder = simulation.LinearSystemBuilderMSWells()
        (sizes, d_nnz, o_nnz) = self.lsbuilder.get_non_zeros()
        n_rowl, n_rowg = sizes
        n_coll, n_colg = sizes

        assert d_nnz.shape == (n_rowl,)
        assert o_nnz.shape == (n_rowl,)

        self.A.createAIJ(size=((n_rowl, n_rowg), (n_coll, n_colg)), nnz=(d_nnz, o_nnz))
        self.x.createMPI((n_rowl, n_rowg))
        self.RHS.createMPI((n_rowl, n_rowg))

    def dump_ascii(self, basename="", comm=PETSc.COMM_WORLD):
        """
        Writes the linear system (Matrix, solution and RHS) in three different files in ASCII format

        :param basename: common part of the file names
        :comm: MPI communicator
        """
        makeviewer = PETSc.Viewer().createASCII

        def dump_item(item, name):
            viewer = makeviewer(f"{basename}/{name}.dat", "w", comm)
            item.view(viewer)
            viewer.destroy

        mpi.master_print(">> Linear system dump")
        dump_item(self.A, "A")
        dump_item(self.RHS, "RHS")
        dump_item(self.x, "x")

    def dump_binary(self, basename="", comm=PETSc.COMM_WORLD):
        """
        Writes the linear system (Matrix, solution and RHS) in three different files in binary format

        :param basename: common part of the file names
        :comm: MPI communicator
        """
        makeviewer = PETSc.Viewer().createBinary

        def dump_item(item, name):
            viewer = makeviewer(f"{basename}/{name}.dat", "w", comm)
            item.view(viewer)
            viewer.destroy

        def dump_part_data(self):
            with open(f"{basename}/part_data.txt", "w") as f:
                # Clearing the file if it already exists
                pass
            with open(f"{basename}/part_data.txt", "a") as f:
                mpi.synchronize()
                if mpi.is_on_master_proc:
                    f.write(f"Number of procs : {mpi.communicator().Get_size()}\n")
                    f.write(f"Block size : {self.lsbuilder.get_block_size()}\n")
                mpi.synchronize()
                f.write(
                    f"\nProc rank : {mpi.proc_rank}\n"
                    f"Global index of first row : {self.lsbuilder.get_rowstart(mpi.proc_rank)}\n"
                    f"Local number of rows : {self.lsbuilder.get_non_zeros()[0][0]}\n"
                    f"Local number of wells : {self.lsbuilder.get_n_wells()}\n"
                )

        mpi.master_print(">> Linear system dump")
        dump_part_data(self)
        dump_item(self.A, "A")
        dump_item(self.RHS, "RHS")
        dump_item(self.x, "x")

    def check_residual_norm(self):

        """
        Displays the residual norm (1-Norm, 2-Norm and infinity norm) for convergence check
        """

        y = self.RHS.duplicate()
        self.RHS.copy(y)  # y = b
        y.scale(-1.0)  # y = -b
        self.A.multAdd(self.x, y, y)  # y = Ax-b
        mpi.master_print("Linear solution check ||Ax-b||")
        norm = y.norm(PETSc.NormType.NORM_1)
        mpi.master_print("  1-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_2)
        mpi.master_print("  2-Norm       ", norm)
        norm = y.norm(PETSc.NormType.NORM_INFINITY)
        mpi.master_print("  Infinity norm", norm)

    def set_from_jacobian(self):

        self.lsbuilder.set_AMPI(self.A)
        self.lsbuilder.set_RHS(self.RHS)
