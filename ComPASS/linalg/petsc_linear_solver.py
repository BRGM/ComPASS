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
import sys
import yaml
from yaml.loader import SafeLoader

# This is a class for a yaml file with all the information we need to save
class data_part(yaml.YAMLObject):
    yaml_tag = "-------"

    def __init__(self, proc_rank, proc_data):
        self.proc = proc_rank
        self.proc_data = proc_data


class PetscLinearSystem:
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

        self.lsbuilder = simulation.LinearSystemBuilder()
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

            with open(f"{basename}/part_data.yaml", "w") as f:
                # Clearing the file if it already exists
                pass
            with open(f"{basename}/part_data.yaml", "a") as f:
                mpi.synchronize()
                Number_of_procs = "None"
                Block_size = "None"
                if mpi.is_on_master_proc:
                    yaml.dump(
                        {
                            "Number_of_procs": mpi.communicator().Get_size(),
                            "Block_size": self.lsbuilder.get_block_size(),
                        },
                        f,
                        default_flow_style=False,
                    )
                mpi.synchronize()
                yaml.dump(
                    data_part(
                        mpi.proc_rank,
                        {
                            "global_index_of_first_row": self.lsbuilder.get_rowstart(
                                mpi.proc_rank
                            ),
                            "local_number_of_rows": self.lsbuilder.get_non_zeros()[0][
                                0
                            ],
                            "local_number_of_wells": self.lsbuilder.get_n_wells(),
                        },
                    ),
                    f,
                    default_flow_style=False,
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


class PetscIterativeSolver(IterativeSolver):

    """
    A structure that holds an iterative PETSc KSP Object to solve the linear system
    """

    def __init__(
        self,
        linear_system,
        settings,
        pc=None,
        comm=PETSc.COMM_WORLD,
    ):

        self.activate_cpramg = None
        self.ksp = PETSc.KSP().create(comm=comm)
        self.settings = settings
        self.ksp.setOperators(linear_system.A, linear_system.A)
        self.pc = pc if pc is not None else CPRAMG(linear_system)
        self.ksp.setPC(self.pc)
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        self.ksp.setFromOptions()
        settings = IterativeSolverSettings(
            self.method, self.tolerance, self.max_iterations, settings.restart_size
        )
        super().__init__(linear_system, settings)

    def set_settings(self, settings):
        self.method = settings.method
        self.tolerance = settings.tolerance
        self.max_iterations = settings.max_iterations
        if self.method == "gmres":
            self.restart_size = settings.restart_size

    settings = property(
        fget=lambda self: IterativeSolverSettings(
            self.method,
            self.tolerance,
            self.max_iterations,
            None,
        ),
        fset=lambda self, value: self.set_settings(value),
        doc="Iterative solver settings",
    )
    method = property(
        fget=lambda self: self.ksp.getType(),
        fset=lambda self, type: self.ksp.setType(type),
        doc="Iterative method used for linear solving",
    )
    tolerance = property(
        fget=lambda self: self.ksp.rtol,
        fset=lambda self, value: self.ksp.setTolerances(rtol=value),
        doc="Relative decrease in the residual norm required for convergence",
    )
    max_iterations = property(
        fget=lambda self: self.ksp.max_it,
        fset=lambda self, value: self.ksp.setTolerances(max_it=value),
        doc="Maximum number of iterations accepted before convergence failure",
    )
    restart_size = property(
        # fget will raise AttributeError because restart value cannot be retrieved through petsc4py API",
        fset=lambda self, value: self.ksp.setGMRESRestart(value),
        doc="Number of iterations at which GMRES restarts",
    )

    def solve(self):

        self.ksp.setConvergenceHistory(self.max_iterations + 1)
        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        self.ksp_reason = self.ksp.getConvergedReason()
        self.nit = self.ksp.getIterationNumber()
        self.residual_history = self.ksp.history

        if self.ksp_reason < 0:  # Negative reason corresponds to convergence failure
            self.number_of_unsuccessful_iterations += self.nit
            raise IterativeSolverFailure(self.ksp_reason, self.nit)
        else:
            self.number_of_successful_iterations += self.nit

        return self.linear_system.x, self.nit

    def write_history(self, basename=""):
        """
        Writes the linear solver residual history in a file

        :param basename: base part of the file path
        :comm: MPI communicator
        """

        def dump_residual_log(self):
            with open(f"{basename}/solver_log.txt", "w") as f:
                f.write(f"{self}\n\n")
                f.write(f"Iteration, Residual norm\n")
                for it, residual in enumerate(self.residual_history):
                    f.write(f"{it:<9}, {residual:.13e}\n")

        mpi.on_master_proc(dump_residual_log)(self)

    def __str__(self):
        return f"{super().__str__()}\n   petsc4py new implementation\n   Settings : {self.settings}\n   Preconditioner : {self.pc}"


class PetscDirectSolver(DirectSolver):
    """
    A structure that holds a direct PETSc KSP Object to solve the linear system
    """

    def __init__(
        self,
        linear_system,
        comm=PETSc.COMM_WORLD,
    ):

        super().__init__(linear_system)
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.linear_system.A, self.linear_system.A)
        self.ksp.setType("preonly")
        self.ksp.getPC().setType("lu")
        self.ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)

    def solve(self):

        self.ksp.solve(self.linear_system.RHS, self.linear_system.x)
        self.ksp_reason = self.ksp.getConvergedReason()
        self.nit = 1

        if self.ksp_reason < 0:
            raise DirectSolverFailure(self.ksp_reason)

        return self.linear_system.x, self.nit

    def __str__(self):
        return f"{super().__str__()}\n   petsc4py new implementation"
