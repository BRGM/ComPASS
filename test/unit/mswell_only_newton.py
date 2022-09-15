import numpy as np
import ComPASS
from ComPASS.utils.units import *
import ComPASS.io.mesh as io
import ComPASS.dump_wells as dw
import ComPASS.mpi as mpi
from mpi4py import MPI
import os
from ComPASS._kernel import get_kernel

from ComPASS.linalg.petsc_linear_solver_mswell import PetscLinearSystemMSWell
from ComPASS.linalg.solver import IterativeSolverSettings
from ComPASS.linalg.petsc_linear_solver import PetscDirectSolver, PetscIterativeSolver

from ComPASS.distributed_system_mswells import DistributedSystemMSWells
from ComPASS.ghosts.synchronizer import Synchronizer

################################################################################
class MSWellsNewtonParam:
    def __init__(
        self,
        t0,
        tf,
        dt,
        dtmax,
        newtonMaxIter=40,
        TotalMaxIter=2000,
        crit_resrel=1.0e-9,
        crit_dxmax=1.0e-10,
        verbose=False,
    ):
        #########################################################################
        # Newton and solver parameters
        self.t0 = t0
        self.tf = tf
        self.dt = dt
        self.dtmax = dtmax
        self.newtonMaxIter = 40
        self.TotalMaxIter = TotalMaxIter
        self.crit_resrel = crit_resrel
        self.crit_dxmax = crit_dxmax
        self.verbose = verbose


##########################################################################################


class MSWellsNewton:
    def __init__(self, newton_parameters, simulation):
        self.kernel = get_kernel()
        self.residuals = self.kernel.ResidualsMSWells()
        self.increments = self.kernel.NewtonIncrements()
        self.increments.init()
        self.distro_system = DistributedSystemMSWells(self.kernel)
        self.linear_system = None
        self.newton_parameters = newton_parameters
        self.simulation = simulation
        self.dxmax = 1.0
        # Postprocessing
        if mpi.is_on_master_proc:
            self.f_reswell = open("reswell_c.val", "wt")

    def increment(self):
        ghosts_synchronizer = Synchronizer(self.distro_system, True)
        x = self.linear_system.x
        ghosts_synchronizer.synchronize(x)
        ghosts_synchronizer.retrieve_solution(self.increments)
        self.kernel.IncPrimSecdMSWells_PrimToSecd(self.increments)
        relaxation = self.kernel.Newton_compute_relaxation_mswells(self.increments)
        # relaxation = 1.0
        self.dxmax = self.kernel.Newton_get_last_max_inc_mswells()
        self.kernel.IncCV_NewtonIncrement(self.increments, relaxation)
        if self.newton_parameters.verbose:
            mpi.master_print("   Relaxation factor/dxmax:", relaxation, self.dxmax)
        self.kernel.NN_flash_all_control_volumes()

    def init_iteration(self):
        self.kernel.IncPrimSecdMSWells_compute()
        self.kernel.LoisThermoHydroMSWells_compute()  # it must be compute again after VSHydroMSWells_init
        self.kernel.VSHydroMSWells_compute()

    def solve(self):

        # Perform Additional Initializations
        self.kernel.mswells_copy_states_from_reservoir()
        #######################################################################
        # For reinit (also modify ResiduMSWells_allocate)
        # also do not call self.kernel.VSHydroMSWells_init() below
        # self.kernel.IncCVMSWells_init_from_file ()
        #######################################################################
        self.kernel.mswells_init_edge_data()
        self.kernel.IncPrimSecdMSWells_compute()
        self.kernel.LoisThermoHydroMSWells_compute()
        self.kernel.mswells_init_leaf_data()
        self.kernel.VSHydroMSWells_init()  # Computes initial pressure, don't call this when reinit from file

        # time iteration
        t = self.newton_parameters.t0
        dt = self.newton_parameters.dt
        total_iters = 0
        while t < self.newton_parameters.tf:

            self.kernel.IncCV_SaveIncPreviousTimeStep()
            self.init_iteration()
            self.kernel.ResiduMSWells_reset_history()
            self.kernel.ResiduMSWells_compute(dt)
            self.kernel.JacobianMSWells_ComputeJacSm(dt, True)  # Perform alignment
            res0 = 0
            resrel = 1.0
            mpi.master_print("------- Time iteration at t,dt:", t, dt, "-------")
            t = t + dt
            for iteration in range(self.newton_parameters.newtonMaxIter):
                ##if nprocs == 1 and total_iters == 465:
                ##        print("MPI rank:", mpi.proc_rank)
                ##        print("Process id:", pid)
                ##        self.kernel.init_wait_for_debug(True)

                pcsp = np.copy(self.simulation.mswell_node_states().p)
                pcsT = np.copy(self.simulation.mswell_node_states().T)
                pcsS = np.copy(self.simulation.mswell_node_states().S)
                total_iters = total_iters + 1

                self.linear_system = PetscLinearSystemMSWell(self.kernel)
                self.linear_system.set_from_jacobian()
                # Save Jacobian, Residual data
                # self.kernel.JacobianMSWells_print_LA_info_to_file(dt,total_iters-1)#directive _DEBUG_JAC_MSWELLS_ needs to be set to 1 in fortran file
                ##    #######################################################################
                ##    # Direct Solver
                linear_solver = PetscDirectSolver(self.linear_system)
                ##    # settings = IterativeSolverSettings(
                ##    #    "gmres",
                ##    #    1e-6,
                ##    #    150,
                ##    #    30,
                ##    # )
                ##    # linear_solver=  PetscIterativeSolver(self.linear_system,settings)
                ##    #######################################################################
                linear_solver.solve()
                # Save Petsc lsystem data
                # self.linear_system.dump_ascii(".")
                if self.newton_parameters.verbose:
                    mpi.master_print(
                        "   Newton ", (iteration + 1), "\t total_it:", total_iters
                    )
                self.increment()
                res_norm = np.linalg.norm(self.residuals.own_mswell_nodes, 1, axis=0)
                global_res_norm = np.zeros(res_norm.shape[0], dtype=np.double)
                MPI.COMM_WORLD.Allreduce(res_norm, global_res_norm, MPI.SUM)
                if iteration == 0:
                    res0 = np.sum(global_res_norm)
                else:
                    resrel = np.sum(global_res_norm) / res0
                if self.newton_parameters.verbose:
                    mpi.master_print(
                        "           residuals",
                        global_res_norm,
                        "res0",
                        res0,
                        "rel",
                        resrel,
                    )

                if mpi.is_on_master_proc:
                    self.f_reswell.write(
                        str(total_iters)
                        + str("\t")
                        + str(np.sum(global_res_norm))
                        + str("\t")
                        + str(resrel)
                        + str("\n")
                    )
                self.init_iteration()
                self.kernel.ResiduMSWells_compute(dt)
                self.kernel.JacobianMSWells_ComputeJacSm(dt, True)  # Perform alignment
                mpi.synchronize()
                if (
                    resrel < self.newton_parameters.crit_resrel
                    or self.dxmax < self.newton_parameters.crit_dxmax
                ):
                    break

            if (iteration == (self.newton_parameters.newtonMaxIter - 1)) and (
                resrel > self.newton_parameters.crit_resrel
            ):
                mpi.master_print("--------------------------------------------")
                mpi.master_print("------ Warning:  Newton not converged!!-------")
                t = t - dt
                dt = 0.5 * dt
                self.kernel.IncCV_LoadIncPreviousTimeStep()
                if total_iters > (self.newton_parameters.TotalMaxIter - 1):
                    exit()
            else:
                mpi.master_print(
                    "--- ** Succed, t,dt, nllew_iter, tnew_iter, resrel, :",
                    t - dt,
                    dt,
                    iteration + 1,
                    total_iters,
                    resrel,
                )
                if self.newton_parameters.verbose:
                    allreduce = mpi.communicator().allreduce
                    mpi.master_print(
                        "max p variation",
                        allreduce(
                            np.fabs(
                                self.simulation.mswell_node_states().p - pcsp
                            ).max(),
                            mpi.MPI.MAX,
                        ),
                    )
                    mpi.master_print(
                        "max T variation",
                        allreduce(
                            np.fabs(
                                self.simulation.mswell_node_states().T - pcsT
                            ).max(),
                            mpi.MPI.MAX,
                        ),
                    )
                    mpi.master_print(
                        "max S variation",
                        allreduce(
                            np.fabs(
                                self.simulation.mswell_node_states().S - pcsS
                            ).max(),
                            mpi.MPI.MAX,
                        ),
                    )
                dt = min(1.1 * dt, self.newton_parameters.dtmax)
                self.kernel.IncCVMSWells_print_info_to_file()  # directive _DEBUG_INCCVM_MSWELLS_ needs to be set to 1 in fortran file
                self.kernel.JacobianMSWells_print_IP_info_to_file(
                    t
                )  # directive _DEBUG_JAC_IP_MSWELLS_ needs to be set to 1 in fortran file
            if mpi.is_on_master_proc:
                self.f_reswell.write(str("\n"))

    #########################################################################
