#
# This file is part of ComPASS.
#
# ComPASS is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
import numpy as np
from . import mpi
from .mpi import MPI as MPI
from .newton_convergence import Legacy as LegacyConvergence
from ._kernel import get_kernel
from .simulation import AlignmentMethod
from .linalg.factory import linear_solver
from .linalg.exceptions import LinearSolverFailure
from .callbacks import InterruptTrigger
from .messages import error
from .linalg import PETSc

NewtonStatus = namedtuple("NewtonStatus", ["newton_iterations", "linear_iterations"])
NewtonLoopTick = namedtuple(
    "NewtonLoopTick", ["timeloop_tick", "current_dt", "newton", "iteration"]
)


class NewtonFailure(Exception):
    def __init__(self, status):
        self.status = status


class IterationExhaustion(NewtonFailure):
    def __init__(self, status):
        super().__init__(status)


def dump_start_info(iteration):
    # if (commRank == 0) then
    #    if (iteration == 1) then
    #        do i = 1, size(fd)
    #            j = fd(i)
    #            write (j, *) ""
    #            write (j, '(A)', advance='no') "     *Residu init conv:  "
    #            do k = 1, NbCompThermique
    #            write (j, '(A,ES12.5)', advance='no') "  ", NewtonResConvInit(k)
    #            end do
    #            write (j, *) ""
    #            write (j, '(A,ES12.5)') "     *Residu init clos:    ", NewtonResClosInit
    #            write (j, *) ""
    #        end do
    #    end if
    #    do i = 1, size(fd)
    #        j = fd(i)
    #        write (j, '(A,I4)', advance="no") "     *Newton iter:", iteration
    #        write (j, '(A,E18.10)', advance="no") "      Newton res norm:", NewtonResNormRel(iteration)
    #    end do
    # end if
    pass


class Newton:
    """
    A structure that manages the Newton loop.
    """

    def __init__(
        self,
        simulation,
        tol,
        maxit,
        lsolver,
        maxit_before_decrease=None,
        convergence_scheme=None,
        iteration_callbacks=None,
        failure_callbacks=None,
    ):
        """
        :param simulation: The simulation object.
        :param tol: The tolerance used for convergence.
        :param maxit: The maximum number of newton iteration.
        :param maxit_before_decrease: after this number of iterations, residual must decrease or failure
        :param lsolver: The linear solver to be used (cf. :class:`ComPASS.petsc.LinearSolver`).
        :param convergence_scheme: The convergence scheme to be used (defaults to :class:`ComPASS.newton_convergence.Legacy`).
        """
        # FIXME: this is transitory
        self.simulation = simulation
        self.lsolver = lsolver
        self.tolerance = tol
        self.maximum_number_of_iterations = maxit
        maxit_before_decrease = maxit_before_decrease or maxit
        self.maxit_before_decrease = max(1, maxit_before_decrease)
        self.failures = 0
        self.lsolver_iterations = []
        self.number_of_useless_iterations = 0
        self.number_of_succesful_iterations = 0
        self.relative_residuals = None
        self.increments = get_kernel().NewtonIncrements()
        self.increments.init()
        if not convergence_scheme:
            self.convergence_scheme = LegacyConvergence(simulation)
        self.check_well_errors_at_convergence = False
        self.iteration_callbacks = iteration_callbacks or ()
        self.failure_callbacks = failure_callbacks or ()

        # when the values are incremented thanks to the solution
        # X^{m+1} = X^{m} + alpha * dX, alpha <= 1
        # alpha is determined to have an impact on the different unknowns
        # less or equal to incre_obj_u
        self.incre_obj_P = 5.001e6
        self.incre_obj_T = 20.0
        self.incre_obj_C = 1.0
        self.incre_obj_S = 0.2

        self.iterations = 0

        # create the timeloop_log_callback object, will contains info to store convegence data
        # store it to avoid erasing the file at each standard time loop call
        # (erase it only wheen creating a new newton object)
        self.timeloop_log_callback = None

    def reset_loop(self):
        kernel = get_kernel()
        kernel.IncCV_LoadIncPreviousTimeStep()

    def init_iteration(self):
        kernel = get_kernel()
        # Enforce Dirichlet values
        kernel.DirichletContribution_update()
        # phase pressures are needed in IncPrimSecd_update_secondary_dependencies (fugacities)
        kernel.LoisThermoHydro_compute_phase_pressures()
        # Update only Well Pressures (well pressure drops are kept constant here)
        kernel.IncCVWells_UpdateWellPressures()
        #        mpi.master_print('init iteration - compute thermo')
        # Update local jacobian contributions (closure laws)
        kernel.IncPrimSecd_update_secondary_dependencies()  # FIXME: this is needed to update globals used in LoisThermoHydro_compute
        kernel.LoisThermoHydro_compute_phase_pressures_derivatives()  # needs to be done IncPrimSecd_update_secondary_dependencies (needs NumIncTotalPrim)
        kernel.LoisThermoHydro_compute()
        #        mpi.master_print('init iteration - compute fluxes')
        kernel.Flux_DarcyFlux_Cell()
        kernel.Flux_DarcyFlux_Frac()
        if kernel.has_energy_transfer_enabled():
            kernel.Flux_FourierFlux_Cell()
            kernel.Flux_FourierFlux_Frac()

    def Jacobian_ComputeJacSm_alignment(self):
        kernel = get_kernel()
        alignment_method = self.simulation.alignment
        if alignment_method == AlignmentMethod.inverse_diagonal:
            kernel.Jacobian_Alignment_diag()
        elif alignment_method == AlignmentMethod.manual:
            kernel.Jacobian_Alignment_man()
        else:
            error("Unknown alignment method.")

    def Jacobian_ComputeJacSm(self, dt):
        kernel = get_kernel()
        kernel.Jacobian_JacBigA_BigSm(dt)
        kernel.Jacobian_Regularization()
        kernel.Jacobian_Schur()
        self.Jacobian_ComputeJacSm_alignment()

    def increment(self, x):
        kernel = get_kernel()
        #        mpi.master_print('increment variables')
        ghosts_synchronizer = self.simulation.info.ghosts_synchronizer
        assert ghosts_synchronizer is not None
        ghosts_synchronizer.synchronize(x)
        # mpi.master_print('retrieve solutions')
        ghosts_synchronizer.retrieve_solution(self.increments)
        # mpi.master_print('nodes increment shape', self.increments.nodes().shape)
        # mpi.master_print(self.increments.nodes())
        #        mpi.master_print('lois thermo hydro')
        # mpi.master_print('cells increment shape', self.increments.cells().shape)
        # mpi.master_print('before\n', self.increments.cells())
        kernel.Jacobian_GetSolCell(self.increments)
        # mpi.master_print('after\n', self.increments.cells())
        kernel.IncPrimSecd_PrimToSecd(self.increments)
        Pmax = self.incre_obj_P
        Tmax = self.incre_obj_T
        Cmax = self.incre_obj_C
        Smax = self.incre_obj_S
        relaxation = kernel.Newton_compute_relaxation(
            self.increments, Pmax, Tmax, Cmax, Smax
        )
        # if relaxation < 1.0:
        #     mpi.master_print("relaxation:", relaxation)

        kernel.IncCV_NewtonIncrement(self.increments, relaxation)
        kernel.DirichletContribution_update()
        # phase pressures are needed in IncPrimSecd_update_secondary_dependencies (fugacities)
        kernel.LoisThermoHydro_compute_phase_pressures()
        #        mpi.master_print('flash all volumes')
        kernel.NN_flash_all_control_volumes()

    def loop(self, timeloop_tick, dt, display_contributions=False):
        if PETSc.Sys.getVersion() >= (3, 18):
            # PETSc garbage collection was introduced in version 3.18
            # https://petsc.org/release/changes/318
            # FIXME: Ideally, we should not call garbage collection by ourselves!
            #        Yet, when using CPR-AMG preconditioner with Hypre
            #        garbage collection is not managed correctly and memory
            #        footprint increases until overflow
            #        cf. issue https://gitlab.inria.fr/compass/v4/ComPASS/-/issues/572.
            #        The manual call to PETSC.garbage_cleanup() is a collective call
            #        made once for each Newton loop.
            #        Our assumption is that this is not too expensive
            #        and currently protects us from potential overflows which do not
            #        seem restricted to Hypre bindings
            #        (e.g. https://gitlab.com/petsc/petsc/-/issues/1309).
            PETSc.garbage_cleanup()
        kernel = get_kernel()
        convergence_scheme = self.convergence_scheme
        assert convergence_scheme is not None
        relative_residuals = []
        self.relative_residuals = relative_residuals
        lsolver = self.lsolver
        self.lsolver_iterations = []
        self.status = NewtonStatus(1, self.lsolver_iterations)
        self.init_iteration()
        # CHECKME: does this need to be done after newton_init_iteration?
        kernel.Residu_reset_history()
        kernel.Residu_compute(dt)
        convergence_scheme.reset_references(dt)
        mpi.master_print(
            f"timestep mass residual(s) scale(s): "
            f"{' '.join(f'{x:.5g}' for x in convergence_scheme.reference_pv[:-1])}\n"
            f"             energy residual scale: {convergence_scheme.reference_pv[-1]:.5g}\n"
            f"  closure equations residual scale: {convergence_scheme.reference_closure:.5g}"
        )

        for iteration in range(self.maximum_number_of_iterations):

            newton_tick = NewtonLoopTick(timeloop_tick, dt, self, iteration + 1)
            self.Jacobian_ComputeJacSm(dt)
            lsolver.linear_system.set_from_jacobian()

            try:
                x, nit = lsolver.solve()
            except LinearSolverFailure as e:
                self.number_of_useless_iterations += iteration + 1
                raise e

            self.lsolver_iterations += [nit]
            self.increment(x)
            self.init_iteration()
            kernel.Residu_compute(dt)
            relative_residuals.append(
                convergence_scheme.relative_norm(display_contributions)
            )
            mpi.master_print(
                f"Newton {iteration + 1: 3d} residuals: ",
                f"mass {' '.join(f'{x:12.5g}' for x in convergence_scheme.latest_pv_norms[:-1])} ",
                f"energy {convergence_scheme.latest_pv_norms[-1]:12.5g}",
                f" closure {convergence_scheme.latest_closure_norm:12.5g} "
                if convergence_scheme.latest_closure_norm > 0
                else "",
                f"relative (conv. crit.) {relative_residuals[-1]:12.5g}",
            )
            self.status = NewtonStatus(iteration + 1, self.lsolver_iterations)
            for callback in self.iteration_callbacks:
                callback(newton_tick)
            if relative_residuals[-1] < self.tolerance:
                if self.check_well_errors_at_convergence:
                    self.check_well_residuals()
                self.number_of_succesful_iterations += iteration + 1
                return self.status
            # after a certain amount of iterations, test that the residu decreases
            # otherwise failure
            if iteration > self.maxit_before_decrease:
                if relative_residuals[-1] / relative_residuals[-2] > 1.0:
                    self.number_of_useless_iterations += iteration + 1
                    raise IterationExhaustion(self.status)

        self.number_of_useless_iterations += iteration + 1
        raise IterationExhaustion(self.status)

    def check_well_residuals(self):
        simulation = self.simulation
        residuals = self.convergence_scheme.residuals
        compute_error = lambda a: 0 if len(a) == 0 else np.fabs(a).max()

        def well_errors(well_data, well_residuals):
            nb_own_wells = well_residuals.shape[0]

            def well_errors_category(code):
                category = np.array(
                    [data.operating_code == code for data in well_data[:nb_own_wells]]
                )
                return compute_error(well_residuals[category])

            max_flowrate_error = well_errors_category("f")
            max_wellhead_error = well_errors_category("p")
            return max_flowrate_error, max_wellhead_error

        all_well_errors = np.vstack(
            [
                well_errors(list(simulation.injectors_data()), residuals.own_injectors),
                well_errors(list(simulation.producers_data()), residuals.own_producers),
            ]
        )
        all_well_errors = np.max(all_well_errors, axis=0)
        global_well_errors = np.zeros(2, dtype=np.double)
        MPI.COMM_WORLD.Allreduce(all_well_errors, global_well_errors, MPI.MAX)
        if mpi.is_on_master_proc:
            print("Maximum well errors:")
            dq, dp = global_well_errors
            print(f"  {dq} for well operatinf on flowrate")
            print(f"  {dp} for well operatinf on pressure")


def default_Newton(simulation, tolerance=1e-5, max_iterations=8, lsolver=None):
    if lsolver is None:
        # The default lsolver is a legacy iterative solver
        # which uses the CPR-AMG preconditioner
        lsolver = linear_solver(
            simulation,
            legacy=False,
            direct=False,
            from_options=True,
        )
    return Newton(simulation, tolerance, max_iterations, lsolver)
