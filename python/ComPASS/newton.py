#
# This file is part of kernel.
#
# kernel.is free software: you can redistribute it and/or modify it under both the terms
# of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
# and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
#

from collections import namedtuple
import ComPASS
import ComPASS.mpi as mpi


class LinearSolver:
    
    def __init__(self, tol, maxit, restart=None):
        self.reset(tol, maxit, restart)
        self.failures = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0

    def reset(self, tol, maxit, restart=None):
        self.tolerance = tol
        self.maximum_number_of_iterations = maxit
        self.restart_iteration = maxit if restart is None else restart
        ComPASS.kernel.SolvePetsc_Ksp_configuration(
           self.tolerance, self.maximum_number_of_iterations,
           self.restart_iteration
       )

NewtonStatus = namedtuple(
    'NewtonStatus', 
    ['newton_iterations', 'linear_iterations']
)

class NewtonFailure(Exception):
    def __init__(self, status):
        self.status = status

class KspFailure(NewtonFailure):
    def __init__(self, status, reason):
        super().__init__(status)
        self.reason = reason

class IterationExhaustion(NewtonFailure):
    def __init__(self, status):
        super().__init__(status)


def dump_start_info(iteration):
    #if (commRank == 0) then
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
    #end if
    pass


class Newton:
    
    def __init__(self, tol, maxit, lsolver):
        self.tolerance = tol
        self.maximum_number_of_iterations = maxit
        self.lsolver = lsolver
        self.failures = 0
        self.number_of_useful_linear_iterations = 0
        self.number_of_succesful_iterations = 0
        self.number_of_useless_iterations = 0
        self.relative_residuals = None
        self.increments = ComPASS.kernel.NewtonIncrements()
        self.increments.init()

    def reset_loop(self):
        kernel = ComPASS.kernel
        assert ComPASS.kernel
        kernel.IncCV_LoadIncPreviousTimeStep()
        kernel.IncCVWells_PressureDrop()

    def init_iteration(self):
        kernel = ComPASS.kernel
        assert ComPASS.kernel
        kernel.DirichletContribution_update()
        kernel.IncCVWells_PressureDrop()
#        mpi.master_print('init iteration - compute thermo')
        kernel.LoisThermoHydro_compute()
#        mpi.master_print('init iteration - compute fluxes')
        kernel.Flux_DarcyFlux_Cell()
        kernel.Flux_DarcyFlux_Frac()
        if kernel.has_energy_transfer_enabled():
            kernel.Flux_FourierFlux_Cell()
            kernel.Flux_FourierFlux_Frac()
    
    def increment(self):
        kernel = ComPASS.kernel
        assert ComPASS.kernel
#        mpi.master_print('increment variables')
        kernel.SolvePetsc_Sync()
        # mpi.master_print('retrieve solutions')
        kernel.SolvePetsc_GetSolNodeFracWell(self.increments)
        # mpi.master_print('nodes increment shape', self.increments.nodes().shape)
        # mpi.master_print(self.increments.nodes())
#        mpi.master_print('lois thermo hydro')
        # mpi.master_print('cells increment shape', self.increments.cells().shape)
        # mpi.master_print('before\n', self.increments.cells())
        kernel.Jacobian_GetSolCell(self.increments)
        # mpi.master_print('after\n', self.increments.cells())
        kernel.LoisThermoHydro_PrimToSecd(self.increments)
        relaxation = kernel.IncCVReservoir_NewtonRelax(self.increments)
        # mpi.master_print('relaxation:', relaxation) 
        kernel.IncCV_NewtonIncrement(self.increments, relaxation)
        kernel.DirichletContribution_update()
#        mpi.master_print('flash all volumes')
        kernel.NN_flash_all_control_volumes()
    
    def loop(self, dt, simulation_context):
        kernel = ComPASS.kernel
        assert ComPASS.kernel
        relative_residuals = []
        self.relative_residuals = relative_residuals
        initial_conservation_residual =  kernel.CTVector()
        conservation_residual = kernel.CTVector()
        initial_closure_residual = None
        lsolver = self.lsolver
        nb_lsolver_iterations = 0
        total_lsolver_iterations = 0
        self.init_iteration()
        # CHECKME: does this need to be done after newton_init_iteration?
        kernel.Residu_reset_history()
        kernel.Residu_compute(dt)
        kernel.Residu_RelativeNorm_local_conservation(conservation_residual)
        kernel.Residu_RelativeNorm_initial_conservation(
            dt, conservation_residual, initial_conservation_residual
        )
        closure_residual = kernel.Residu_RelativeNorm_local_closure()
        initial_closure_residual = kernel.Residu_RelativeNorm_initial_closure(closure_residual)
        mpi.master_print('initial residuals (reference)',
                initial_conservation_residual.as_array(), initial_closure_residual)
        mpi.master_print('                    residuals',
                    conservation_residual.as_array(), closure_residual)
        for iteration in range(self.maximum_number_of_iterations):
            kernel.Jacobian_ComputeJacSm(dt)
            kernel.SolvePetsc_SetUp()
            ksp_status = kernel.SolvePetsc_KspSolve()
            # mpi.master_print('KSP status', ksp_status) 
            if not ComPASS.activate_direct_solver:
                nb_lsolver_iterations = kernel.SolvePetsc_KspSolveIterationNumber()
                # mpi.master_print('with', nb_lsolver_iterations, 'linear iterations')
                total_lsolver_iterations+= nb_lsolver_iterations
                mpi.master_print(
                    'linear iterations:', 
                    kernel.SolvePetsc_Ksp_iterations()
                )
            if ksp_status<0:
                lsolver.failures+= 1
                if not ComPASS.activate_direct_solver:
                    lsolver.number_of_useless_iterations+= nb_lsolver_iterations
                raise KspFailure(
                    NewtonStatus(iteration, total_lsolver_iterations),
                    ksp_status,
                )
            # kernel.SolvePetsc_check_solution()
            lsolver.number_of_succesful_iterations+= nb_lsolver_iterations
            self.increment()
            self.init_iteration()
            kernel.Residu_compute(dt)
            kernel.Residu_RelativeNorm_local_conservation(conservation_residual)
            closure_residual = kernel.Residu_RelativeNorm_local_closure()
            relative_residuals.append(
                kernel.Residu_compute_relative_norm(
                    initial_conservation_residual, conservation_residual,
                    initial_closure_residual, closure_residual,
                )
            )
            mpi.master_print('Newton % 3d          residuals'%(iteration + 1),
                        conservation_residual.as_array(), closure_residual,
                        'rel', relative_residuals[-1])
            if relative_residuals[-1] < self.tolerance:
                return NewtonStatus(iteration+1, total_lsolver_iterations)
        mpi.master_print('Newton relative residuals:')
        for i, r in enumerate(relative_residuals):
            mpi.master_print('%02d: %15.9e' % (i, r))
        raise IterationExhaustion(NewtonStatus(iteration, total_lsolver_iterations))
    

