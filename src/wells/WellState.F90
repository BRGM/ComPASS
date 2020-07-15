!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module WellState

   use CommonMPI, only: CommonMPI_abort
   use DefModel, only: &
      NbPhase, NbComp, MCP, NumPhasePresente_ctx, NbPhasePresente_ctx
   use IncCVWells, only: PerfoWellProd, NodebyWellProdLocal, NodebyWellInjLocal
   use IncCVReservoir, only: IncNode
   use LoisThermoHydro, only: &
      DensitemolaireKrViscoCompNode, &
      DensitemolaireKrViscoEnthalpieNode
   use Thermodynamics, only: f_Enthalpie
   use MeshSchema, only: NodeDatabyWellProdLocal
   use IncPrimSecd, only: IncPrimSecd_compPrim_nodes

   implicit none

   ! FIXME: the following globals are bound to disapear
   ! molar/energy fluxes for production well(s)
   double precision, allocatable, dimension(:, :) :: summolarFluxProd !< Molar flux for production well: summolarFluxProd(component,node_well)
   double precision, allocatable, dimension(:) :: sumnrjFluxProd !< Energy flux for production well: sumnrjFluxProd(node_well)

   public :: &
      WellState_allocate, &
      WellState_free, &
      WellState_FlowrateWellProd, &
      WellState_solve_for_temperature

contains

   subroutine WellState_allocate
      integer :: Nnz
      Nnz = NodebyWellProdLocal%Pt(NodebyWellProdLocal%Nb + 1)
      allocate (summolarFluxProd(NbComp, Nnz))
      allocate (sumnrjFluxProd(Nnz))
   end subroutine WellState_allocate

   subroutine WellState_free()

      deallocate (summolarFluxProd)
      deallocate (sumnrjFluxProd)

   end subroutine WellState_free

   !> \brief Compute the flowrate at each node of the production well:
   !! Fill the CSR vectors summolarFluxProd and sumnrjFluxProd.
   !! Use PerfoWellProd\%Pression and IncNode\%Pression.
   subroutine WellState_FlowrateWellProd

      integer :: k, s, nums, icp, m, mph, sparent
      double precision :: Pws, Ps, WIDws
      double precision:: Flux_ks(NbComp), perf_inflow
#ifdef _THERMIQUE_
      double precision:: FluxT_ks
#endif

#ifndef NDEBUG
      if (.not. allocated(summolarFluxProd)) &
         call CommonMPI_abort("Work array should be allocated.")
      if (.not. allocated(sumnrjFluxProd)) &
         call CommonMPI_abort("Work array should be allocated.")
#endif

      summolarFluxProd = 0.d0
      sumnrjFluxProd = 0.d0

      do k = 1, NodebyWellProdLocal%Nb

#ifndef NDEBUG
         if (NodeDatabyWellProdLocal%Val(NodebyWellProdLocal%Pt(k + 1))%PtParent /= -1) &
            call CommonMPI_abort("Wrong head")
#endif

         ! parents_idxs are greater that their sons_idxs
         do s = NodebyWellProdLocal%Pt(k) + 1, NodebyWellProdLocal%Pt(k + 1)

            nums = NodebyWellProdLocal%Num(s)
            sparent = NodeDatabyWellProdLocal%Val(s)%PtParent
            Pws = PerfoWellProd(s)%Pression ! P_{w,s}
            Ps = IncNode(nums)%Pression ! P_s
            WIDws = NodeDatabyWellProdLocal%Val(s)%WID

#ifndef NDEBUG
            if (WIDws < 0.d0) &
               call CommonMPI_abort("Well index is negative.")
#endif

            perf_inflow = WIDws*(Ps - Pws)

            if (perf_inflow > 0.d0) then
               Flux_ks = 0.d0
#ifdef _THERMIQUE_
               FluxT_ks = 0.d0
#endif
               do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
                  mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)
                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i
                        Flux_ks(icp) = Flux_ks(icp) + DensiteMolaireKrViscoCompNode(icp, m, nums)*perf_inflow
                     end if
                  end do
#ifdef _THERMIQUE_
                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieNode(m, nums)*perf_inflow
#endif
               end do
               summolarFluxProd(:, s) = summolarFluxProd(:, s) + Flux_ks(:)
#ifdef _THERMIQUE_
               sumnrjFluxProd(s) = sumnrjFluxProd(s) + FluxT_ks
#endif
            end if

            !now update the parent node, except the well root node
            if (sparent /= -1) then ! head node if sparent = -1
               summolarFluxProd(:, sparent) = summolarFluxProd(:, sparent) + summolarFluxProd(:, s)
#ifdef _THERMIQUE_
               sumnrjFluxProd(sparent) = sumnrjFluxProd(sparent) + sumnrjFluxProd(s)
#endif
#ifndef NDEBUG
            else
               if (s /= NodebyWellProdLocal%Pt(k + 1)) &
                  call CommonMPI_abort("Head found in the middle of the well.")
#endif
            end if

         end do

      end do

   end subroutine WellState_FlowrateWellProd

   ! FIXME: we should change the API to pass E/n instead of E and n
#ifdef NDEBUG
   pure &
#endif
      subroutine WellState_solve_for_temperature(Phase, E, p, T, n, converged, ResT, newton_tol_in, newton_itermax_in)

      integer, intent(in) :: Phase
      double precision, intent(in) :: E !< total energy
      double precision, intent(in) :: p !< pressure
      double precision, intent(inout) :: T !< initial temperature and result
      double precision, intent(in) :: n !< total number of moles
      logical, intent(inout) :: converged
      double precision, intent(inout):: ResT
      double precision, intent(in), optional :: newton_tol_in
      integer, intent(in), optional :: newton_itermax_in
      double precision :: H, dHdP, dHdT, dHdC(NbComp), dHdS(NbPhase), newton_tol
      double precision :: dummyCi(NbComp), dummySat(NbPhase) ! not used by f_Enthalpie !
      integer :: i, newton_itermax

      newton_tol = 1.0d-5
      if (present(newton_tol_in)) newton_tol = newton_tol_in
      newton_itermax = 1000
      if (present(newton_itermax_in)) newton_itermax = newton_itermax_in
      converged = .false.
      ResT = 1E+10
      do i = 1, newton_itermax
         call f_Enthalpie(Phase, p, T, dummyCi, dummySat, H, dHdP, dHdT, dHdC, dHdS)
         ResT = E - H*n ! residual
         if (abs(ResT) < newton_tol) then
            converged = .true.
            exit
         else
            T = T + ResT/(dHdT*n)
         end if
      end do

   end subroutine WellState_solve_for_temperature

end module WellState
