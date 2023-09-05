!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

#include "gravity_average_rho.def"

module Jacobian
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use InteroperabilityStructures
   use mpi
   use CommonType
   use CommonMPI
   use DefModel
   use LoisThermoHydro
   use NumbyContext
   use Physics
   use Newton
   use SchemeParameters
   use IncCVReservoir
   use IncCVReservoirTypes
   use IncCVWells
   use VAGFrac
   use NumbyContext
   use IncPrimSecd
   use MeshSchema
   use Flux
   use Residu
   use MeshSchemaMSWells
   use JacobianMSWells, &
      MSWells_JacA => JacA, MSWells_Sm => Sm
   use IncCVMSWells
#else
   use iso_c_binding, only: c_double, c_int, c_bool
   use InteroperabilityStructures, only: &
      csr_block_matrix_wrapper, retrieve_csr_block_matrix, &
      cpp_array_wrapper_dim2, retrieve_double_array_dim2
   use mpi, only: MPI_Abort
   use CommonType, only: CSRArray2dble, CSR
   use CommonMPI, only: &
      commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use DefModel, only: &
      NumPhasePresente_ctx, NbPhasePresente_ctx, LIQUID_PHASE, &
      NbComp, NbPhase, NbCompThermique, NbContexte, MCP, aligmat, NbIncTotalPrimMax, &
      phase_can_be_present

#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      GAS_PHASE, LIQUID_PHASE, DIPHASIC_CONTEXT
#endif

   use LoisThermoHydro, only: &
      DensiteMolaireKrViscoCompWellInj, DensiteMolaireKrViscoEnthalpieWellInj, &
      DensiteMolaireKrViscoCompNode, DensiteMolaireKrViscoCompCell, DensiteMolaireKrViscoCompFrac, &
      DensiteMolaireKrViscoEnthalpieNode, DensiteMolaireKrViscoEnthalpieCell, DensiteMolaireKrViscoEnthalpieFrac, &
      divDensiteMolaireKrViscoCompNode, divDensiteMolaireKrViscoCompCell, divDensiteMolaireKrViscoCompFrac, &
      divDensiteMolaireKrViscoEnthalpieNode, divDensiteMolaireKrViscoEnthalpieCell, divDensiteMolaireKrViscoEnthalpieFrac, &
      divTemperatureNode, divTemperatureCell, divTemperatureFrac, &
      SmTemperatureNode, SmTemperatureCell, SmTemperatureFrac, &
      DensiteMassiqueNode, DensiteMassiqueCell, DensiteMassiqueFrac, &
      SmDensiteMassiqueNode, SmDensiteMassiqueCell, SmDensiteMassiqueFrac, &
      SmPressionNode, SmPressionCell, SmPressionFrac, &
#ifdef _WITH_FREEFLOW_STRUCTURES_
      divFreeFlowMolarFlowrateNode, SmFreeFlowMolarFlowrateNode, &
      FreeFlowMolarFlowrateCompNode, divFreeFlowMolarFlowrateCompNode, SmFreeFlowMolarFlowrateCompNode, &
      FreeFlowHmCompNode, divFreeFlowHmCompNode, SmFreeFlowHmCompNode, &
      divFreeFlowHTTemperatureNetRadiationNode, SmFreeFlowHTTemperatureNetRadiationNode, &
      FreeFlowMolarFlowrateEnthalpieNode, divFreeFlowMolarFlowrateEnthalpieNode, SmFreeFlowMolarFlowrateEnthalpieNode, &
      AtmEnthalpieNode, &
#endif
      SmDensiteMolaireKrViscoCompNode, SmDensiteMolaireKrViscoCompCell, SmDensiteMolaireKrViscoCompFrac, &
      SmDensiteMolaireKrViscoEnthalpieNode, SmDensiteMolaireKrViscoEnthalpieCell, SmDensiteMolaireKrViscoEnthalpieFrac, &
      SmDensiteMolaireEnergieInterneSatNode, SmDensiteMolaireEnergieInterneSatCell, SmDensiteMolaireEnergieInterneSatFrac, &
      divDensiteMassiqueNode, divDensiteMassiqueCell, divDensiteMassiqueFrac, &
      divDensiteMolaireKrViscoEnthalpieWellInj, divDensiteMolaireKrViscoCompWellInj, &
      divDensiteMolaireSatCompNode, divDensiteMolaireSatCompCell, divDensiteMolaireSatCompFrac, &
      divDensiteMolaireEnergieInterneSatNode, divDensiteMolaireEnergieInterneSatCell, divDensiteMolaireEnergieInterneSatFrac, &
      divPressionNode, divPressionCell, divPressionFrac, &
      divPhasePressureNode, divPhasePressureCell, divPhasePressureFrac, &
      SmDensiteMolaireSatComp, divSaturationNode, divSaturationCell, divSaturationFrac

   use NumbyContext, only: &
      NbCompCtilde_ctx, NumCompCtilde_ctx, NbIncPTC_ctx

   use Physics, only: gravity, CpRoche

   use Newton, only: Newton_increments_pointers, Newton_increments, Newton_pointers_to_values
   use SchemeParameters, only: eps

   use IncCVReservoir, only: &
      IncNode, IncCell, IncFrac

   use IncCVReservoirTypes, only: TYPE_IncCVReservoir

   use IncCVWells, only: &
      PerfoWellInj, DataWellInjLocal, &
      PerfoWellProd, PerfoWellInj
   use VAGFrac, only: &
      TkLocal_Darcy, TkLocal_Fourier, TkFracLocal_Darcy, TkFracLocal_Fourier, &
      VolDarcy, &
      PoroVolDarcy, &
      PoroVolFourier, &
      Poro_1VolFourier

   use NumbyContext, only: &
      NbIncTotalPrim_ctx, NbIncTotalPrim_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   use IncPrimSecd, only: &
      NumIncTotalPrimNode, NumIncTotalPrimCell, NumIncTotalPrimFrac

   use MeshSchema, only: &
      IdNodeLocal, &
#ifdef _WITH_FREEFLOW_STRUCTURES_
      IsFreeflowNode, AtmState, &
#endif
      NodebyCellLocal, FracbyCellLocal, NodebyFaceLocal, &
      NodebyWellProdLocal, NodeDatabyWellProdLocal, &
      NodebyWellInjLocal, NodeDatabyWellInjLocal, &
      NbCellLocal_Ncpus, NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbFracOwn_Ncpus, NbWellProdLocal_Ncpus, NbFaceOwn_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellInjOwn_Ncpus, NbWellProdOwn_Ncpus, DataWellProdLocal, &
      CellbyNodeOwn, NodebyFracOwn, NodebyNodeOwn, FracbyNodeOwn, FracbyFracOwn, &
      WellInjbyNodeOwn, WellProdbyNodeOwn, CellbyFracOwn, &
      NbFracCellMax, NbNodeCellMax, NbNodeFaceMax, &
      FracToFaceLocal, XNodeLocal, XCellLocal, XFaceLocal, &
      SurfFreeFlowLocal, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus, &
      NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbMSWellNodeOwn_Ncpus, NbMSWellNodeLocal_Ncpus, &
      MSWellbyNodeOwn

   use Flux, only: FluxDarcyKI, FluxDarcyFI

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY
   use Flux, only: epsilon_avrho
#endif

   use Residu, only: &
      ResiduNode, ResiduCell, ResiduFrac, &
      ResiduWellInj, ResiduWellProd

   use MeshSchemaMSWells, only: &
      NbMSWellLocal, NbMSWellNodeOwn, &
      NbMSWellNodeLocal, &
      IncIdxNodebyMSWellLocal
   use JacobianMSWells, only: &
      MSWells_JacA => JacA, MSWells_Sm => Sm
   use IncCVMSWells, only: &
      IncMSWell, IncCVMSWells_compute_coupling

#endif
   implicit none

   !            node, frac, cell, wellinj, wellprod, mswellnodes local
   !           | a11,  a12, a13, a14, a15,  a16 | node own
   ! JacBigA = | a21,  a22, a23, 0,   0,     0  | frac own
   !           | a31,  a32, a33, 0,   0,     0  | cell (own and ghost)
   !           ! a41,  0,   0,   a44, 0,     0  | wellinj own
   !           ! a51,  0,   0,   0,   a55,   0  | wellprod own
   !           ! a61 , 0,   0,   0,   0,   a66  | mswellnodes own

   !            node, frac, wellinj, wellprod, mswellnodes
   ! JacA =    | A11, A12, A13, A14,  A15 | node own
   !           | A21, A22, 0,   0,     0  | frac own
   !           | A31, 0,   A33, 0,     0  | wellinj own,
   !           | A41, 0,   0,   A44,   0  | wellprod own
   !           | A51, 0,   0,   0,    A55 | mswellnodes own

   type(CSRArray2dble), public, target :: JacBigA
   type(CSRArray2dble), public, target :: JacA

   ! right hand side before and after Schur complement
   real(c_double), pointer, dimension(:, :), public :: &
      bigSm, Sm

   ! pivot(:,k): pivot indices for cell k
   integer, allocatable, dimension(:, :), private :: pivot

   ! JacBigA is a sparse matrix.
   ! for an element in JacBigA, we need to known its num in JacBigA%Val.
   ! the following vectors are used for this purpose
   integer, allocatable, dimension(:), private :: csrK, csrSR

   !Offsets for mswell_nodes
   !where does the index of the  MSWell block matrix starts
   integer mswells_JacBigA_row_off, mswells_JacA_row_off
   integer mswells_JacBigA_num_off, mswells_JacA_num_off
   !Offsets for local mswell_nodes variable index enumeration (to be used for cols index)
   integer mswells_JacBigA_unk_lcl_off, mswells_JacA_unk_lcl_off

   public :: &
      Jacobian_StrucJacBigA, & !< non-zero structure of Jacobian before Schur
      Jacobian_StrucJacA, & !< non-zero structure of Jacobian after Schur
      Jacobian_JacBigA_BigSm, & !< Jacobian and right hand side
      Jacobian_Schur, & !< Schur complement
      Jacobian_Regularization, &
      Jacobian_Alignment_diag, &
      Jacobian_Alignment_man, &
      Jacobian_free

   private :: &
      ! the fill of Jacobian/ right hand side is decomposed into three subroutines
      !   (1) Jacobian_JacBigA_BigSm_accumulation: term n_k(X_j^n)
      !   (2) Jacobian_JacBigA_BigSm_cell: loop of cell
      !      (2.1) nodes in cell
      !      (2.2) fracs in cell
      !   (3) Jacobian_JacBigA_BigSm_frac: loop of frac
      Jacobian_JacBigA_BigSm_accumulation, &  ! (1)
      Jacobian_JacBigA_BigSm_cell, &  ! (2)
      Jacobian_JacBigA_BigSm_frac, &  ! (3)
      Jacobian_JacBigA_BigSm_wellinj, &  !
      Jacobian_JacBigA_BigSm_wellprod, &  !
      Jacobian_JacBigA_BigSm_mswells, &    !Initialize mswells states if the mswell is closed
      Jacobian_JacBigA_BigSm_coupling_mswells, &
      !
      ! In (2), we compute div( DensiteMolaire*Kr*Visco*Comp)*V_{k,s}
      !                  + DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
      !                  + div( DensiteMolaire*Kr*Visco*Enthalpie*FluxFourier ) if thermique
      !         V_{k,s} is Darcy Flux: k is cell, s is nodes/fracs; (2.1)
      !                                k is frac, s is nodes (2.2)
      !
      ! To compute DensiteMolaire*Kr*Visco*Comp*div(V_{k,s}),
      !    first div(pho) in the loop of k,s
      !    then  div(V_{k,s}) in the loop of k,s
      !    last  DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
      !
      ! div(DensiteMolaire*Kr*Visco*Comp)*V_{k,s}
      Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux, &  ! (k is cell, s is node) or (k is cell, s is frac) or (k is frac, s is node)
      ! DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
      Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux, &
      Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_mph, &
      ! div(V_{k,s})
      Jacobian_divDarcyFlux, &
      Jacobian_divDarcyFlux_mph, &
      ! div(rho)
      Jacobian_divrho_gravity, &
      ! For thermique, we compute div(FluxFourier)
      ! then div( DensiteMolaire*Kr*Visco*Enthalpie*FluxFourier )
      Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux, &
      Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_mph, &
      ! div(FourierFlux)
      Jacobian_divFourierFlux, &
#ifdef _WITH_FREEFLOW_STRUCTURES_
      ! div(FreeFlow)
      Jacobian_JacBigA_BigSm_FF_node, &
      Jacobian_divMolarFreeFlow_node, &
      Jacobian_divThermalFreeFlow_node, &
#endif
      Jacobian_RowCol_KSR, &
      Jacobian_RowCol_FR, &
      Jacobian_Regularization_row, &
      Jacobian_Alignment_diag_row, &
      Jacobian_Alignment_man_row, &
      Jacobian_JacBigA_locate_frac_row

contains

   subroutine retrieve_jacobian(A) &
      bind(C, name="retrieve_jacobian")
      type(csr_block_matrix_wrapper), intent(inout) :: A

      call retrieve_csr_block_matrix(JacA, A)

   end subroutine retrieve_jacobian

   subroutine retrieve_right_hand_side(a) &
      bind(C, name="retrieve_right_hand_side")
      type(cpp_array_wrapper_dim2), intent(inout) :: a

      call retrieve_double_array_dim2(Sm, a)

   end subroutine retrieve_right_hand_side

   subroutine retrieve_big_jacobian(A) &
      bind(C, name="retrieve_big_jacobian")
      type(csr_block_matrix_wrapper), intent(inout) :: A

      call retrieve_csr_block_matrix(JacBigA, A)

   end subroutine retrieve_big_jacobian

   subroutine retrieve_big_right_hand_side(a) &
      bind(C, name="retrieve_big_right_hand_side")
      type(cpp_array_wrapper_dim2), intent(inout) :: a

      call retrieve_double_array_dim2(bigSm, a)

   end subroutine retrieve_big_right_hand_side

   subroutine dump_jacobian(specific_row, specific_col)

      integer, optional, intent(in) :: specific_row, specific_col
      integer :: s, start

      do s = 1, NbNodeOwn_Ncpus(commRank + 1)
         if (present(specific_row)) then
            if (s /= specific_row) cycle
         endif
         call dump_row(s, specific_col)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1)
      do s = 1, NbFracOwn_Ncpus(commRank + 1)
         if (present(specific_row)) then
            if (s + start /= specific_row) cycle
         endif
         call dump_row(s + start, specific_col)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)
      do s = 1, NbCellLocal_Ncpus(commRank + 1)
         if (present(specific_row)) then
            if (s + start /= specific_row) cycle
         endif
         call dump_row(s + start, specific_col)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
              + NbCellLocal_Ncpus(commRank + 1)
      do s = 1, NbWellInjOwn_Ncpus(commRank + 1)
         if (present(specific_row)) then
            if (s + start /= specific_row) cycle
         endif
         call dump_row(s + start, specific_col)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
              + NbCellLocal_Ncpus(commRank + 1) + NbWellInjOwn_Ncpus(commRank + 1)
      do s = 1, NbWellProdOwn_Ncpus(commRank + 1)
         if (present(specific_row)) then
            if (s + start /= specific_row) cycle
         endif
         call dump_row(s + start, specific_col)
      end do

   end subroutine dump_jacobian

   subroutine dump_row(s, specific_col)

      integer, intent(in) :: s
      integer, optional, intent(in) :: specific_col
      integer :: i, j, n
      double precision :: a

      do n = JacBigA%Pt(s) + 1, JacBigA%Pt(s + 1)
         if (present(specific_col)) then
            if (JacBigA%Num(n) /= specific_col) cycle
         endif
         write (*, *) 'sites', s, JacBigA%Num(n), 'nonzero', n, 'on proc', commRank
         do j = 1, NbCompThermique
            do i = 1, NbCompThermique
               a = JacBigA%Val(i, j, n)
               if (abs(a) <= 0.d0 .or. abs(a - 1.d0) < 1e-20) then
                  write (*, "(I15)", advance="no") int(a)
               else
                  write (*, "(E15.7)", advance="no") a
               end if
            end do
            write (*, *)
         end do
      end do

   end subroutine dump_row

   !> \brief Fill BigSm with the residual values
   !!
   !!   If the node is dir, we suppose that the residu is 0,
   !!     in other words, we set its right hand side as 0.
   !!
   !!   The Jacobian corresponding to dir nodes are Id matrix,
   !!     setting bigSm(dirichlet)=0 or not has no influence mathematically,
   !!     however, it could have influence for linear solver
   !!     since the values could be very large 10**7.
   !!     It is observed when debugging!!!
   subroutine Jacobian_init_RHS_from_residual()

      integer :: j, start

      do j = 1, NbNodeOwn_Ncpus(commRank + 1)
         if (IdNodeLocal(j)%P == "d") then
            bigSm(:, j) = 0.d0
         else
            bigSm(:, j) = -ResiduNode(:, j)
         end if
      end do

      start = NbNodeOwn_Ncpus(commRank + 1)
      do j = 1, NbFracOwn_Ncpus(commRank + 1)
         bigSm(:, j + start) = -ResiduFrac(:, j)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)
      do j = 1, NbCellLocal_Ncpus(commRank + 1)
         bigSm(:, j + start) = -ResiduCell(:, j)
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
              + NbCellLocal_Ncpus(commRank + 1)
      do j = 1, NbWellInjOwn_Ncpus(commRank + 1)
         bigSm(1, j + start) = -ResiduWellInj(j)   ! Pressure
         bigSm(2:NbCompThermique, j + start) = 0.d0 ! not used
      end do

      start = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
              + NbCellLocal_Ncpus(commRank + 1) + NbWellInjOwn_Ncpus(commRank + 1)
      do j = 1, NbWellProdOwn_Ncpus(commRank + 1)
         bigSm(1, j + start) = -ResiduWellProd(j)  ! Pressure
         bigSm(2:NbCompThermique, j + start) = 0.d0 ! not used
      end do

   end subroutine Jacobian_init_RHS_from_residual

   !> \brief fill Jacobian and right hand side before Schur: main subroutine
   !> 1. init right hand side
   !> 2. div prim and Sm for term \f$n_k(X_j^n)\f$
   !> 3.1 loop of cell
   !> 3.2 loop of frac
   !> 3.3 loop of well inj
   !> 3.4 loop of well prod
   !> 3.5 loop of FreeFlow Nodes
   subroutine Jacobian_JacBigA_BigSm(Delta_t) &
      bind(C, name="Jacobian_JacBigA_BigSm")

      real(c_double), intent(in) :: Delta_t

      JacBigA%Val = 0.d0

      call Jacobian_init_RHS_from_residual

      call Jacobian_JacBigA_BigSm_accumulation(Delta_t)

      call Jacobian_JacBigA_BigSm_cell

      call Jacobian_JacBigA_BigSm_frac

      call Jacobian_JacBigA_BigSm_wellinj

      call Jacobian_JacBigA_BigSm_wellprod

#ifdef _WITH_FREEFLOW_STRUCTURES_
      call Jacobian_JacBigA_BigSm_FF_node
#endif
      !loop of mswells
      call Jacobian_JacBigA_BigSm_mswells
      if (IncCVMSWells_compute_coupling()) then
         !Compute contribution  with coupling
         call Jacobian_JacBigA_BigSm_coupling_mswells
      endif

      call Jacobian_JacBigA_BigSm_dirichlet_nodes

      ! to print the jacobian, uncomment the following line
      ! call dump_jacobian

   end subroutine Jacobian_JacBigA_BigSm

   ! FIXME: if we had the same number of dof in rows and cols we would not need col argument
   subroutine Jacobian_JacBigA_BigSm_find_diagonal_nz(J, row, col, nz)
      type(CSRArray2dble), intent(in) :: J
      integer, intent(in) :: row
      integer, intent(in) :: col
      integer, intent(out) :: nz

#ifndef NDEBUG
      logical :: found_diagonal = .false.
#endif

      do nz = J%Pt(row) + 1, J%Pt(row + 1)
         if (J%Num(nz) == col) then
#ifndef NDEBUG
            found_diagonal = .true.
#endif
            exit
         endif
      end do

#ifndef NDEBUG
      if (nz /= J%Pt(row) + findloc(J%Num((J%Pt(row) + 1):J%Pt(row + 1)), col, 1)) &
         call CommonMPI_abort("Jacobian_JacBigA_BigSm_find_diagonal_nz something went wrong")
      if (.not. found_diagonal) &
         call CommonMPI_abort("Jacobian_JacBigA_BigSm_find_diagonal_nz failed")
#endif

   end subroutine Jacobian_JacBigA_BigSm_find_diagonal_nz

   subroutine Jacobian_JacBigA_BigSm_dirichlet_nodes

      integer :: i, s, nz

#ifndef _THERMIQUE_
      call CommonMPI_abort("Jacobian_JacBigA_BigSm_dirichlet_nodes assumes thermal transfer is on.")
#endif

      do s = 1, NbNodeOwn_Ncpus(commRank + 1)
         if (IdNodeLocal(s)%T == "d") then
            call Jacobian_JacBigA_BigSm_find_diagonal_nz(JacBigA, s, s, nz)
            if (IdNodeLocal(s)%P == "d") then
               JacBigA%Val(:, :, JacBigA%Pt(s) + 1:JacBigA%Pt(s + 1)) = 0.d0
               do i = 1, NbCompThermique
                  JacBigA%Val(i, i, nz) = 1.d0
               end do
               BigSm(:, s) = 0.d0
            else
               JacBigA%Val(:, NbCompThermique, JacBigA%Pt(s) + 1:JacBigA%Pt(s + 1)) = 0.d0
               JacBigA%Val(NbCompThermique, NbCompThermique, nz) = 1.d0
               BigSm(NbCompThermique, s) = 0.d0
            end if
         else
            if (IdNodeLocal(s)%P == "d") &
               call CommonMPI_abort("Dirichlet condition on P only are not handled!")
         endif
      end do

   end subroutine Jacobian_JacBigA_BigSm_dirichlet_nodes

   !> \brief Sub-subroutine of Jacobian_JacBigA_BigSm for term n_k(X_j^n)
   subroutine Jacobian_JacBigA_BigSm_accumulation(Delta_t)

      double precision, intent(in) :: Delta_t

      integer :: k, rowk, colk, i, icp, nz, j, m, mph

      ! 2.1 div prim: n_k(X_j^n), k is node
      do k = 1, NbNodeOwn_Ncpus(commRank + 1) ! node

         rowk = k
         colk = k
         call Jacobian_JacBigA_BigSm_find_diagonal_nz(JacBigA, rowk, colk, nz)

         ! loop of components in Ctilde, n_k is an unknown independent
         do i = 1, NbCompCtilde_ctx(IncNode(k)%ic)
            icp = NumCompCtilde_ctx(i, IncNode(k)%ic)

            j = NbIncTotalPrim_ctx(IncNode(k)%ic) + i
            JacBigA%val(j, icp, nz) = VolDarcy%nodes(k)/Delta_t
         end do

         ! loop of component, n_k is not an unknown independent
         do m = 1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k, k is node
            mph = NumPhasePresente_ctx(m, IncNode(k)%ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  do j = 1, NbIncTotalPrim_ctx(IncNode(k)%ic)
                     JacBigA%Val(j, icp, nz) = JacBigA%Val(j, icp, nz) &
                                               + divDensiteMolaireSatCompNode(j, icp, mph, k)*PoroVolDarcy%nodes(k)/Delta_t
                  end do

                  bigSm(icp, rowk) = bigSm(icp, rowk) &
                                     - SmDensiteMolaireSatComp%nodes(icp, mph, k)*PoroVolDarcy%nodes(k)/Delta_t

               end if
            end do

         end do ! end of loop n_k

#ifdef _THERMIQUE_

         do m = 1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k, k is node
            mph = NumPhasePresente_ctx(m, IncNode(k)%ic)

            do j = 1, NbIncTotalPrim_ctx(IncNode(k)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                + PoroVolFourier%nodes(k)/Delta_t &
                                                *divDensiteMolaireEnergieInterneSatNode(j, mph, k)
            end do

            bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                      - PoroVolFourier%nodes(k)*SmDensiteMolaireEnergieInterneSatNode(mph, k)/Delta_t
         end do

         do j = 1, NbIncTotalPrim_ctx(IncNode(k)%ic)
            JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                             + Poro_1VolFourier%nodes(k)*CpRoche*divTemperatureNode(j, k)/Delta_t
         end do

         bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                   - Poro_1VolFourier%nodes(k)*CpRoche*SmTemperatureNode(k)/Delta_t

#endif

      end do ! end of div prim n_k, node

      ! 2.2 div prim: n_k(X_j^n), k is frac
      do k = 1, NbFracOwn_Ncpus(commRank + 1) ! node

         rowk = k + NbNodeOwn_Ncpus(commRank + 1) ! row of frac k in vector (node own, frac own, cell)
         colk = k + NbNodeLocal_Ncpus(commRank + 1)
         call Jacobian_JacBigA_BigSm_find_diagonal_nz(JacBigA, rowk, colk, nz)

         ! loop of composants in Ctilde, n_k is an unknown independent
         do i = 1, NbCompCtilde_ctx(IncFrac(k)%ic)
            icp = NumCompCtilde_ctx(i, IncFrac(k)%ic)

            j = NbIncTotalPrim_ctx(IncFrac(k)%ic) + i
            JacBigA%Val(j, icp, nz) = VolDarcy%fractures(k)/Delta_t
         end do

         ! loop of composants, n_k is not an unknown independent
         do m = 1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is node
            mph = NumPhasePresente_ctx(m, IncFrac(k)%ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                     JacBigA%Val(j, icp, nz) = JacBigA%Val(j, icp, nz) &
                                               + divDensiteMolaireSatCompFrac(j, icp, mph, k)*PoroVolDarcy%fractures(k)/Delta_t
                  end do

                  bigSm(icp, rowk) = bigSm(icp, rowk) &
                                     - SmDensiteMolaireSatComp%fractures(icp, mph, k)*PoroVolDarcy%fractures(k)/Delta_t
               end if
            end do

         end do ! end of loop phase

#ifdef _THERMIQUE_

         do m = 1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is frac
            mph = NumPhasePresente_ctx(m, IncFrac(k)%ic)

            do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = &
                  JacBigA%Val(j, NbComp + 1, nz) &
                  + PoroVolFourier%fractures(k)*divDensiteMolaireEnergieInterneSatFrac(j, mph, k)/Delta_t
            end do

            bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                      - PoroVolFourier%fractures(k)*SmDensiteMolaireEnergieInterneSatFrac(mph, k)/Delta_t
         end do

         do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
            JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                             + Poro_1VolFourier%fractures(k)*CpRoche*divTemperatureFrac(j, k)/Delta_t
         end do

         bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                   - Poro_1VolFourier%fractures(k)*CpRoche*SmTemperatureFrac(k)/Delta_t
#endif

      end do ! end of div prim n_k, frac

      ! 2.3 div prim: n_k(X_j^n), k is cell
      do k = 1, NbCellLocal_Ncpus(commRank + 1) ! cell

         rowk = k + NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) ! row of cell k in vector (node own, frac own, cell)
         colk = k + NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1)
         call Jacobian_JacBigA_BigSm_find_diagonal_nz(JacBigA, rowk, colk, nz)

         ! loop of composants in Ctilde, n_k is an unknown independent
         do i = 1, NbCompCtilde_ctx(IncCell(k)%ic)
            icp = NumCompCtilde_ctx(i, IncCell(k)%ic)

            j = NbIncTotalPrim_ctx(IncCell(k)%ic) + i
            JacBigA%Val(j, icp, nz) = VolDarcy%cells(k)/Delta_t
         end do

         ! loop of composants, n_k is not an unknown independent
         do m = 1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k, k is cell
            mph = NumPhasePresente_ctx(m, IncCell(k)%ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                     JacBigA%Val(j, icp, nz) = JacBigA%Val(j, icp, nz) &
                                               + divDensiteMolaireSatCompCell(j, icp, mph, k)*PoroVolDarcy%cells(k)/Delta_t
                  end do

                  bigSm(icp, rowk) = bigSm(icp, rowk) &
                                     - SmDensiteMolaireSatComp%cells(icp, mph, k)*PoroVolDarcy%cells(k)/Delta_t

               end if
            end do

         end do ! end of loop Q_k

#ifdef _THERMIQUE_

         do m = 1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k, k is cell
            mph = NumPhasePresente_ctx(m, IncCell(k)%ic)

            do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                + PoroVolFourier%cells(k)*divDensiteMolaireEnergieInterneSatCell(j, mph, k)/Delta_t
            end do

            bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                      - PoroVolFourier%cells(k)*SmDensiteMolaireEnergieInterneSatCell(mph, k)/Delta_t

         end do

         do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
            JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                             + Poro_1VolFourier%cells(k)*CpRoche*divTemperatureCell(j, k)/Delta_t
         end do

         bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                   - Poro_1VolFourier%cells(k)*CpRoche*SmTemperatureCell(k)/Delta_t
#endif

      end do ! end of div prim n_k, cell

   end subroutine Jacobian_JacBigA_BigSm_accumulation

   !> \brief Loop over the cells and node by cell to compute the Jacobian
   !! loop of cell, index is k
   !! {
   !!   1. loop of nodes s in cell k
   !!      {
   !!        1.1 conservation equations for cell control volumes
   !!            div( B * DarcyFlux_{k,s}), k is cell, s is node;
   !!               where B is DensiteMolaire*Kr/Visco or DensiteMolaire/Visco*Enthalpie
   !!            this term has three contributions to Jacobian:
   !!               A_kr (k is row cell, r is col cell k)
   !!               A_kr (k is row cell, r is col node s)
   !!               A_kr (k is row cell, r is col nodes/fracs in cell k)
   !!        1.2 conservation equations for (own) cells control volumes
   !!            A_sk, k is cell, s is node own;
   !!      }
   !!
   !!   2. loop of fracs in cell k
   !!      {
   !!        A_ks, k is cell, s is frac;
   !!        A_sk, k is cell, s is frac own;
   !!      }
   !! }
   subroutine Jacobian_JacBigA_BigSm_cell

      ! div prims and Sm from term div(DensiteMolaire*Kr/Visco*Comp)*FluxDarcy
      double precision :: &
         divK1(NbIncTotalPrimMax, NbComp), & ! k for cell, represent k in paper
         divS1(NbIncTotalPrimMax, NbComp), & ! s for node/frac, represent s in paper
         Sm1(NbComp)

      ! div prims and Sm from term DensiteMolaire*Kr/Visco*Comp*div(FluxDarcy)
      ! three contributions from this item
      double precision :: &
         divK2(NbIncTotalPrimMax, NbComp), & ! A_kr, k is row cell, r is col cell k
         divS2(NbIncTotalPrimMax, NbComp), & ! A_kr, k is row cell, r is col col node s
         divR2(NbIncTotalPrimMax, NbComp, NbNodeCellMax + NbFracCellMax), & ! A_kr, k is row cell, r is nodes/fracs in cell k
         Sm2(NbComp)

      ! div prim of Darcy flux
      double precision :: &
         divrho_k(NbIncTotalPrimMax, NbPhase), &
         divrho_s(NbIncTotalPrimMax, NbPhase), &
         Smrho_k(NbPhase), &
         Smrho_s(NbPhase), &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), & !
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r(NbIncTotalPrimMax, NbPhase, NbNodeCellMax + NbFracCellMax), &
         SmDarcyFlux(NbPhase)

#ifdef _THERMIQUE_
      ! div prim of Fourier flux
      double precision :: &
         divFourierFlux_k(NbIncTotalPrimMax), &
         divFourierFlux_r(NbIncTotalPrimMax, NbNodeCellMax + NbFracCellMax), & ! r represent s' in paper
         SmFourierFlux

      ! div prims and Sm from term div(DensiteMolaire*Kr*Enthalpie/Visco*FluxDarcy)
      double precision :: &
         divEgK(NbIncTotalPrimMax), & ! K for cell, represent k in paper
         divEgS(NbIncTotalPrimMax), & ! S for node/frac, represent s in paper
         divEgR(NbIncTotalPrimMax, NbNodeCellMax + NbFracCellMax), & ! R for node/frac, represent s' in paper
         SmEg
#endif

      integer :: k, s, nums, sf, r, numr, rf, i, j
      integer :: rows, cols, colr, nz
      integer :: m
      integer :: NbNodeCell, NbFracCell

      integer :: rowK, colk, &                    ! row/col of cell k in JacBigA
                 rowSR(NbNodeCellMax + NbFracCellMax), & ! rows (in JacBigA) of nodes/frac in cell k
                 colSR(NbNodeCellMax + NbFracCellMax)    ! cols (in JacBigA) of nodes/frac in cell k

      ! main loop of cell
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         nbNodeCell = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         nbFracCell = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)

         ! (rowk, colk) is row/col of cell k in JacBigA
         ! (rowSR, colSR) are rows/cols (in JacBigA) of nodes/frac connected to cell k
         ! computes mapping between local connectivity and the position of connected dofs in the Jacobian matrix
         call Jacobian_RowCol_KSR(k, nbNodeCell, nbFracCell, &
                                  rowk, colk, rowSR, colSR)

         ! csrK(:) = 0
         do m = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
            csrK(JacBigA%Num(m)) = m - JacBigA%Pt(rowk)
         end do

         ! two parts: 1. s is node; 2. s is frac

         ! 1. s is node
         do s = 1, NbNodeCell
            nums = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + s) ! num node of s

            ! compute \sum{P_i \cap Q_{k/s} } div(DensiteMolaire*Kr/Visco) * DarcyFlux : k is cell, s is node
            ! this term has two contributions: A_kr, r is k -> divK1
            !                                  A_kr, r is s -> divS1
            call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux( &
               NbIncTotalPrim_ctx, &
               IncCell(k)%ic, IncNode(nums)%ic, FluxDarcyKI(:, s, k), &
               divDensiteMolaireKrViscoCompCell(:, :, :, k), SmDensiteMolaireKrViscoCompCell(:, :, k), &
               divDensiteMolaireKrViscoCompNode(:, :, :, nums), SmDensiteMolaireKrViscoCompNode(:, :, nums), &
               divK1, divS1, Sm1)

            ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
            call Jacobian_divrho_gravity( &
               IncCell(k), divSaturationCell(:, :, k), &
               DensiteMassiqueCell(:, k), divDensiteMassiqueCell(:, :, k), SmDensiteMassiqueCell(:, k), &
               divrho_k, Smrho_k, &
               IncNode(nums), divSaturationNode(:, :, nums), &
               DensiteMassiqueNode(:, nums), divDensiteMassiqueNode(:, :, nums), SmDensiteMassiqueNode(:, nums), &
               divrho_s, Smrho_s)

            ! compute div Darcy flux
            call Jacobian_divDarcyFlux( &
               NbIncTotalPrim_ctx, phase_can_be_present, & ! does not depend on k or s
               IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncNode(nums)%ic, XCellLocal(3, k), &
               XNodeLocal(3, :), XFaceLocal(3, :), gravity, & ! does not depend on k or s
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               TkLocal_Darcy(k)%pt(s, :), &
               divPressionCell(:, k), divPhasePressureCell(:, :, k), divrho_k, &
               SmPressionCell(k), Smrho_k, &
               divrho_s, Smrho_s, &
               divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
               divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

            ! compute DensiteMolaire*Kr/Visco * div(Flux) using div(DarcyFlux)
            call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncNode(nums)%ic, FluxDarcyKI(:, s, k), &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               DensiteMolaireKrViscoCompCell(:, :, k), &
               DensiteMolaireKrViscoCompNode(:, :, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

            ! compute div (DensiteMolaire*Enthalpie/Visco * DarcyFlux) using div(DarcyFlux),
            ! k is cell, s is node
            call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncNode(nums)%ic, FluxDarcyKI(:, s, k), &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               DensiteMolaireKrViscoEnthalpieCell(:, k), &
               divDensiteMolaireKrViscoEnthalpieCell(:, :, k), &
               SmDensiteMolaireKrViscoEnthalpieCell(:, k), &
               DensiteMolaireKrViscoEnthalpieNode(:, nums), &
               divDensiteMolaireKrViscoEnthalpieNode(:, :, nums), &
               SmDensiteMolaireKrViscoEnthalpieNode(:, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

            ! compute div FluxFourier
            call Jacobian_divFourierFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               TkLocal_Fourier(k)%pt(s, :), &
               divTemperatureCell(:, k), SmTemperatureCell(k), &
               divTemperatureNode, SmTemperatureNode, &
               divTemperatureFrac, SmTemperatureFrac, &
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

            ! divK, divS, divR to JacBigA
            ! row of cell k then row of node s

            ! A_kk
            nz = JacBigA%Pt(rowk) + csrK(colk)
            do i = 1, NbComp
               do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                  JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divK1(j, i) + divK2(j, i)
               end do
            end do

#ifdef _THERMIQUE_

            do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                + divEgK(j) + divFourierFlux_k(j)
            end do
#endif

            ! A_ks
            cols = colSR(s)
            nz = JacBigA%Pt(rowk) + csrK(cols)
            do i = 1, NbComp
               do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divS1(j, i) + divS2(j, i)
               end do
            end do

#ifdef _THERMIQUE_

            ! ps. divFourierFlux_s = 0
            do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) + divEgS(j)
            end do
#endif

            ! A_kr, r is node
            do r = 1, NbNodeCell ! r represent s' in paper, r is node
               numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + r)

               colr = colSR(r) ! col
               nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divR2(j, i, r)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   + divEgR(j, r) + divFourierFlux_r(j, r)
               end do
#endif
            end do

            ! A_kr, r is frac
            do r = 1, NbFracCell ! r represent s' in paper, r is frac
               numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + r) ! numr is frac num
               rf = r + NbNodeCell

               colr = colSR(rf)
               nz = JacBigA%Pt(rowk) + csrK(colr)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divR2(j, i, rf)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   + divEgR(j, rf) + divFourierFlux_r(j, rf)
               end do
#endif
            end do

            ! Sm
            bigSm(1:NbComp, rowk) = bigSm(1:NbComp, rowk) &
                                    - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_

            bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                      - SmEg - SmFourierFlux
#endif

            ! we consider the jacobian rows corresponding to node s when node s is own
            if (IdNodeLocal(nums)%Proc == "o") then

               ! csrSR(:) = 0
               rows = rowSR(s)
               do m = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                  csrSR(JacBigA%Num(m)) = m - JacBigA%Pt(rows)
               end do

               ! A_sk, s is node, k is cell
               nz = JacBigA%Pt(rows) + csrSR(colk)

               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divK1(j, i) - divK2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_
               do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   - divEgK(j) - divFourierFlux_k(j)
               end do
#endif

               ! A_ss, s is node
               cols = colSR(s)
               nz = JacBigA%Pt(rows) + csrSR(cols)

               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) &
                                             - divS1(j, i) - divS2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_
               ! ps. divFourierFlux_s=0
               do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   - divEgS(j)
               end do
#endif

               ! A_sr, s is node, r is node
               do r = 1, NbNodeCell ! r represent s' in paper
                  numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + r)

                  colr = colSR(r)
                  nz = JacBigA%Pt(rows) + csrSR(colr)

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divR2(j, i, r)
                     end do
                  end do

#ifdef _THERMIQUE_
                  do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgR(j, r) - divFourierFlux_r(j, r)
                  end do
#endif
               end do

               ! A_sr, s is node, r is frac
               do r = 1, NbFracCell ! r represent s' in paper
                  numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + r)
                  rf = r + NbNodeCell

                  colr = colSR(rf)
                  nz = JacBigA%Pt(rows) + csrSR(colr)

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divR2(j, i, rf)
                     end do
                  end do

#ifdef _THERMIQUE_
                  do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgR(j, rf) - divFourierFlux_r(j, rf)
                  end do
#endif
               end do

               bigSm(1:NbComp, rows) = bigSm(1:NbComp, rows) &
                                       + Sm1(1:NbComp) + Sm2(1:NbComp)

#ifdef _THERMIQUE_
               bigSm(NbComp + 1, rows) = bigSm(NbComp + 1, rows) &
                                         + SmEg + SmFourierFlux
#endif

            end if

         end do ! end of row s

         ! r (s' in paper) will be taken care of "later" by the same loop when dealing with s

         ! 2. s is frac
         do s = 1, NbFracCell
            nums = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + s) ! nums is frac num

            sf = s + NbNodeCell

            ! compute div(DensiteMolaire*Kr/Visco) * DarcyFlux when k is cell, s is frac
            call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux( &
               NbIncTotalPrim_ctx, &
               IncCell(k)%ic, IncFrac(nums)%ic, FluxDarcyKI(:, sf, k), &
               divDensiteMolaireKrViscoCompCell(:, :, :, k), SmDensiteMolaireKrViscoCompCell(:, :, k), &
               divDensiteMolaireKrViscoCompFrac(:, :, :, nums), SmDensiteMolaireKrViscoCompFrac(:, :, nums), &
               divK1, divS1, Sm1)

            ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
            call Jacobian_divrho_gravity( & ! cell <-> frac
               IncCell(k), divSaturationCell(:, :, k), &
               DensiteMassiqueCell(:, k), divDensiteMassiqueCell(:, :, k), SmDensiteMassiqueCell(:, k), &
               divrho_k, Smrho_k, &
               IncFrac(nums), divSaturationFrac(:, :, nums), &
               DensiteMassiqueFrac(:, nums), divDensiteMassiqueFrac(:, :, nums), SmDensiteMassiqueFrac(:, nums), &
               divrho_s, Smrho_s)

            ! compute div Darcy flux
            call Jacobian_divDarcyFlux( &
               NbIncTotalPrim_ctx, phase_can_be_present, & ! does not depend on k or s
               IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncFrac(nums)%ic, XCellLocal(3, k), &
               XNodeLocal(3, :), XFaceLocal(3, :), gravity, & ! does not depend on k or s
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               TkLocal_Darcy(k)%pt(sf, :), &
               divPressionCell(:, k), divPhasePressureCell(:, :, k), divrho_k, &
               SmPressionCell(k), Smrho_k, &
               divrho_s, Smrho_s, &
               divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
               divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

            ! compute DensiteMolaire*Kr/Visco * div(Flux) using div(DarcyFlux)
            call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncFrac(nums)%ic, FluxDarcyKI(:, sf, k), &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               DensiteMolaireKrViscoCompCell(:, :, k), &
               DensiteMolaireKrViscoCompFrac(:, :, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

            ! compute div (DensiteMolaire*Enthalpie/Visco * DarcyFlux) using div(DarcyFlux)
            call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, IncFrac(nums)%ic, FluxDarcyKI(:, sf, k), &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               DensiteMolaireKrViscoEnthalpieCell(:, k), &
               divDensiteMolaireKrViscoEnthalpieCell(:, :, k), &
               SmDensiteMolaireKrViscoEnthalpieCell(:, k), &
               DensiteMolaireKrViscoEnthalpieFrac(:, nums), &
               divDensiteMolaireKrViscoEnthalpieFrac(:, :, nums), &
               SmDensiteMolaireKrViscoEnthalpieFrac(:, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

            ! compute div FluxFourier
            call Jacobian_divFourierFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncCell(k)%ic, &
               NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + 1:NodebyCellLocal%Pt(k + 1)), &
               FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + 1:FracbyCellLocal%Pt(k + 1)), &
               TkLocal_Fourier(k)%pt(sf, :), &
               divTemperatureCell(:, k), SmTemperatureCell(k), &
               divTemperatureNode, SmTemperatureNode, & ! does not depend on k or s
               divTemperatureFrac, SmTemperatureFrac, & ! does not depend on k or s
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

            ! divK, divS, divR to JacBigA
            ! row of cell k then row of frac s (nums)

            ! A_kk
            nz = JacBigA%Pt(rowk) + csrK(colk)
            do i = 1, NbComp
               do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                  JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divK1(j, i) + divK2(j, i)
               end do
            end do

#ifdef _THERMIQUE_

            do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                + divEgK(j) + divFourierFlux_k(j)
            end do
#endif

            ! A_ks, s is frac
            cols = colSR(sf)
            nz = JacBigA%Pt(rowk) + csrK(cols)
            do i = 1, NbComp
               do j = 1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                  JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divS1(j, i) + divS2(j, i)
               end do
            end do

#ifdef _THERMIQUE_

            ! ps. divFourierFlux_s = 0
            do j = 1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
               JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) + divEgS(j)
            end do
#endif

            ! A_kr, r is node
            do r = 1, NbNodeCell ! r represent s' in paper, r is node
               numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + r)

               colr = colSR(r) ! col
               nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divR2(j, i, r)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   + divEgR(j, r) + divFourierFlux_r(j, r)
               end do
#endif
            end do

            ! A_kr, r is frac
            do r = 1, NbFracCell ! r represent s' in paper, r is frac
               numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + r) ! numr is frac num
               rf = r + NbNodeCell

               colr = colSR(rf)
               nz = JacBigA%Pt(rowk) + csrK(colr)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divR2(j, i, rf)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   + divEgR(j, rf) + divFourierFlux_r(j, rf)
               end do
#endif
            end do

            ! Sm
            bigSm(1:NbComp, rowk) = bigSm(1:NbComp, rowk) &
                                    - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_
            bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                      - SmEg - SmFourierFlux
#endif

            ! frac own
            ! its corresponding face num <= NbFaceOwn
            if (FracToFaceLocal(nums) <= NbFaceOwn_Ncpus(commRank + 1)) then

               ! csrSR(:) = 0
               rows = rowSR(sf)

               do m = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                  csrSR(JacBigA%Num(m)) = m - JacBigA%Pt(rows)
               end do

               ! A_sk, s is frac (own), k is cell
               nz = JacBigA%Pt(rows) + csrSR(colk)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divK1(j, i) - divK2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   - divEgK(j) - divFourierFlux_k(j)
               end do
#endif

               ! A_ss, s is frac (own)
               cols = colSR(sf)
               nz = JacBigA%Pt(rows) + csrSR(cols)

               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divS1(j, i) - divS2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_

               ! ps. divFourierFlux_s=0
               do j = 1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   - divEgS(j)
               end do
#endif

               ! A_sr, s is frac (own), r is node
               do r = 1, NbNodeCell ! r represent s' in paper
                  numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + r)

                  colr = colSR(r)
                  nz = JacBigA%Pt(rows) + csrSR(colr)

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divR2(j, i, r)
                     end do
                  end do

#ifdef _THERMIQUE_
                  do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgR(j, r) - divFourierFlux_r(j, r)
                  end do
#endif
               end do

               ! A_sr, s is frac (own), r is frac
               do r = 1, NbFracCell ! r represent s' in paper
                  numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + r)
                  rf = r + NbNodeCell

                  colr = colSR(rf)
                  nz = JacBigA%Pt(rows) + csrSR(colr)

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divR2(j, i, rf)
                     end do
                  end do

#ifdef _THERMIQUE_

                  do j = 1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgR(j, rf) - divFourierFlux_r(j, rf)
                  end do
#endif
               end do

               ! Sm
               bigSm(1:NbComp, rows) = bigSm(1:NbComp, rows) &
                                       + Sm1(1:NbComp) + Sm2(1:NbComp)

#ifdef _THERMIQUE_

               bigSm(NbComp + 1, rows) = bigSm(NbComp + 1, rows) &
                                         + SmEg + SmFourierFlux
#endif

            end if ! end of frac won

         end do ! end of s frac

      end do ! end of cell k

   end subroutine Jacobian_JacBigA_BigSm_cell

   subroutine Jacobian_JacBigA_locate_frac_row(row)

      integer, intent(in) :: row

      integer :: rowk, colk, &      ! row/col of frac k in JacBigA
                 rowSR(NbNodeFaceMax), & ! rows (in JacBigA) of nodes in frac k
                 colSR(NbNodeFaceMax)    ! cols (in JacBigA) of nodes in frac k

      integer :: k, fk, i
      integer :: nbNodeFrac

      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         fk = FracToFaceLocal(k) ! fk is face num
         nbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)

         call Jacobian_RowCol_FR(k, nbNodeFrac, &
                                 rowk, colk, rowSR, colSR)

         if (row >= rowk .and. row < rowk + NbCompThermique) then
            print *, "row:", row, "frac:", k, "face:", fk
         end if

         do i = 1, NbNodeFaceMax
            if (row >= rowSR(i) .and. row < rowSR(i) + NbCompThermique) then
               print *, "row:", row, "frac:", k, "face:", fk, "node", i, rowSR
            end if
         end do

      end do

   end subroutine Jacobian_JacBigA_locate_frac_row

   ! loop of frac, index is k
   ! {
   !   loop of nodes in frac k
   !   {
   !     A_ks, k is frac, s is node; (rowk/colk is row/col of cell k)
   !     A_sk, k is frac, s is node own;
   !   }
   ! }
   subroutine Jacobian_JacBigA_BigSm_frac

      ! div prims and Sm from term div(DensiteMolaire*Kr/Visco*Comp)*FluxDarcy
      double precision :: &
         divK1(NbIncTotalPrimMax, NbComp), & ! K for frac, represent k in paper
         divS1(NbIncTotalPrimMax, NbComp), & ! S for node, represent s in paper
         Sm1(NbComp)

      ! div prims and Sm from term DensiteMolaire*Kr/Visco*Comp*div(FluxDarcy)
      double precision :: &
         divK2(NbIncTotalPrimMax, NbComp), &
         divS2(NbIncTotalPrimMax, NbComp), &
         divR2(NbIncTotalPrimMax, NbComp, NbNodeFaceMax), & ! R for node/frac, represent s' in paper
         Sm2(NbComp)

      double precision :: &
         divrho_k(NbIncTotalPrimMax, NbPhase), &
         divrho_s(NbIncTotalPrimMax, NbPhase), &
         Smrho_k(NbPhase), &
         Smrho_s(NbPhase), &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r(NbIncTotalPrimMax, NbPhase, NbNodeFaceMax), & ! r represen
         SmDarcyFlux(NbPhase)

#ifdef _THERMIQUE_

      ! div FluxFourier
      double precision :: &
         divFourierFlux_k(NbIncTotalPrimMax), &
         !  divFourierFlux_s(NbIncTotalPrimMax), &
         divFourierFlux_r(NbIncTotalPrimMax, NbNodeFaceMax), & ! r represent s' in paper
         SmFourierFlux

      ! div prims and Sm from term div(DensiteMolaire*Kr*Enthalpie/Visco*FluxDarcy)
      double precision :: &
         divEgK(NbIncTotalPrimMax), & ! K for frac, represent k in paper
         divEgS(NbIncTotalPrimMax), & ! S for node, represent s in paper
         divEgR(NbIncTotalPrimMax, NbNodeFaceMax), & ! R for node/frac, represent s' in paper
         SmEg
#endif

      integer :: rowK, colk, &      ! row/col of frac k in JacBigA
                 rowSR(NbNodeFaceMax), & ! rows (in JacBigA) of nodes in frac k
                 colSR(NbNodeFaceMax)    ! cols (in JacBigA) of nodes in frac k

      integer :: k, s, nums, r, numr, fk, i, j
      integer :: rows, cols, colr, nz
      integer :: m
      integer :: nbNodeFrac
      integer(c_int), dimension(0) :: empty_array

      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         fk = FracToFaceLocal(k) ! fk is face num
         nbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)

         call Jacobian_RowCol_FR(k, nbNodeFrac, &
                                 rowk, colk, rowSR, colSR)

         ! csrK(:) = 0
         do m = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
            csrK(JacBigA%Num(m)) = m - JacBigA%Pt(rowk) ! CHECKME: m: 1 -> NbCompThermique ?
         end do

         do s = 1, nbNodeFrac
            nums = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + s)

            ! compute div(DensiteMolaire*Kr/Visco) * DarcyFlux : k is frac, s is node
            call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux( &
               NbIncTotalPrim_ctx, &
               IncFrac(k)%ic, IncNode(nums)%ic, FluxDarcyFI(:, s, k), &
               divDensiteMolaireKrViscoCompFrac(:, :, :, k), SmDensiteMolaireKrViscoCompFrac(:, :, k), &
               divDensiteMolaireKrViscoCompNode(:, :, :, nums), SmDensiteMolaireKrViscoCompNode(:, :, nums), &
               divK1, divS1, Sm1)

            ! compute div Darcy flux
            ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
            call Jacobian_divrho_gravity( & ! frac <-> node
               IncFrac(k), divSaturationFrac(:, :, k), &
               DensiteMassiqueFrac(:, k), divDensiteMassiqueFrac(:, :, k), SmDensiteMassiqueFrac(:, k), &
               divrho_k, Smrho_k, &
               IncNode(nums), divSaturationNode(:, :, nums), &
               DensiteMassiqueNode(:, nums), divDensiteMassiqueNode(:, :, nums), SmDensiteMassiqueNode(:, nums), &
               divrho_s, Smrho_s)

            call Jacobian_divDarcyFlux( &
               NbIncTotalPrim_ctx, phase_can_be_present, & ! does not depend on k or s
               IncNode, IncFrac, & ! does not depend on k or s
               IncFrac(k)%ic, IncNode(nums)%ic, XFaceLocal(3, fk), &
               XNodeLocal(3, :), XFaceLocal(3, :), gravity, & ! does not depend on k or s
               NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + 1:NodebyFaceLocal%Pt(fk + 1)), &
               empty_array, &
               TkFracLocal_Darcy(k)%pt(s, :), &
               divPressionFrac(:, k), divPhasePressureFrac(:, :, k), divrho_k, &
               SmPressionFrac(k), Smrho_k, &
               divrho_s, Smrho_s, &
               divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
               divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

            ! compute DensiteMolaire*Kr/Visco * div(Flux) using div(DarcyFlux)
            call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncFrac(k)%ic, IncNode(nums)%ic, FluxDarcyFI(:, s, k), &
               NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + 1:NodebyFaceLocal%Pt(fk + 1)), &
               empty_array, &
               DensiteMolaireKrViscoCompFrac(:, :, k), &
               DensiteMolaireKrViscoCompNode(:, :, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

            ! compute div (DensiteMolaire*Enthalpie/Visco * DarcyFlux) using div(DarcyFlux)
            ! reuse what's been computed by Jacobian_divDarcyFlux frac/node
            call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncFrac(k)%ic, IncNode(nums)%ic, FluxDarcyFI(:, s, k), &
               NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + 1:NodebyFaceLocal%Pt(fk + 1)), &
               empty_array, &
               DensiteMolaireKrViscoEnthalpieFrac(:, k), &
               divDensiteMolaireKrViscoEnthalpieFrac(:, :, k), &
               SmDensiteMolaireKrViscoEnthalpieFrac(:, k), &
               DensiteMolaireKrViscoEnthalpieNode(:, nums), &
               divDensiteMolaireKrViscoEnthalpieNode(:, :, nums), &
               SmDensiteMolaireKrViscoEnthalpieNode(:, nums), &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

            ! compute div FluxFourier
            call Jacobian_divFourierFlux( &
               NbIncTotalPrim_ctx, IncNode, IncFrac, & ! does not depend on k or s
               IncFrac(k)%ic, &
               NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + 1:NodebyFaceLocal%Pt(fk + 1)), &
               empty_array, &
               TkFracLocal_Fourier(k)%pt(s, :), &
               divTemperatureFrac(:, k), SmTemperatureFrac(k), &
               divTemperatureNode, SmTemperatureNode, & ! does not depend on k or s
               divTemperatureFrac, SmTemperatureFrac, & ! does not depend on k or s
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

            ! row with frac own
            if (k <= NbFracOwn_Ncpus(commRank + 1)) then

               ! A_kk
               nz = JacBigA%Pt(rowk) + csrK(colk)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divK1(j, i) + divK2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                   + divEgK(j) + divFourierFlux_k(j)
               end do
#endif

               ! A_ks
               cols = colSR(s)
               nz = JacBigA%Pt(rowk) + csrK(cols)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divS1(j, i) + divS2(j, i)
                  end do
               end do

#ifdef _THERMIQUE_

               ! ps. divFourierFlux_s = 0
               do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) + divEgS(j)
               end do
#endif

               ! A_ks'
               do r = 1, NbNodeFrac ! r represent s' in paper, r is node

                  numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + r)

                  colr = colSR(r) ! col
                  nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divR2(j, i, r)
                     end do
                  end do

#ifdef _THERMIQUE_

                  do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      + divEgR(j, r) + divFourierFlux_r(j, r)
                  end do
#endif
               end do

               ! Sm
               bigSm(1:NbComp, rowk) = bigSm(1:NbComp, rowk) &
                                       - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_

               bigSm(NbComp + 1, rowk) = bigSm(NbComp + 1, rowk) &
                                         - SmEg - SmFourierFlux
#endif

            end if

            ! node s is own
            if (IdNodeLocal(nums)%Proc == "o") then

               ! csrSR(:) = 0
               rows = rowSR(s)
               do m = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                  csrSR(JacBigA%Num(m)) = m - JacBigA%Pt(rows)
               end do

               ! A_sk
               nz = JacBigA%Pt(rows) + csrSR(colk)

               if (IdNodeLocal(nums)%P /= "d") then

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divK1(j, i) - divK2(j, i)
                     end do
                  end do
               end if

#ifdef _THERMIQUE_

               if (IdNodeLocal(nums)%T /= "d") then

                  do j = 1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgK(j) - divFourierFlux_k(j)
                  end do
               end if
#endif

               ! A_ss
               cols = colSR(s)
               nz = JacBigA%Pt(rows) + csrSR(cols)

               if (IdNodeLocal(nums)%P /= "d") then

                  do i = 1, NbComp
                     do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                        JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divS1(j, i) - divS2(j, i)
                     end do
                  end do

               end if

#ifdef _THERMIQUE_

               if (IdNodeLocal(nums)%T /= "d") then

                  ! ps. divFourierFlux_s=0
                  do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                     JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                      - divEgS(j)
                  end do
               end if
#endif

               ! A_ss', s' is node
               do r = 1, NbNodeFrac ! r represent s' in paper
                  numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + r)

                  colr = colSR(r)
                  nz = JacBigA%Pt(rows) + csrSR(colr)

                  if (IdNodeLocal(nums)%P /= "d") then
                     do i = 1, NbComp
                        do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                           JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) - divR2(j, i, r)
                        end do
                     end do
                  end if

#ifdef _THERMIQUE_

                  if (IdNodeLocal(nums)%T /= "d") then
                     do j = 1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                        JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) &
                                                         - divEgR(j, r) - divFourierFlux_r(j, r)
                     end do
                  end if
#endif
               end do

               ! Sm
               if (IdNodeLocal(nums)%P /= "d") then
                  bigSm(1:NbComp, rows) = bigSm(1:NbComp, rows) &
                                          + Sm1(1:NbComp) + Sm2(1:NbComp)
               end if

#ifdef _THERMIQUE_

               if (IdNodeLocal(nums)%T /= "d") then
                  bigSm(NbComp + 1, rows) = bigSm(NbComp + 1, rows) &
                                            + SmEg + SmFourierFlux
               end if
#endif

            end if

         end do ! s in frac k
      end do ! frac k

   end subroutine Jacobian_JacBigA_BigSm_frac

   ! loop of injection well
   ! (qw - qmol) * (Pwmax - Pw) = - qmol * (Pwmax-Pw) + qw * (Pwmax - Pw)
   !   where qw = sum_{s} sum_{i} q_{w,s,i}
   subroutine Jacobian_JacBigA_BigSm_wellinj

      integer :: k, rowk, colk, s, nums, rows, cols, m, icp, nz
      double precision :: Tws, Ps_Pws, Ts, WIDws, WIFws
      real(c_double) :: dP_w(NbComp), dP_s(NbIncTotalPrimMax, NbComp) ! NbComp mass balance equation
      real(c_double) :: dP_ER_w, der_ER_s(NbIncTotalPrimMax) ! energy balance equation
      logical :: something_is_injected, well_is_closed

      nz = -1

      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)

         something_is_injected = .false.
         ! CHECKME: we could directly jump to regularization
         well_is_closed = (DataWellInjLocal(k)%IndWell == 'c')

         ! A_kk, k is well
         rowk = k + NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
                + NbCellLocal_Ncpus(commRank + 1)
         colk = k + NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1) &
                + NbCellLocal_Ncpus(commRank + 1)

         ! assembly equations of own wells
         if (k <= NbWellInjOwn_Ncpus(commRank + 1)) then
            do m = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
               csrK(JacBigA%Num(m)) = m - JacBigA%Pt(rowk)
            end do
         end if

         ! Pwmax_Pw = DataWellInjLocal(k)%PressionMax - IncPressionWellInj(k) ! P_w^{max}-P_w

         ! Two cases:
         !   if indwell=="p": well equation is Pwmax-Pw = 0
         !   if indwell=="f": well equation is qw-qmol = 0, where qw = sum_{s} sum_{i} q_{w,s,i}

         ! Two steps:
         !   Step 1. well equation if indwell=="p"
         !   Step 2. loop of nodes of well k
         !           - well equation if indwell=="f"
         !           - reservoir equation

         ! Step 1. well equation: Pwmax - Pw = 0
         if (DataWellInjLocal(k)%IndWell == 'p') then

            if (k <= NbWellInjOwn_Ncpus(commRank + 1)) then ! own injection well
               nz = JacBigA%Pt(rowk) + csrK(colk)
               JacBigA%Val(1, 1, nz) = -1.d0
            end if
            something_is_injected = .true.
         end if

         ! Step 2. nodes in well k
         do s = NodebyWellInjLocal%Pt(k) + 1, NodebyWellInjLocal%Pt(k + 1)
            nums = NodebyWellInjLocal%Num(s) ! num node

            Ps_Pws = IncNode(nums)%phase_pressure(LIQUID_PHASE) - PerfoWellInj(s)%Pression ! P_s - P_{w,s}
            Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
            Ts = IncNode(nums)%Temperature    ! T_s

            WIDws = NodeDatabyWellInjLocal%Val(s)%WID
            WIFws = NodeDatabyWellInjLocal%Val(s)%WIF

            dP_w = 0.d0
            dP_s = 0.d0
#ifdef _THERMIQUE_
            dP_ER_w = 0.d0
            der_ER_s = 0.d0
#endif

            if (IdNodeLocal(nums)%Proc == "o") then ! ps. this node can be in the boundary
               rows = nums
               do m = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                  csrSR(JacBigA%Num(m)) = m - JacBigA%Pt(rows)
               end do
            end if

            cols = nums
            if ((Ps_Pws < 0.d0) .AND. (.NOT. well_is_closed)) then ! if >0, this term is 0

               something_is_injected = .true.

               do icp = 1, NbComp
                  dP_w(icp) = divDensiteMolaireKrViscoCompWellInj(icp, s)*WIDws*Ps_Pws &
                              - DensiteMolaireKrViscoCompWellInj(icp, s)*WIDws
                  dP_s(:, icp) = DensiteMolaireKrViscoCompWellInj(icp, s)*WIDws &
                                 *(divPressionNode(:, nums) + divPhasePressureNode(:, LIQUID_PHASE, nums))
               end do

#ifdef _THERMIQUE_
               dP_ER_w = divDensiteMolaireKrViscoEnthalpieWellInj(s)*WIDws*Ps_Pws &
                         - DensiteMolaireKrViscoEnthalpieWellInj(s)*WIDws
               der_ER_s = DensiteMolaireKrViscoEnthalpieWellInj(s)*WIDws &
                          *(divPressionNode(:, nums) + divPhasePressureNode(:, LIQUID_PHASE, nums))
#endif

               if (k <= NbWellInjOwn_Ncpus(commRank + 1)) then ! own injection well

                  if (DataWellInjLocal(k)%IndWell == 'f') then

                     ! A_kk, k is own injection well
                     nz = JacBigA%Pt(rowk) + csrK(colk)
                     do icp = 1, NbComp
                        JacBigA%Val(1, 1, nz) = JacBigA%Val(1, 1, nz) + dP_w(icp)
                     end do

                     ! A_ks, k is own injection well, s is node
                     nz = JacBigA%Pt(rowk) + csrK(cols)
                     do icp = 1, NbComp
                        JacBigA%Val(:, 1, nz) = JacBigA%Val(:, 1, nz) + dP_s(:, icp)
                     end do
                  end if
               end if

               if (IdNodeLocal(nums)%Proc == "o") then ! ps. this node can be in the boundary

                  ! Ask, s is node, k is injection well
                  nz = JacBigA%Pt(rows) + csrSR(colk)
                  JacBigA%Val(1, 1:NbComp, nz) = JacBigA%Val(1, 1:NbComp, nz) + dP_w(:) ! term q_{w,s,i}, derivative of P_w

#ifdef _THERMIQUE_
                  JacBigA%Val(1, NbComp + 1, nz) = JacBigA%Val(1, NbComp + 1, nz) + dP_ER_w
#endif

                  ! Ass, s is node
                  nz = JacBigA%Pt(rows) + csrSR(cols)
                  JacBigA%Val(:, 1:NbComp, nz) = JacBigA%Val(:, 1:NbComp, nz) + dP_s(:, 1:NbComp) ! term q_{w,s,i}, derivative of P_w

#ifdef _THERMIQUE_
                  JacBigA%Val(:, NbComp + 1, nz) = JacBigA%Val(:, NbComp + 1, nz) + der_ER_s
#endif
               end if

            end if

!           ! (3) WIF_{w,s} (T_s - T_{w,s})
! #ifdef _THERMIQUE_
!           if(IdNodeLocal(nums)%Proc=="o") then ! ps. this node can be in the boundary

!              ! Ass, s is node
!              nz = JacBigA%Pt(rows) + csrSR(cols)
!              JacBigA%Val(:,NbComp+1,nz) = JacBigA%Val(:,NbComp+1,nz) + divTemperatureNode(:,nums)

!              bigSm(NbComp+1,nums) = bigSm(NbComp+1,nums) - WIFws * SmTemperatureNode(nums)
!           end if
! #endif
         end do

         !Make sure Jacobian is not singular
         if (well_is_closed .OR. ((DataWellInjLocal(k)%IndWell == 'f') .AND. (.NOT. something_is_injected))) then

            if (k <= NbWellInjOwn_Ncpus(commRank + 1)) then ! own injection well
               nz = JacBigA%Pt(rowk) + csrK(colk)
               JacBigA%Val(1, 1, nz) = 1.d0
            end if

#ifdef COMPASS_LOG_WELL_INFO
            write (*, *) '[Well Monitoring]  Jacobian regularization has been performed, since nothing is injected for well', &
               k, 'on proc', commRank + 1
#endif
         end if

      end do !k-well loop

   end subroutine Jacobian_JacBigA_BigSm_wellinj

   ! loop of production well
   ! (qmol - qw) * (Pw - Pwmin) = qmol * (Pw-Pwmin) - qw * (Pw - Pwmin)
   ! where qw = sum_{s} sum_{i} q_{w,s,i}
   subroutine Jacobian_JacBigA_BigSm_wellprod

      integer :: k, rowk, colk, s, nums, rows, cols, m, mph, n, icp, nz
      double precision :: Pws, Ps, WIDws, WIFws, Ps_Pws
      real(c_double) :: dP_w(NbComp), dP_s(NbIncTotalPrimMax, NbComp) ! NbComp mass balance equation
      real(c_double) :: dP_ER_w, der_ER_s(NbIncTotalPrimMax) ! energy balance equation
      logical :: something_is_produced, well_is_own, well_is_closed, well_node_is_active

      nz = -1

      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)
         something_is_produced = .false.
         well_is_own = (k <= NbWellProdOwn_Ncpus(commRank + 1))
         ! CHECKME: we could directly jump to regularization
         well_is_closed = (DataWellProdLocal(k)%IndWell == 'c')

         ! A_kk, k is well
         rowk = k + NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1) &
                + NbCellLocal_Ncpus(commRank + 1) + NbWellInjOwn_Ncpus(commRank + 1)
         colk = k + NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1) &
                + NbCellLocal_Ncpus(commRank + 1) + NbWellInjLocal_Ncpus(commRank + 1)

         ! assembly equations for own wells
         if (well_is_own) then
            do n = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
               csrK(JacBigA%Num(n)) = n - JacBigA%Pt(rowk)
            end do
         end if

         ! Pw_Pwmin = IncPressionWellProd(k) - DataWellProdLocal(k)%PressionMin ! P_w - P_w^{min}

         ! Two cases:
         !   if indwell=="p": well equation is Pw-Pwmin = 0
         !   if indwell=="f": well equation is qmol-qw = 0, where qw = sum_{s} sum_{i} q_{w,s,i}

         ! Step 1. well equation if indwell=="p"
         ! Step 2. loop of nodes of well k
         !         - well equation if indwell=="f"
         !         - reservoir equation

         ! Step 1. well equation: Pw - Pwmin = 0
         if (DataWellProdLocal(k)%IndWell == 'p') then
            something_is_produced = .true.
            if (well_is_own) then
               nz = JacBigA%Pt(rowk) + csrK(colk)
               JacBigA%Val(1, 1, nz) = 1.d0
            end if
         end if

         ! Step 2.

         ! nodes of well k
         do s = NodebyWellProdLocal%Pt(k) + 1, NodebyWellProdLocal%Pt(k + 1)
            nums = NodebyWellProdLocal%Num(s) ! nums is node num
            Pws = PerfoWellProd(s)%Pression   ! P_{w,s}
            WIDws = NodeDatabyWellProdLocal%Val(s)%WID ! WID_{w,s}
#ifdef _THERMIQUE_
            WIFws = NodeDatabyWellProdLocal%Val(s)%WIF ! WIF_{w,s}
#endif
            dP_w = 0.d0
            dP_s = 0.d0
#ifdef _THERMIQUE_
            dP_ER_w = 0.d0
            der_ER_s = 0.d0
#endif
            well_node_is_active = .false.

            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)
               Ps = IncNode(nums)%phase_pressure(mph) ! IncNode(nums)%Pression       ! P_s
               Ps_Pws = Ps - Pws
               if ((Ps_Pws > 0.d0) .AND. (.NOT. well_is_closed)) then ! if Ps_Pws < 0 then this term is zero
                  something_is_produced = .true.
                  well_node_is_active = .true.
                  ! derivative of
                  !   sum_{Q_s \cap P_i} q_{w,s,i}
                  !   sum_{Q_s \cap P_i} q_{w,s,e}
                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i
                        ! p^w_s = p^w + {\Delta p}^(n-1)_s
                        ! q^{i \mapsto w}_{i,s} = M^{\alpha}_i WI_s (p_s - p^w_s)
                        dP_w(icp) = dP_w(icp) - DensiteMolaireKrViscoCompNode(icp, mph, nums)*WIDws
                        ! No capillary pressure in the well
                        dP_s(:, icp) = dP_s(:, icp) + divDensiteMolaireKrViscoCompNode(:, icp, mph, nums)*WIDws*Ps_Pws
                        dP_s(:, icp) = dP_s(:, icp) + DensiteMolaireKrViscoCompNode(icp, mph, nums)*WIDws &
                                       *(divPressionNode(:, nums) + divPhasePressureNode(:, mph, nums))
                     end if
                  end do
#ifdef _THERMIQUE_
                  dP_ER_w = dP_ER_w - DensiteMolaireKrViscoEnthalpieNode(mph, nums)*WIDws
                  der_ER_s = der_ER_s + divDensiteMolaireKrViscoEnthalpieNode(:, mph, nums)*WIDws*Ps_Pws
                  der_ER_s = der_ER_s + DensiteMolaireKrViscoEnthalpieNode(mph, nums)*WIDws &
                             *(divPressionNode(:, nums) + divPhasePressureNode(:, mph, nums))
#endif
               end if
            end do ! end loop over phases (Q_s)

            if (well_node_is_active) then

               if (well_is_own) then ! own production well
                  if (DataWellProdLocal(k)%IndWell == 'f') then
                     ! A_kk, k is own production well
                     nz = JacBigA%Pt(rowk) + csrK(colk)
                     do icp = 1, NbComp
                        JacBigA%Val(1, 1, nz) = JacBigA%Val(1, 1, nz) - dP_w(icp) ! term -\sum_{i} q_{w,s,i}, derivative of P_w
                     end do

                     ! A_ks, k is own production well, s is node
                     nz = JacBigA%Pt(rowk) + csrK(nums)
                     do icp = 1, NbComp
                        JacBigA%Val(:, 1, nz) = JacBigA%Val(:, 1, nz) - dP_s(:, icp) ! term -\sum_{i} q_{w,s,i}, derivative of node s
                     end do
                  end if
               end if

               if (IdNodeLocal(nums)%Proc == "o") then ! node own, this node can not be in the boundary (why ? boundary -> dirichlet ?)

                  rows = nums
                  cols = nums
                  do n = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                     csrSR(JacBigA%Num(n)) = n - JacBigA%Pt(rows)
                  end do

                  ! Ask, s is node, k is production well
                  nz = JacBigA%Pt(rows) + csrSR(colk)
                  do icp = 1, NbComp
                     JacBigA%Val(1, icp, nz) = JacBigA%Val(1, icp, nz) + dP_w(icp) ! term q_{w,s,i}, derivative of P_w
                  end do
#ifdef _THERMIQUE_
                  JacBigA%Val(1, NbComp + 1, nz) = JacBigA%Val(1, NbComp + 1, nz) + dP_ER_w
#endif

                  ! Ass, s is node
                  nz = JacBigA%Pt(rows) + csrSR(cols)
                  do icp = 1, NbComp
                     JacBigA%Val(:, icp, nz) = JacBigA%Val(:, icp, nz) + dP_s(:, icp) ! term q_{w,s,i}, derivative of node s
                  end do
#ifdef _THERMIQUE_
                  JacBigA%Val(:, NbComp + 1, nz) = JacBigA%Val(:, NbComp + 1, nz) + der_ER_s
#endif
               end if

            end if

         end do ! end of nodes of well k

         !Make sure Jacobian is not singular
         if (well_is_closed .OR. ((DataWellProdLocal(k)%IndWell == 'f') .AND. (.NOT. something_is_produced))) then
            if (well_is_own) then
               nz = JacBigA%Pt(rowk) + csrK(colk)
               JacBigA%Val(1, 1, nz) = 1.d0
            end if
#ifdef COMPASS_LOG_WELL_INFO
            write (*, *) '[Well Monitoring]:  Jacobian regularization has been performed, since nothing is produced for well', &
               k, 'on proc', commRank + 1
#endif
         end if

      end do !k-well loop

   end subroutine Jacobian_JacBigA_BigSm_wellprod

   ! This subroutine only adds to JacBigA and Sm the contribution from
   ! mswells without coupling with the resevoir

   subroutine Jacobian_JacBigA_BigSm_mswells
      integer :: k, s, unk_idx_s, glb_unk_idx_s, num_s, Nz, Nz_local
      integer :: i, j, l, JacRow_end_idx_s, JacRow_begin_idx_s
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, NbMSWellLocal

         !For all nodes at the well, from queue to head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) !numbering in mesh

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx @ MSWells_JacA
            glb_unk_idx_s = unk_idx_s + mswells_JacBigA_row_off  !unkown idx row @ JacBigA

            !Set default state for mswells if it is closed. TODO: this should be done elsewhere
            if (DataMSWellLocal(k)%IndWell == 'c') IncMSWell(s)%coats%ic = 1 !Default context

            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            if (DataMSWellLocal(k)%IndWell == 'c') then !the mswell is closed
               !TODO: This is already done in JacobianMSWells_JacA_Sm_prod, but we are still separating
               !      the calling of that function and the present function at the python level

               !Nz: Diagonal idx in MSWells_JacA & JacBigA (Recall we built the diagonals terms first for each row in a66)
               !TODO: we can search for Nz manually to make it more robust
               Nz = JacBigA%Pt(glb_unk_idx_s) + 1 + 1 !need to sum the offset col of a61

               if (NbComp .eq. 1) then
                  !Put the Identity matrix
                  JacBigA%Val(1, 1, Nz) = 1.d0
#ifdef _THERMIQUE_
                  JacBigA%Val(NbComp + 1, NbComp + 1, Nz) = 1.d0

#endif
               else if (NbComp .eq. 2) then
                  !This matrix has to be the inverse as the one used for the manual alignment which depends on the default context being used above
                  JacBigA%Val(1, 1, Nz) = 1.d0
                  JacBigA%Val(1, 3, Nz) = -1.d0
                  JacBigA%Val(2, 3, Nz) = 1.d0
                  JacBigA%Val(3, 2, Nz) = 1.d0

               else
#ifndef NDEBUG
                  call CommonMPI_abort( &
                     "Jacobian: mswells closing feature only supported for components no bigger than 2")
#endif

               end if

               cycle !continue s-loop
            else !the mswell is not closed. Then just copy entries from MSWells only  RHS & block-matrix

               bigSm(:, glb_unk_idx_s) = MSWells_Sm(:, unk_idx_s)

               JacRow_begin_idx_s = MSWells_JacA%Pt(unk_idx_s)     !First column at the row of node s
               JacRow_end_idx_s = MSWells_JacA%Pt(unk_idx_s + 1)   !Last column at the row of node s
               !Sum MSWells contribution
               do l = 1, JacRow_end_idx_s - JacRow_begin_idx_s

                  Nz_local = MSWells_JacA%Pt(unk_idx_s) + l
                  Nz = JacBigA%Pt(glb_unk_idx_s) + l + 1 !need to sum the offset col of a61

                  do j = 1, NbIncTotalPrimMax
                     do i = 1, NbComp

                        JacBigA%Val(j, i, Nz) = JacBigA%Val(j, i, Nz) + MSWells_JacA%Val(j, i, Nz_local)
                     end do
#ifdef _THERMIQUE_

                     JacBigA%Val(j, NbComp + 1, Nz) = JacBigA%Val(j, NbComp + 1, Nz) + MSWells_JacA%Val(j, NbComp + 1, Nz_local)
#endif
                  end do
               end do
            end if

         end do ! end of loop s

      end do! end of loop NbMSWellLocal

   end subroutine Jacobian_JacBigA_BigSm_mswells

   ! mswells with coupling with the resevoir
   ! No coupling at the head of mswell
   subroutine Jacobian_JacBigA_BigSm_coupling_mswells
      integer :: k, s, unk_idx_s, glb_unk_idx_s, glb_unk_col_idx_s, nums, nz
      integer :: rows, cols, m, mph, n, icp
      integer :: i, j, l
      double precision :: Pws, Ps, WIDws, WIFws, Ps_Pws
      real(c_double) :: dP_w(NbComp), dP_s(NbIncTotalPrimMax, NbComp) ! NbComp mass balance equation
      real(c_double) :: dP_ER_w, der_ER_s(NbIncTotalPrimMax) ! energy balance equation
      logical :: something_is_produced

#ifndef ComPASS_DIPHASIC_CONTEXT
      if (NbMSWellLocal > 0) then
         call CommonMPI_abort( &
            "In function Jacobian_JacBigA_BigSm_coupling_mswells: Multi-segmented wells are only implemented for diphasic physics!")
      endif
#else

      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         !For all nodes at the well, ** from queue to the node before the head **
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            nums = NodebyMSWellLocal%Num(s) !numbering in mesh

            if (IdNodeLocal(nums)%Proc == "g") cycle ! node not own

            something_is_produced = .false.
            Pws = IncMSWell(s)%coats%Pression ! P_{w,s}
            WIDws = NodeDatabyMSWellLocal%Val(s)%WID ! WID_{w,s}
#ifdef _THERMIQUE_
            WIFws = NodeDatabyMSWellLocal%Val(s)%WIF ! WIF_{w,s}
#endif
            dP_w = 0.d0
            dP_s = 0.d0
#ifdef _THERMIQUE_
            dP_ER_w = 0.d0
            der_ER_s = 0.d0
#endif

            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)
               Ps = IncNode(nums)%phase_pressure(mph) ! IncNode(nums)%Pression       ! P_s
               Ps_Pws = Ps - Pws
               if ((Ps_Pws > 0.d0)) then ! if Ps_Pws < 0 then this term is zero
                  something_is_produced = .true.
                  ! derivative of
                  !   sum_{Q_s \cap P_i} q_{w,s,i}
                  !   sum_{Q_s \cap P_i} q_{w,s,e}
                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i
                        ! p^w_s = p^w + {\Delta p}^(n-1)_s
                        ! q^{i \mapsto w}_{i,s} = M^{\alpha}_i WI_s (p_s - p^w_s)
                        dP_w(icp) = dP_w(icp) - DensiteMolaireKrViscoCompNode(icp, mph, nums)*WIDws
                        ! No capillary pressure in the well
                        dP_s(:, icp) = dP_s(:, icp) + divDensiteMolaireKrViscoCompNode(:, icp, mph, nums)*WIDws*Ps_Pws
                        dP_s(:, icp) = dP_s(:, icp) + DensiteMolaireKrViscoCompNode(icp, mph, nums)*WIDws &
                                       *(divPressionNode(:, nums) + divPhasePressureNode(:, mph, nums))
                     end if
                  end do
#ifdef _THERMIQUE_
                  dP_ER_w = dP_ER_w - DensiteMolaireKrViscoEnthalpieNode(mph, nums)*WIDws
                  der_ER_s = der_ER_s + divDensiteMolaireKrViscoEnthalpieNode(:, mph, nums)*WIDws*Ps_Pws
                  der_ER_s = der_ER_s + DensiteMolaireKrViscoEnthalpieNode(mph, nums)*WIDws &
                             *(divPressionNode(:, nums) + divPhasePressureNode(:, mph, nums))
#endif
               end if
            end do ! end loop over phases (Q_s)

            if (something_is_produced) then

               unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx @ MSWells_JacA
               glb_unk_idx_s = unk_idx_s + mswells_JacBigA_row_off  !unkown row-idx @ JacBigA
               glb_unk_col_idx_s = unk_idx_s + mswells_JacBigA_unk_lcl_off  !unkown col-idx @ JacBigA

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !MSWell node eq.

               !The first col-idx is a61
               nz = JacBigA%Pt(glb_unk_idx_s) + 1
               do icp = 1, NbComp
                  JacBigA%Val(:, icp, nz) = JacBigA%Val(:, icp, nz) - dP_s(:, icp) ! term q_{w,s,i}, derivative of node s
               end do
#ifdef _THERMIQUE_
               JacBigA%Val(:, NbComp + 1, nz) = JacBigA%Val(:, NbComp + 1, nz) - der_ER_s
#endif

               !The second col-idx is a66, and the first idx is the diagonal
               nz = JacBigA%Pt(glb_unk_idx_s) + 2
               do icp = 1, NbComp
                  JacBigA%Val(1, icp, nz) = JacBigA%Val(1, icp, nz) - dP_w(icp) ! term q_{w,s,i}, derivative of P_w
               end do
#ifdef _THERMIQUE_
               JacBigA%Val(1, NbComp + 1, nz) = JacBigA%Val(1, NbComp + 1, nz) - dP_ER_w
#endif

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Reservoir eq.

               rows = nums
               cols = nums

               do n = JacBigA%Pt(rows) + 1, JacBigA%Pt(rows + 1)
                  csrSR(JacBigA%Num(n)) = n - JacBigA%Pt(rows)
               end do

               ! Ask, s is node, k is production mswell
               nz = JacBigA%Pt(rows) + csrSR(glb_unk_col_idx_s)
               do icp = 1, NbComp
                  JacBigA%Val(1, icp, nz) = JacBigA%Val(1, icp, nz) + dP_w(icp) ! term q_{w,s,i}, derivative of P_w
               end do
#ifdef _THERMIQUE_
               JacBigA%Val(1, NbComp + 1, nz) = JacBigA%Val(1, NbComp + 1, nz) + dP_ER_w
#endif

               ! Ass, s is node
               nz = JacBigA%Pt(rows) + csrSR(cols)
               do icp = 1, NbComp
                  JacBigA%Val(:, icp, nz) = JacBigA%Val(:, icp, nz) + dP_s(:, icp) ! term q_{w,s,i}, derivative of node s
               end do
#ifdef _THERMIQUE_
               JacBigA%Val(:, NbComp + 1, nz) = JacBigA%Val(:, NbComp + 1, nz) + der_ER_s
#endif

            end if !something_is_produced

         end do !s
      end do !k
#endif

   end subroutine Jacobian_JacBigA_BigSm_coupling_mswells

#ifdef _WITH_FREEFLOW_STRUCTURES_
   ! loop of node, index is nums
   ! 1.1 div
   subroutine Jacobian_JacBigA_BigSm_FF_node
      ! div prims and Sm from FreeFlow term
      double precision :: &
         divS3(NbIncTotalPrimMax, NbComp), & ! s for node, represent s in paper
         Sm3(NbComp), &
         divTFF(NbIncTotalPrimMax), &
         SmTFF
      integer :: nums, i, j, nz

      do nums = 1, NbNodeOwn_Ncpus(commRank + 1)

         if (IsFreeflowNode(nums)) then ! loop over freeflow dof only

            ! compute the contribution of the freeflow
            call Jacobian_divMolarFreeFlow_node( &
               NbIncTotalPrim_ctx, &  ! does not depend on k or s
               AtmState(nums)%Comp, IncNode(nums), SurfFreeFlowLocal(nums), &
               divFreeFlowMolarFlowrateCompNode(:, :, :, nums), SmFreeFlowMolarFlowrateCompNode(:, :, nums), &
               divFreeFlowMolarFlowrateNode(:, :, nums), SmFreeFlowMolarFlowrateNode(:, nums), &
               divFreeFlowHmCompNode(:, :, :, nums), SmFreeFlowHmCompNode(:, :, nums), &
               divS3, Sm3)

#ifdef _THERMIQUE_
            ! compute the thermal contribution of the freeflow
            call Jacobian_divThermalFreeFlow_node( &
               NbIncTotalPrim_ctx, &  ! does not depend on k or s
               IncNode(nums), SurfFreeFlowLocal(nums), AtmEnthalpieNode(:, nums), &
               divFreeFlowMolarFlowrateEnthalpieNode(:, :, nums), &
               SmFreeFlowMolarFlowrateEnthalpieNode(:, nums), &
               divFreeFlowMolarFlowrateNode(:, :, nums), &
               SmFreeFlowMolarFlowrateNode(:, nums), &
               divFreeFlowHTTemperatureNetRadiationNode(:, nums), &
               SmFreeFlowHTTemperatureNetRadiationNode(nums), &
               divTFF, SmTFF)
#endif

            ! the diagonal element (s,s) is JacBigA(nz)
            do i = JacBigA%Pt(nums) + 1, JacBigA%Pt(nums + 1)
               if (JacBigA%Num(i) == nums) then
                  nz = i
                  exit
               end if
            end do

            if (IdNodeLocal(nums)%P /= "d") then
               ! JacBigA%Val(:,:,nz)
               do i = 1, NbComp
                  do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                     JacBigA%Val(j, i, nz) = JacBigA%Val(j, i, nz) + divS3(j, i)
                  end do
               end do
               ! Sm
               bigSm(1:NbComp, nums) = bigSm(1:NbComp, nums) - Sm3(1:NbComp)

            end if ! Dirichlet node

#ifdef _THERMIQUE_
            if (IdNodeLocal(nums)%T /= "d") then
               ! ps. divFourierFlux_s=0
               do j = 1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  JacBigA%Val(j, NbComp + 1, nz) = JacBigA%Val(j, NbComp + 1, nz) + divTFF(j)
               end do
               ! Sm
               bigSm(NbComp + 1, nums) = bigSm(NbComp + 1, nums) - SmTFF
            end if ! Dirichlet node
#endif

         endif ! FreeFlow node

      enddo ! node nums

   end subroutine Jacobian_JacBigA_BigSm_FF_node

   ! Derivatives of the FreeFlow terms in the molar balance equations
   pure subroutine Jacobian_divMolarFreeFlow_node( &
      NbIncTotalPrim, &
      atm_comp, Incs, area_s, &
      divFreeFlowMolarFlowrateComp_s, SmFreeFlowMolarFlowrateComp_s, &
      divFreeFlowMolarFlowrate_s, SmFreeFlowMolarFlowrate_s, &
      divFreeFlowHmComp_s, SmFreeFlowHmComp_s, &
      divS, Sm0)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), intent(in) :: Incs
      double precision, intent(in) :: &
         atm_comp(NbComp, NbPhase), area_s, &
         divFreeFlowMolarFlowrateComp_s(NbIncTotalPrimMax, NbComp, NbPhase), &
         divFreeFlowMolarFlowrate_s(NbIncTotalPrimMax, NbPhase), &
         divFreeFlowHmComp_s(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmFreeFlowMolarFlowrateComp_s(NbComp, NbPhase), &
         SmFreeFlowMolarFlowrate_s(NbPhase), &
         SmFreeFlowHmComp_s(NbComp, NbPhase)

      double precision, intent(out) :: &
         divS(NbIncTotalPrimMax, NbComp), &
         Sm0(NbComp)

      integer :: m, mph, icp, j

      divS(:, :) = 0.d0
      Sm0(:) = 0.d0

      ! -> divS, node
      do m = 1, NbPhasePresente_ctx(Incs%ic) ! Q_s
         mph = NumPhasePresente_ctx(m, Incs%ic)

         if (Incs%FreeFlow_flowrate(mph) >= 0.d0) then

            ! To understand better, change the order of the loop do m=.. and the loop do icp=..
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i

                  do j = 1, NbIncTotalPrim(Incs%ic)
                     divS(j, icp) = divS(j, icp) + area_s*( &
                                    divFreeFlowMolarFlowrateComp_s(j, icp, m) + &
                                    divFreeFlowHmComp_s(j, icp, m))
                  end do

                  ! Sm0
                  Sm0(icp) = Sm0(icp) + area_s*( &
                             SmFreeFlowMolarFlowrateComp_s(icp, m) + &
                             SmFreeFlowHmComp_s(icp, m))
               end if
            end do ! end of icp

         else ! Incs%FreeFlow_flowrate(mph)<0.d0
            ! liq phase never enters in this loop because always FreeFlow_flowrate(liq)>=0.d0

            ! To understand better, change the order of the loop do m=.. and the loop do icp=..
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i

                  do j = 1, NbIncTotalPrim(Incs%ic)
                     divS(j, icp) = divS(j, icp) + area_s*( & ! only gas phase because FreeFlow_flowrate(liq)>=0.d0
                                    divFreeFlowMolarFlowrate_s(j, m)*atm_comp(icp, mph) + & ! atm_comp is a constant, no derivative wrt atm_comp
                                    divFreeFlowHmComp_s(j, icp, m))
                  end do

                  Sm0(icp) = Sm0(icp) + area_s*( &
                             SmFreeFlowMolarFlowrate_s(m)*atm_comp(icp, mph) + &
                             SmFreeFlowHmComp_s(icp, m))
               end if
            end do ! end of icp
         endif

      enddo

   end subroutine Jacobian_divMolarFreeFlow_node

   ! Derivatives of the FreeFlow terms in the energy balance equation
   pure subroutine Jacobian_divThermalFreeFlow_node( &
      NbIncTotalPrim, &
      Incs, area_s, AtmEnthalpie_s, &
      divFreeFlowMolarFlowrateEnthalpie_s, SmFreeFlowMolarFlowrateEnthalpie_s, &
      divFreeFlowMolarFlowrate_s, SmFreeFlowMolarFlowrate_s, &
      divFreeFlowHTTemperatureNetRadiation_s, SmFreeFlowHTTemperatureNetRadiation_s, &
      divS, Sm0)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), intent(in) :: Incs
      double precision, intent(in) :: &
         area_s, &
         divFreeFlowMolarFlowrateEnthalpie_s(NbIncTotalPrimMax, NbPhase), &
         divFreeFlowMolarFlowrate_s(NbIncTotalPrimMax, NbPhase), &
         divFreeFlowHTTemperatureNetRadiation_s(NbIncTotalPrimMax), &
         SmFreeFlowMolarFlowrateEnthalpie_s(NbPhase), &
         SmFreeFlowMolarFlowrate_s(NbPhase), &
         SmFreeFlowHTTemperatureNetRadiation_s, &
         AtmEnthalpie_s(NbPhase)

      double precision, intent(out) :: &
         divS(NbIncTotalPrimMax), Sm0

      integer :: m, mph, j

      divS = 0.d0
      Sm0 = 0.d0

      ! -> divS, node
      do m = 1, NbPhasePresente_ctx(Incs%ic) ! Q_s
         mph = NumPhasePresente_ctx(m, Incs%ic)

         if (Incs%FreeFlow_flowrate(mph) >= 0.d0) then

            do j = 1, NbIncTotalPrim(Incs%ic)
               divS(j) = divS(j) + area_s* &
                         divFreeFlowMolarFlowrateEnthalpie_s(j, m)
            enddo
            Sm0 = Sm0 + area_s*SmFreeFlowMolarFlowrateEnthalpie_s(m)

         else ! Incs%FreeFlow_flowrate(mph)<0.d0
            ! liq phase never enters in this loop because always FreeFlow_flowrate(liq)>=0.d0
            do j = 1, NbIncTotalPrim(Incs%ic)
               divS(j) = divS(j) + area_s* &
                         divFreeFlowMolarFlowrate_s(j, m)*AtmEnthalpie_s(mph) ! AtmEnthalpie_s is a constant, no derivative
            enddo
            Sm0 = Sm0 + area_s* &
                  SmFreeFlowMolarFlowrate_s(m)*AtmEnthalpie_s(mph)

         endif ! sign of flux
      enddo

      do j = 1, NbIncTotalPrim(Incs%ic)
         divS(j) = divS(j) + area_s* &
                   divFreeFlowHTTemperatureNetRadiation_s(j)
      enddo
      Sm0 = Sm0 + area_s*SmFreeFlowHTTemperatureNetRadiation_s

   end subroutine Jacobian_divThermalFreeFlow_node
#endif

   !> term: \sum{P_i \cap Q_{k or s} } &
   !!          F_i = Comp * DensiteMolaire * PermRel / Viscosite ) * FluxDarcyKI
   !!          (div (Comp * DensiteMolaire * PermRel / Viscosite ) * FluxDarcyKI)
   !!       = divK * div(X_k) + divS * div(X_s) + Sm
   !!
   !!       where k is cell, s is node
   !! compute:
   !!       divK, divS, Sm
   !!
   !! if FluxDarcyKI>=0 then
   !!    \sum{P_i \cap Q_{k} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
   !! else
   !!    \sum{P_i \cap Q_{s} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
   !!
   !! it is equivalent to:
   !!    \sum{P_i \cap Q_{k} \cap FluxDarcyKI>=0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
   !!    \sum{P_i \cap Q_{s} \cap FluxDarcyKI<0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
   pure subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux( &
      NbIncTotalPrim, & ! does not depend on k or s
      ic_k, ic_s, &
      DarcyFlux, &
      divDensiteMolaireKrViscoComp_k, SmDensiteMolaireKrViscoComp_k, &
      divDensiteMolaireKrViscoComp_s, SmDensiteMolaireKrViscoComp_s, &
      divK, divS, Sm0)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      integer(c_int), intent(in) :: ic_k, ic_s
      double precision, intent(in) :: &
         DarcyFlux(NbPhase), &
         divDensiteMolaireKrViscoComp_k(NbIncTotalPrimMax, NbComp, NbPhase), &
         divDensiteMolaireKrViscoComp_s(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmDensiteMolaireKrViscoComp_k(NbComp, NbPhase), &
         SmDensiteMolaireKrViscoComp_s(NbComp, NbPhase)

      double precision, intent(out) :: &
         divK(NbIncTotalPrimMax, NbComp), &
         divS(NbIncTotalPrimMax, NbComp), &
         Sm0(NbComp)

      integer :: m, mph, icp, j

      divK(:, :) = 0.d0
      divS(:, :) = 0.d0
      Sm0(:) = 0.d0

      ! -> divK, upwind is k, divS=0
      do m = 1, NbPhasePresente_ctx(ic_k) ! Q_k
         mph = NumPhasePresente_ctx(m, ic_k)

         if (DarcyFlux(mph) >= 0.d0) then

            ! To understand better, change the order of the loop do m=.. and the loop do icp=...
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i

                  do j = 1, NbIncTotalPrim(ic_k)
                     divK(j, icp) = divK(j, icp) + &
                                    divDensiteMolaireKrViscoComp_k(j, icp, mph)*DarcyFlux(mph)

                  end do

                  Sm0(icp) = Sm0(icp) + &
                             SmDensiteMolaireKrViscoComp_k(icp, mph)*DarcyFlux(mph)
               end if
            end do ! end of icp

         end if
      end do ! end of Q_k

      ! -> divS, upwind is s, divK=0
      do m = 1, NbPhasePresente_ctx(ic_s) ! Q_s
         mph = NumPhasePresente_ctx(m, ic_s)

         if (DarcyFlux(mph) < 0.d0) then

            ! To understand better, change the order of the loop do m=.. and the loop do icp=..
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i

                  do j = 1, NbIncTotalPrim(ic_s)
                     divS(j, icp) = divS(j, icp) + &
                                    divDensiteMolaireKrViscoComp_s(j, icp, mph)*DarcyFlux(mph)
                  end do

                  Sm0(icp) = Sm0(icp) + SmDensiteMolaireKrViscoComp_s(icp, mph)*DarcyFlux(mph)

               end if
            end do ! end of icp

         end if

      end do ! end of Q_s

   end subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux

   ! term: \sum{P_i \cap (Q_k \cap Q_s)} &
   !           DensiteMolaire * PermRel / Viscosite * Comp * div(FluxDarcyKI)
   !       = ( divK * div(X_k)
   !         + divS * div(X_s)
   !         + \sum_r div(X_r) * div(X_r)
   !         + Sm
   ! compute
   !        divK, divS, divR, Sm
   pure subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux( &
      NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
      ic_k, ic_s, DarcyFlux, Nodebyk, Fracbyk, &
      DensiteMolaireKrViscoComp_k, &
      DensiteMolaireKrViscoComp_s, &
      divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
      divK, divS, divR, Sm0)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_k, ic_s
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, intent(in) :: &
         DarcyFlux(NbPhase), &
         DensiteMolaireKrViscoComp_k(NbComp, NbPhase), &
         DensiteMolaireKrViscoComp_s(NbComp, NbPhase), &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(in) :: divDarcyFlux_r ! r represen

      double precision, intent(out) :: &
         divK(NbIncTotalPrimMax, NbComp), &
         divS(NbIncTotalPrimMax, NbComp), &
         Sm0(NbComp)
      double precision, dimension(:, :, :), intent(out) :: divR

      integer :: m, mph
      integer :: NbNode_in_k, NbFrac_in_k

      divK(:, :) = 0.d0
      divS(:, :) = 0.d0
      divR(:, :, :) = 0.d0
      Sm0(:) = 0.d0

      ! number of nodes/fracs in k
      NbNode_in_k = size(Nodebyk)
      NbFrac_in_k = size(Fracbyk)

      ! sum_{P_i \cap Q_{k} \cap FluxDarcyKI>=0}
      do m = 1, NbPhasePresente_ctx(ic_k) ! Q_k
         mph = NumPhasePresente_ctx(m, ic_k)

         if (DarcyFlux(mph) >= 0.d0) then
            call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_mph( &
               NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
               mph, ic_k, ic_s, &
               NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
               DensiteMolaireKrViscoComp_k, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK, divS, divR, Sm0)
         end if
      end do ! end of Q_k

      ! sum_{P_i \cap Q_{s} \cap DarcyFlux<0}
      do m = 1, NbPhasePresente_ctx(ic_s) ! Q_s
         mph = NumPhasePresente_ctx(m, ic_s)

         if (DarcyFlux(mph) < 0.d0) then
            call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_mph( &
               NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
               mph, ic_k, ic_s, &
               NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
               DensiteMolaireKrViscoComp_s, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK, divS, divR, Sm0)
         end if

      end do ! end of Q_s

   end subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux

   pure subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_mph( &
      NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s, FaceToFracLocal is missing
      mph, &
      ic_k, ic_s, NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
      DensiteMolaireKrViscoComp_up, &
      divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
      divK, divS, divR, Sm0)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_k, ic_s
      integer, intent(in) :: mph, NbNode_in_k, NbFrac_in_k
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, intent(in) :: &
         DensiteMolaireKrViscoComp_up(NbComp, NbPhase), &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(in) :: divDarcyFlux_r ! r represen

      double precision, intent(inout) :: &
         divK(NbIncTotalPrimMax, NbComp), &
         divS(NbIncTotalPrimMax, NbComp), &
         Sm0(NbComp)
      double precision, dimension(:, :, :), intent(inout) :: divR

      integer :: j, icp, r, numr, rf

      do icp = 1, NbComp ! P_i
         if (MCP(icp, mph) == 1) then

            do j = 1, NbIncTotalPrim(ic_k) ! divK
               divK(j, icp) = divK(j, icp) &
                              + divDarcyFlux_k(j, mph)*DensiteMolaireKrViscoComp_up(icp, mph)
            end do

            do j = 1, NbIncTotalPrim(ic_s) ! divS
               divS(j, icp) = divS(j, icp) &
                              + divDarcyFlux_s(j, mph)*DensiteMolaireKrViscoComp_up(icp, mph)
            end do

            ! divR for r is node in dof(k)
            do r = 1, NbNode_in_k
               numr = Nodebyk(r)

               do j = 1, NbIncTotalPrim(IncNode(numr)%ic)
                  divR(j, icp, r) = divR(j, icp, r) &
                                    + divDarcyFlux_r(j, mph, r)*DensiteMolaireKrViscoComp_up(icp, mph)
               end do
            end do

            ! divR for r is frac in dof(k)
            do r = 1, NbFrac_in_k
               numr = Fracbyk(r)
               rf = r + NbNode_in_k

               do j = 1, NbIncTotalPrim(IncFrac(numr)%ic)
                  divR(j, icp, rf) = divR(j, icp, rf) &
                                     + divDarcyFlux_r(j, mph, rf)*DensiteMolaireKrViscoComp_up(icp, mph)
               end do
            end do

            Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph)*DensiteMolaireKrViscoComp_up(icp, mph)
         end if
      end do ! end of icp

   end subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_mph

   pure subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux( &
      NbIncTotalPrim, & ! does not depend on k or s
      IncNode, IncFrac, ic_k, ic_s, DarcyFlux, &
      Nodebyk, Fracbyk, &
      DensiteMolaireKrViscoEnthalpie_k, &
      divDensiteMolaireKrViscoEnthalpie_k, &
      SmDensiteMolaireKrViscoEnthalpie_k, &
      DensiteMolaireKrViscoEnthalpie_s, &
      divDensiteMolaireKrViscoEnthalpie_s, &
      SmDensiteMolaireKrViscoEnthalpie_s, &
      divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
      divEgK, divEgS, divEgR, SmEg)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_k, ic_s
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, intent(in) :: &
         DarcyFlux(NbPhase), &
         DensiteMolaireKrViscoEnthalpie_k(NbPhase), &
         DensiteMolaireKrViscoEnthalpie_s(NbPhase), &
         divDensiteMolaireKrViscoEnthalpie_k(NbIncTotalPrimMax, NbPhase), &
         divDensiteMolaireKrViscoEnthalpie_s(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaireKrViscoEnthalpie_k(NbPhase), &
         SmDensiteMolaireKrViscoEnthalpie_s(NbPhase), &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(in) :: divDarcyFlux_r

      double precision, intent(out) :: &
         divEgK(NbIncTotalPrimMax), &
         divEgS(NbIncTotalPrimMax), &
         SmEg
      double precision, dimension(:, :), intent(out) :: divEgR

      integer :: m, mph
      integer :: NbNode_in_k, NbFrac_in_k

      divEgK(:) = 0.d0
      divEgS(:) = 0.d0
      divEgR(:, :) = 0.d0
      SmEg = 0.d0

      ! number of nodes/fracs in k
      NbNode_in_k = size(Nodebyk)
      NbFrac_in_k = size(Fracbyk)

      ! upwind is k, DensiteMolaireKrViscoEnthalpie(k)*DarcyFlux
      do m = 1, NbPhasePresente_ctx(ic_k) ! Q_k
         mph = NumPhasePresente_ctx(m, ic_k)

         if (DarcyFlux(mph) >= 0.d0) then
            call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_mph( &
               NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
               mph, ic_k, ic_s, DarcyFlux, &
               NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
               DensiteMolaireKrViscoEnthalpie_k, &
               divDensiteMolaireKrViscoEnthalpie_k, &
               SmDensiteMolaireKrViscoEnthalpie_k, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

         end if
      end do

      ! upwind is s
      do m = 1, NbPhasePresente_ctx(ic_s) ! Q_s
         mph = NumPhasePresente_ctx(m, ic_s)

         if (DarcyFlux(mph) < 0.d0) then
            call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_mph( &
               NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
               mph, ic_s, ic_k, DarcyFlux, &
               NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
               DensiteMolaireKrViscoEnthalpie_s, &
               divDensiteMolaireKrViscoEnthalpie_s, &
               SmDensiteMolaireKrViscoEnthalpie_s, &
               divDarcyFlux_s, divDarcyFlux_k, divDarcyFlux_r, SmDarcyFlux, &
               divEgS, divEgK, divEgR, SmEg)
         end if
      end do

   end subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux

   pure subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_mph( &
      NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
      mph, ic_up, ic_down, DarcyFlux, & ! miss FaceToFracLocal which will disapear
      NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
      DensiteMolaireKrViscoEnthalpie_up, &
      divDensiteMolaireKrViscoEnthalpie_up, &
      SmDensiteMolaireKrViscoEnthalpie_up, &
      divDarcyFlux_up, divDarcyFlux_down, divDarcyFlux_r, SmDarcyFlux, &
      divEg_up, divEg_down, divEgR, SmEg)

      integer, intent(in) :: &
         NbIncTotalPrim(NbContexte), &
         mph, NbNode_in_k, NbFrac_in_k
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_up, ic_down
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, intent(in) :: &
         DarcyFlux(NbPhase), &
         DensiteMolaireKrViscoEnthalpie_up(NbPhase), &
         divDensiteMolaireKrViscoEnthalpie_up(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaireKrViscoEnthalpie_up(NbPhase), &
         divDarcyFlux_up(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_down(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(in) :: divDarcyFlux_r

      double precision, intent(inout) :: &
         divEg_up(NbIncTotalPrimMax), &
         divEg_down(NbIncTotalPrimMax), &
         SmEg
      double precision, dimension(:, :), intent(inout) :: divEgR

      integer :: j, r, numr, rf

      ! divEg_up
      do j = 1, NbIncTotalPrim(ic_up)
         divEg_up(j) = divEg_up(j) &
                       + divDensiteMolaireKrViscoEnthalpie_up(j, mph)*DarcyFlux(mph) &
                       + divDarcyFlux_up(j, mph)*DensiteMolaireKrViscoEnthalpie_up(mph)
      end do

      ! divEg_down
      do j = 1, NbIncTotalPrim(ic_down)
         divEg_down(j) = divEg_down(j) &
                         + divDarcyFlux_down(j, mph)*DensiteMolaireKrViscoEnthalpie_up(mph)
      end do

      ! divEgR
      do r = 1, NbNode_in_k ! divR for r is node in dof(k)
         numr = Nodebyk(r)

         do j = 1, NbIncTotalPrim(IncNode(numr)%ic)
            divEgR(j, r) = divEgR(j, r) &
                           + divDarcyFlux_r(j, mph, r)*DensiteMolaireKrViscoEnthalpie_up(mph)
         end do
      end do

      ! divEgR, r is frac in dof(k)
      do r = 1, NbFrac_in_k
         numr = Fracbyk(r) ! numr is frac num
         rf = r + NbNode_in_k

         do j = 1, NbIncTotalPrim(IncFrac(numr)%ic)
            divEgR(j, rf) = divEgR(j, rf) &
                            + divDarcyFlux_r(j, mph, rf)*DensiteMolaireKrViscoEnthalpie_up(mph)
         end do
      end do

      SmEg = SmEg &
             + SmDensiteMolaireKrViscoEnthalpie_up(mph)*DarcyFlux(mph) &
             + SmDarcyFlux(mph)*DensiteMolaireKrViscoEnthalpie_up(mph)

   end subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_mph

   ! div: div prim
   ! div (V_{k,s}^alpha ) = divDarcyFlux_k * div(X_k)
   !                      + divDarcyFlux_s * div(X_s)
   !                      + \sum_{r \in V_k} divDarcyFlux_r * div(X_r)
   !                      + SmDarcyFlux
   pure subroutine Jacobian_divDarcyFlux( &
      NbIncTotalPrim, phase_can_be_present, IncNode, IncFrac, & ! does not depend on k or s
      ic_k, ic_s, zk, zNode, zFrac, gravity, & ! XNodeLocal(3, numr)  XFaceLocal(3, numr)
      Nodebyk, Fracbyk, &
      TkLocal_Darcy_ks, &
      divPression_k, divPhasePressure_k, divrho_k, &
      SmPression_k, Smrho_k, &
      divrho_s, Smrho_s, &
      divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
      divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
      divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
      SmDarcyFlux)

      integer, intent(in) :: NbIncTotalPrim(NbContexte)
      logical(c_bool), intent(in) :: phase_can_be_present(NbPhase, NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_k, ic_s
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, dimension(:), intent(in) :: &
         zNode, zFrac, TkLocal_Darcy_ks, &
         SmPressionNode, SmPressionFrac
      double precision, dimension(:, :), intent(in) :: &
         divPressionNode, divPressionFrac
      double precision, dimension(:, :, :), intent(in) :: &
         divPhasePressureNode, divPhasePressureFrac
      double precision, intent(in) :: &
         zk, gravity, &
         divPression_k(NbIncTotalPrimMax), &
         divPhasePressure_k(NbIncTotalPrimMax, NbPhase), &
         divrho_k(NbIncTotalPrimMax, NbPhase), &
         divrho_s(NbIncTotalPrimMax, NbPhase), &
         SmPression_k, &
         Smrho_k(NbPhase), &
         Smrho_s(NbPhase)

      double precision, intent(out) :: &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(out) :: divDarcyFlux_r ! r represent s' in paper

      integer :: NbNode_in_k, NbFrac_in_k

      integer :: r, numr, rf, m, mph
      double precision :: sum_aks, sum_aksgz

      divDarcyFlux_k(:, :) = 0.d0
      divDarcyFlux_s(:, :) = 0.d0
      divDarcyFlux_r(:, :, :) = 0.d0
      SmDarcyFlux(:) = 0.d0

      ! number of nodes/fracs in k
      NbNode_in_k = size(Nodebyk)
      NbFrac_in_k = size(Fracbyk)

      ! sum_aks = \sum_r a_{k,s}^r
      ! sum_aksgz = \sum_r a_{k,s}^r * g * (z_k-z_r)
      sum_aks = 0.d0
      sum_aksgz = 0.d0

      do r = 1, NbNode_in_k ! r is node
         numr = Nodebyk(r) ! numr is node num here

         sum_aks = sum_aks + TkLocal_Darcy_ks(r)
         sum_aksgz = sum_aksgz + TkLocal_Darcy_ks(r)*(zk - zNode(numr))
      end do

      do r = 1, NbFrac_in_k ! r is frac
         numr = FracToFaceLocal(Fracbyk(r)) ! numr is face num here
         rf = r + NbNode_in_k

         sum_aks = sum_aks + TkLocal_Darcy_ks(rf)
         sum_aksgz = sum_aksgz + TkLocal_Darcy_ks(rf)*(zk - zFrac(numr))
      end do

      sum_aksgz = sum_aksgz*gravity

      do m = 1, NbPhasePresente_ctx(ic_k) ! Q_k
         mph = NumPhasePresente_ctx(m, ic_k)

         call Jacobian_divDarcyFlux_mph( &
            NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
            mph, ic_k, ic_s, &
            NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
            TkLocal_Darcy_ks, &
            divPression_k, divPhasePressure_k, divrho_k, &
            SmPression_k, Smrho_k, &
            divrho_s, Smrho_s, &
            divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
            divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
            sum_aks, sum_aksgz, &
            divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
            SmDarcyFlux)
      end do

      do m = 1, NbPhasePresente_ctx(ic_s) ! Q_s
         mph = NumPhasePresente_ctx(m, ic_s)

         if (.not. phase_can_be_present(mph, ic_k)) then ! this phase is not in Q_k

            call Jacobian_divDarcyFlux_mph( &
               NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s
               mph, ic_k, ic_s, &
               NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
               TkLocal_Darcy_ks, &
               divPression_k, divPhasePressure_k, divrho_k, &
               SmPression_k, Smrho_k, &
               divrho_s, Smrho_s, &
               divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
               divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
               sum_aks, sum_aksgz, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
               SmDarcyFlux)
         end if ! end of Id_Qki(mph)
      end do ! end of Q_s

   end subroutine Jacobian_divDarcyFlux

   pure subroutine Jacobian_divDarcyFlux_mph( &
      NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s, FaceToFracLocal is missing
      mph, ic_k, ic_s, &
      NbNode_in_k, Nodebyk, NbFrac_in_k, Fracbyk, &
      TkLocal_Darcy_ks, &
      divPression_k, divPhasePressure_k, divrho_k, &
      SmPression_k, Smrho_k, &
      divrho_s, Smrho_s, &
      divPressionNode, divPhasePressureNode, SmPressionNode, & ! does not depend on k or s
      divPressionFrac, divPhasePressureFrac, SmPressionFrac, & ! does not depend on k or s
      sum_aks, sum_aksgz, &
      divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
      SmDarcyFlux)

      integer, intent(in) :: &
         NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer, intent(in) :: mph, NbNode_in_k, NbFrac_in_k
      integer(c_int), intent(in) :: ic_k, ic_s
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, dimension(:), intent(in) :: &
         TkLocal_Darcy_ks, &
         SmPressionNode, SmPressionFrac
      double precision, dimension(:, :), intent(in) :: &
         divPressionNode, divPressionFrac
      double precision, dimension(:, :, :), intent(in) :: &
         divPhasePressureNode, divPhasePressureFrac
      double precision, intent(in) :: &
         divPression_k(NbIncTotalPrimMax), &
         divPhasePressure_k(NbIncTotalPrimMax, NbPhase), &
         divrho_k(NbIncTotalPrimMax, NbPhase), &
         divrho_s(NbIncTotalPrimMax, NbPhase), &
         SmPression_k, &
         Smrho_k(NbPhase), &
         Smrho_s(NbPhase), &
         sum_aks, sum_aksgz

      double precision, intent(inout) :: &
         divDarcyFlux_k(NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s(NbIncTotalPrimMax, NbPhase), &
         SmDarcyFlux(NbPhase)
      double precision, dimension(:, :, :), intent(inout) :: divDarcyFlux_r ! r represent s' in paper

      integer :: r, numr, rf, j

      ! divDarcyFlux_k
      do j = 1, NbIncTotalPrim(ic_k)
         divDarcyFlux_k(j, mph) = &
            sum_aks*divPression_k(j) &        ! \sum a_{ks}^{s'} P_k
            + sum_aks*divPhasePressure_k(j, mph) & ! \sum a_{ks}^{s'} PressionCap_k
            + sum_aksgz*divrho_k(j, mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
      end do

      ! SmDarcyFlux from k
      SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                         + sum_aks*SmPression_k &
                         + sum_aksgz*Smrho_k(mph)    ! SmPressionCap=0

      ! divDarcyFlux_s
      do j = 1, NbIncTotalPrim(ic_s)
         divDarcyFlux_s(j, mph) = &
            sum_aksgz*divrho_s(j, mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
      end do

      ! SmDarcyFlux from s
      SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                         + sum_aksgz*Smrho_s(mph)

      ! divDarcyFlux_r, r represent s' in paper, r is node
      do r = 1, NbNode_in_k
         numr = Nodebyk(r)

         ! P(mph) = PressionNode + PressionCap(mph)
         do j = 1, NbIncTotalPrim(IncNode(numr)%ic)
            divDarcyFlux_r(j, mph, r) = divDarcyFlux_r(j, mph, r) &
                                        - TkLocal_Darcy_ks(r)*divPressionNode(j, numr) &        ! a_{ks}^{s'} -P_s'
                                        - TkLocal_Darcy_ks(r)*divPhasePressureNode(j, mph, numr)     ! a_{ks}^{s'} -PressionCap_s'
         end do

         ! SmDarcyFlux from r (node)
         SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                            - TkLocal_Darcy_ks(r)*SmPressionNode(numr)    ! -a_{ks}^{s'} * Sm

      end do ! end of r node

      ! divDarcyFlux_r, r represent s' in paper, r is frac
      do r = 1, NbFrac_in_k
         numr = Fracbyk(r) ! numr is frac num here
         rf = r + NbNode_in_k

         do j = 1, NbIncTotalPrim(IncFrac(numr)%ic)
            divDarcyFlux_r(j, mph, rf) = divDarcyFlux_r(j, mph, rf) &
                                         - TkLocal_Darcy_ks(rf)*divPressionFrac(j, numr) &        ! a_{ks}^{s'} -P_s'
                                         - TkLocal_Darcy_ks(rf)*divPhasePressureFrac(j, mph, numr)     ! a_{ks}^{s'} -PressionCap_s'
         end do

         ! SmDarcyFlux from r (frac)
         SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                            - TkLocal_Darcy_ks(rf)*SmPressionFrac(numr)    ! -a_{ks}^{s'} * Sm

      end do ! end of r frac

   end subroutine Jacobian_divDarcyFlux_mph

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES

   subroutine Jacobian_divrho_gravity( &
      X1, dSdXp1, rho1, drhodXp1, rhsrho1, davrhodXp1, rhsavrho1, &
      X2, dSdXp2, rho2, drhodXp2, rhsrho2, davrhodXp2, rhsavrho2)
      type(TYPE_IncCVReservoir), intent(in) :: X1, X2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: dSdXp1, dSdXp2
      double precision, dimension(NbPhase), intent(in) :: rho1, rho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: drhodXp1, drhodXp2
      double precision, dimension(NbPhase), intent(in) :: rhsrho1, rhsrho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(out) :: davrhodXp1, davrhodXp2
      double precision, dimension(NbPhase), intent(out) :: rhsavrho1, rhsavrho2

      integer :: m, mph, n(NbPhase), nmax

      davrhodXp1 = 0.d0
      davrhodXp2 = 0.d0
      rhsavrho1 = 0.d0
      rhsavrho2 = 0.d0
      n = 0

      do m = 1, NbPhasePresente_ctx(X1%ic)
         mph = NumPhasePresente_ctx(m, X1%ic)
         n(mph) = n(mph) + 1
         davrhodXp1(:, mph) = drhodXp1(:, mph)
         rhsavrho1(mph) = rhsrho1(mph)
      end do

      do m = 1, NbPhasePresente_ctx(X2%ic)
         mph = NumPhasePresente_ctx(m, X2%ic)
         n(mph) = n(mph) + 1
         davrhodXp2(:, mph) = drhodXp2(:, mph)
         rhsavrho2(mph) = rhsrho2(mph)
      end do

      do m = 1, NbPhase
         nmax = max(n(m), 1)
         davrhodXp1(:, m) = davrhodXp1(:, m)/nmax
         davrhodXp2(:, m) = davrhodXp2(:, m)/nmax
         rhsavrho1(m) = rhsavrho1(m)/nmax
         rhsavrho2(m) = rhsavrho2(m)/nmax
      enddo

   end subroutine Jacobian_divrho_gravity

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_EXTRAPOLATE

   subroutine Jacobian_divrho_gravity( &
      X1, dSdXp1, rho1, drhodXp1, rhsrho1, davrhodXp1, rhsavrho1, &
      X2, dSdXp2, rho2, drhodXp2, rhsrho2, davrhodXp2, rhsavrho2)
      type(TYPE_IncCVReservoir), intent(in) :: X1, X2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: dSdXp1, dSdXp2
      double precision, dimension(NbPhase), intent(in) :: rho1, rho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: drhodXp1, drhodXp2
      double precision, dimension(NbPhase), intent(in) :: rhsrho1, rhsrho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(out) :: davrhodXp1, davrhodXp2
      double precision, dimension(NbPhase), intent(out) :: rhsavrho1, rhsavrho2

      integer :: k

      davrhodXp1 = 0.d0
      davrhodXp2 = 0.d0
      rhsavrho1 = 0.d0
      rhsavrho2 = 0.d0

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               davrhodXp1(:, k) = 0.5*drhodXp1(:, k)
               davrhodXp2(:, k) = 0.5*drhodXp2(:, k)
               rhsavrho1(k) = 0.5*rhsrho1(k)
               rhsavrho2(k) = 0.5*rhsrho2(k)
            else
               ! FIXME: how to use phase pressure?
               !    call f_VolumetricMassDensity_with_derivatives(k, X2%Pression, X2%Temperature, X2%Comp, rhoext, dP, dT, dC)
               davrhodXp1(:, k) = 0.5*drhodXp1(:, k)
               !   davrhodXp2(:, k) = 0.5*drhodXp2(:, k)
               rhsavrho1(k) = 0.5*rhsrho1(k)
               !   rhsavrho2(k) = 0.5*rhsrho2(k)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) then
               ! FIXME: how to use phase pressure?
               !    call f_VolumetricMassDensity_with_derivatives(k, X1%Pression, X1%Temperature, X1%Comp, rhoext, dP, dT, dC)
               ! davrhodXp1(:, k) = 0.5*drhodXp1(:, k)
               davrhodXp2(:, k) = 0.5*drhodXp2(:, k)
               !   rhsavrho1(k) = 0.5*rhsrho1(k)
               rhsavrho2(k) = 0.5*rhsrho2(k)
            endif
         endif
      end do

   end subroutine Jacobian_divrho_gravity

! ComPASS_GRAVITY_AVERAGE_RHO_EXTRAPOLATE
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_MIN

   subroutine Jacobian_divrho_gravity( &
      X1, dSdXp1, rho1, drhodXp1, rhsrho1, davrhodXp1, rhsavrho1, &
      X2, dSdXp2, rho2, drhodXp2, rhsrho2, davrhodXp2, rhsavrho2)
      type(TYPE_IncCVReservoir), intent(in) :: X1, X2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: dSdXp1, dSdXp2
      double precision, dimension(NbPhase), intent(in) :: rho1, rho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: drhodXp1, drhodXp2
      double precision, dimension(NbPhase), intent(in) :: rhsrho1, rhsrho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(out) :: davrhodXp1, davrhodXp2
      double precision, dimension(NbPhase), intent(out) :: rhsavrho1, rhsavrho2

      integer :: k

      davrhodXp1 = 0.d0
      davrhodXp2 = 0.d0
      rhsavrho1 = 0.d0
      rhsavrho2 = 0.d0

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               if (rho1(k) < rho2(k)) then
                  davrhodXp1(:, k) = drhodXp1(:, k)
                  rhsavrho1(k) = rhsrho1(k)
               else
                  davrhodXp2(:, k) = drhodXp2(:, k)
                  rhsavrho2(k) = rhsrho2(k)
               endif
            else
               davrhodXp1(:, k) = drhodXp1(:, k)
               rhsavrho1(k) = rhsrho1(k)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) then
               davrhodXp2(:, k) = drhodXp2(:, k)
               rhsavrho2(k) = rhsrho2(k)
            endif
         end if
      end do

   end subroutine Jacobian_divrho_gravity

! ComPASS_GRAVITY_AVERAGE_RHO_MIN
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY

   subroutine Jacobian_divrho_gravity( &
      X1, dSdXp1, rho1, drhodXp1, rhsrho1, davrhodXp1, rhsavrho1, &
      X2, dSdXp2, rho2, drhodXp2, rhsrho2, davrhodXp2, rhsavrho2)
      type(TYPE_IncCVReservoir), intent(in) :: X1, X2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: dSdXp1, dSdXp2
      double precision, dimension(NbPhase), intent(in) :: rho1, rho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: drhodXp1, drhodXp2
      double precision, dimension(NbPhase), intent(in) :: rhsrho1, rhsrho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(out) :: davrhodXp1, davrhodXp2
      double precision, dimension(NbPhase), intent(out) :: rhsavrho1, rhsavrho2

      integer :: k
      double precision :: Stot

      davrhodXp1 = 0.d0
      davrhodXp2 = 0.d0
      rhsavrho1 = 0.d0
      rhsavrho2 = 0.d0

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               Stot = X1%Saturation(k) + X2%Saturation(k)
               if (Stot > epsilon_avrho) then
                  davrhodXp1(:, k) = (dSdXp1(:, k)*rho1(k) + X1%Saturation(k)*drhodXp1(:, k))/Stot
                  rhsavrho1(k) = X1%Saturation(k)*rhsrho1(k)/Stot
                  davrhodXp2(:, k) = (dSdXp2(:, k)*rho2(k) + X2%Saturation(k)*drhodXp2(:, k))/Stot
                  rhsavrho2(k) = X2%Saturation(k)*rhsrho2(k)/Stot
               else
                  davrhodXp1(:, k) = 0.5d0*drhodXp1(:, k)
                  rhsavrho1(k) = 0.5d0*rhsrho1(k)
                  davrhodXp2(:, k) = 0.5d0*drhodXp2(:, k)
                  rhsavrho2(k) = 0.5d0*rhsrho2(k)
               end if
            else
               davrhodXp1(:, k) = drhodXp1(:, k)
               rhsavrho1(k) = rhsrho1(k)
            end if
         else
            if (phase_can_be_present(k, X2%ic)) then
               davrhodXp2(:, k) = drhodXp2(:, k)
               rhsavrho2(k) = rhsrho2(k)
            endif
         endif
      end do

   end subroutine Jacobian_divrho_gravity

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS

   subroutine Jacobian_divrho_gravity( &
      X1, dSdXp1, rho1, drhodXp1, rhsrho1, davrhodXp1, rhsavrho1, &
      X2, dSdXp2, rho2, drhodXp2, rhsrho2, davrhodXp2, rhsavrho2)
      type(TYPE_IncCVReservoir), intent(in) :: X1, X2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: dSdXp1, dSdXp2
      double precision, dimension(NbPhase), intent(in) :: rho1, rho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(in) :: drhodXp1, drhodXp2
      double precision, dimension(NbPhase), intent(in) :: rhsrho1, rhsrho2
      double precision, dimension(NbIncTotalPrimMax, NbPhase), intent(out) :: davrhodXp1, davrhodXp2
      double precision, dimension(NbPhase), intent(out) :: rhsavrho1, rhsavrho2

      integer :: k
      double precision :: chi
      double precision, dimension(NbIncTotalPrimMax) :: dchidXp1, dchidXp2

      davrhodXp1 = 0.d0
      davrhodXp2 = 0.d0
      rhsavrho1 = 0.d0
      rhsavrho2 = 0.d0

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               if (X1%Saturation(k) < X2%Saturation(k)) then
                  if (X2%Saturation(k) > 0.d0) then
                     chi = X1%Saturation(k)/X2%Saturation(k)
                     dchidXp1 = dSdXp1(:, k)/X2%Saturation(k)
                     dchidXp2 = -chi*dSdXp2(:, k)/X2%Saturation(k)
                     !  rho(k) = (chi*rho1(k) + rho2(k))/(1 + chi)
                     !  davrhodX(:, k) = (dchidX*rho1(k) + chi * drhodX(:, k)) / (1 + chi) - dchidX * (chi*rho1(k) + rho2(k)) / (1 + chi)**2
                     !  davrhodX(:, k) = (dchidX*rho1(k)*(1+chi) - dchidX * (chi*rho1(k) + rho2(k)) / (1 + chi)**2  + (chi  / (1 + chi)) * drhodX(:, k)
                     davrhodXp1(:, k) = dchidXp1*(rho1(k) - rho2(k))/(1 + chi)**2 + (chi/(1 + chi))*drhodXp1(:, k)
                     davrhodXp2(:, k) = dchidXp2*(rho1(k) - rho2(k))/(1 + chi)**2 + (chi/(1 + chi))*drhodXp2(:, k)
                     rhsavrho1(k) = (chi*rhsrho1(k))/(1 + chi)
                     rhsavrho2(k) = rhsrho2(k)/(1 + chi)
                  endif
               else
                  if (X1%Saturation(k) > 0.d0) then
                     chi = X2%Saturation(k)/X1%Saturation(k)
                     dchidXp1 = -chi*dSdXp1(:, k)/X1%Saturation(k)
                     dchidXp2 = dSdXp2(:, k)/X1%Saturation(k)
                     davrhodXp1(:, k) = dchidXp1*(rho2(k) - rho1(k))/(1 + chi)**2 + (chi/(1 + chi))*drhodXp1(:, k)
                     davrhodXp2(:, k) = dchidXp2*(rho2(k) - rho1(k))/(1 + chi)**2 + (chi/(1 + chi))*drhodXp2(:, k)
                     rhsavrho1(k) = rhsrho1(k)/(1 + chi)
                     rhsavrho2(k) = (chi*rhsrho2(k))/(1 + chi)
                  endif
               endif
            else
               davrhodXp1(:, k) = drhodXp1(:, k)
               rhsavrho1(k) = rhsrho1(k)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) then
               davrhodXp2(:, k) = drhodXp2(:, k)
               rhsavrho2(k) = rhsrho2(k)
            endif
         endif
      end do

   end subroutine Jacobian_divrho_gravity

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS
#endif

   pure subroutine Jacobian_divFourierFlux( &
      NbIncTotalPrim, IncNode, IncFrac, & ! does not depend on k or s, FaceToFracLocal is missing
      ic_k, Nodebyk, Fracbyk, &
      TkLocal_Fourier_ks, &
      divTemperature_k, SmTemperature_k, &
      divTemperatureNode, SmTemperatureNode, & ! does not depend on k or s
      divTemperatureFrac, SmTemperatureFrac, & ! does not depend on k or s
      divFourierFlux_k, & ! sum_{s'} a_{k,s}^s' T_k
      divFourierFlux_r, & ! sum_{s'} a_{k,s}^s' -T_s'
      SmFourierFlux)

      integer, intent(in) :: &
         NbIncTotalPrim(NbContexte)
      type(TYPE_IncCVReservoir), dimension(:), intent(in) :: &
         IncNode, IncFrac
      integer(c_int), intent(in) :: ic_k
      integer(c_int), dimension(:), intent(in) :: Nodebyk, Fracbyk

      double precision, dimension(:, :), intent(in) :: &
         divTemperatureNode, divTemperatureFrac
      double precision, dimension(:), intent(in) :: &
         SmTemperatureNode, SmTemperatureFrac, &
         TkLocal_Fourier_ks
      double precision, intent(in) :: &
         divTemperature_k(NbIncTotalPrimMax), &
         SmTemperature_k

      double precision, intent(out) :: &
         divFourierFlux_k(NbIncTotalPrimMax), &
         SmFourierFlux
      double precision, dimension(:, :), intent(out) :: divFourierFlux_r ! r represent s' in paper

      integer :: j, r, numr, rf
      integer :: NbNode_in_k, NbFrac_in_k

      double precision :: sum_aks

      divFourierFlux_k(:) = 0.d0
      divFourierFlux_r(:, :) = 0.d0
      SmFourierFlux = 0.d0

      ! number of nodes/fracs in k
      NbNode_in_k = size(Nodebyk)
      NbFrac_in_k = size(Fracbyk)

      ! sum_aks = \sum_{r \in dof(k)} a_{k,s}^r
      sum_aks = 0.d0
      do r = 1, NbNode_in_k + NbFrac_in_k
         sum_aks = sum_aks + TkLocal_Fourier_ks(r)
      end do

      ! divFourierfFlux_k
      do j = 1, NbIncTotalPrim(ic_k)
         divFourierFlux_k(j) = sum_aks*divTemperature_k(j)
      end do

      SmFourierFlux = sum_aks*SmTemperature_k

      ! divFourierFlux_r, r represent s' in paper, r is node
      do r = 1, NbNode_in_k
         numr = Nodebyk(r)

         do j = 1, NbIncTotalPrim(IncNode(numr)%ic)
            divFourierFlux_r(j, r) = -TkLocal_Fourier_ks(r)*divTemperatureNode(j, numr)
         end do

         SmFourierFlux = SmFourierFlux &
                         - TkLocal_Fourier_ks(r)*SmTemperatureNode(numr)
      end do

      ! divFourierFlux_r, r represent s' in paper, r is frac
      do r = 1, NbFrac_in_k
         numr = Fracbyk(r) ! numr is frac num here
         rf = r + NbNode_in_k

         do j = 1, NbIncTotalPrim(IncFrac(numr)%ic)
            divFourierFlux_r(j, rf) = -TkLocal_Fourier_ks(rf)*divTemperatureFrac(j, numr)
         end do

         SmFourierFlux = SmFourierFlux &
                         - TkLocal_Fourier_ks(rf)*SmTemperatureFrac(numr)
      end do

   end subroutine Jacobian_divFourierFlux

   ! Regularization of Jacobian
   ! operation on JacBigA, bigSm
   subroutine Jacobian_Regularization() &
      bind(C, name="Jacobian_Regularization")

      integer :: k, rowk, colk

      ! rows of node own
      do k = 1, NbNodeOwn_Ncpus(commRank + 1)

         rowk = k ! row of k
         colk = rowk ! col of diag element in rowk

         ! regularization for node
         call Jacobian_Regularization_row(k, rowk, colk, "n")
      end do

      ! rows of frac own
      do k = 1, NbFracOwn_Ncpus(commRank + 1)

         rowk = k + NbNodeOwn_Ncpus(commRank + 1) ! row of k
         colk = k + NbNodeLocal_Ncpus(commRank + 1) ! col of diag element in rowk

         ! regularization for frac
         call Jacobian_Regularization_row(k, rowk, colk, "f")
      end do

      ! rows of cell
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         rowk = k + NbNodeOwn_Ncpus(commRank + 1) &
                + NbFracOwn_Ncpus(commRank + 1) ! row of k

         colk = k + NbNodeLocal_Ncpus(commRank + 1) &
                + NbFracLocal_Ncpus(commRank + 1)! col of diag element in rowk

         ! regularization for cell
         call Jacobian_Regularization_row(k, rowk, colk, "c")
      end do

   end subroutine Jacobian_Regularization

   ! sub subroutine of Jacobian_Regularization
   ! used for regularization for node/frac/cell
   subroutine Jacobian_Regularization_row(k, rowk, colk, cv)

      ! row and col of diag block element in JacBigA
      integer, intent(in) :: k, rowk, colk
      character, intent(in) :: cv

      integer :: i, icp, iph, j, nz, nzj
      integer :: errcode, Ierr
      double precision :: sumcol

      ! look for the diag: JacBigA(:,:,nz)
      do i = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
         if (JacBigA%Num(i) == colk) then
            nz = i
            exit
         end if
      end do

      do i = 1, NbCompThermique

         ! look for col nul in JacBigA(:,:,nz)
         ! ps. index order of matrix JacBigA(:,:,nz) is (col,row)
         sumcol = 0.d0 ! sum of col i
         do j = 1, NbCompThermique
            sumcol = sumcol + abs(JacBigA%Val(i, j, nz))
         end do

         ! warning
         ! TODO: CHECKME: the following is a magic number
         ! eps*1e-3 for small permeability in the pressure equation
         ! TODO: find an automatic scaling
         if (sumcol < eps*1e-3) then ! col i is null

            ! look for component C_{i}^alpha corresponding to the col i

            if (cv .eq. 'n') then ! node
               j = NumIncTotalPrimNode(i, k)
               icp = 0
               iph = 0

               if (j <= NbIncPTC_ctx(IncNode(k)%ic)) then
                  icp = NumIncPTC2NumIncComp_comp_ctx(j, IncNode(k)%ic)
                  iph = NumIncPTC2NumIncComp_phase_ctx(j, IncNode(k)%ic)
               endif
            else if (cv .eq. 'f') then ! fracs
               j = NumIncTotalPrimFrac(i, k)
               icp = NumIncPTC2NumIncComp_comp_ctx(j, IncFrac(k)%ic)
               iph = NumIncPTC2NumIncComp_phase_ctx(j, IncFrac(k)%ic)
            else if (cv .eq. 'c') then ! cell
               j = NumIncTotalPrimCell(i, k)
               icp = NumIncPTC2NumIncComp_comp_ctx(j, IncCell(k)%ic)
               iph = NumIncPTC2NumIncComp_phase_ctx(j, IncCell(k)%ic)
            else
               print *, ""
               print *, "Regularization error: cv should be node/frac/cell"
               errcode = 51
               call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
            end if

            if (icp .eq. 0) then
               print *, " "
               ! print*, JacBigA%Val(:,1,nz)
               ! print*, JacBigA%Val(:,2,nz)
               print *, "Regularization error: icp=0"
               write (*, '(A,4I8)') "    row/col/i/commRank=", rowk, colk, i, commRank
               errcode = 52
               if (cv .eq. 'n') then ! node
                  write (*, *) "j ", j, "cv ", cv, " ic ", IncNode(k)%ic
               else if (cv .eq. 'f') then ! fracs
                  write (*, *) "j ", j, "cv ", cv, " ic ", IncFrac(k)%ic
               else if (cv .eq. 'c') then ! cell
                  write (*, *) "j ", j, "cv ", cv, " ic ", IncCell(k)%ic
               endif
               write (*, *) "sumcol ", sumcol
               write (*, *) "i", i
               ! write(*,*) "x ", XFaceLocal(:,FracToFaceLocal(k))

               ! write(*,*) "Pression ", IncFrac(k)%Pression
               ! write(*,*) "Cwater: ", IncFrac(k)%Comp(:,1)
               ! write(*,*) "Coil  : ", IncFrac(k)%Comp(:,2)
               ! write(*,*) "Satu  : ", IncFrac(k)%Saturation(:)

               call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
            endif

            ! row icp in JacBigA is set to be zero
            do j = JacBigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1)
               nzj = JacBigA%Num(j)
               JacBigA%Val(:, icp, nzj) = 0.d0
            end do

            ! JacBigA(i,icp,nz)=1
            JacBigA%Val(i, icp, nz) = 1.d0

            ! bigSm(i,rowk) = 0
            bigSm(i, rowk) = 0.d0
         end if
      end do

   end subroutine Jacobian_Regularization_row

   ! Alignment of Jacobian: diag method
   ! operation on JacA, Sm
   subroutine Jacobian_Alignment_diag() &
      bind(C, name="Jacobian_Alignment_diag")

      integer :: k, rowk, colk

      if (NbMSWellLocal > 0) then
#ifndef NDEBUG
         call CommonMPI_abort( &
            "Jacobian: Diagonal alignment is not supported for mswells")
#endif
      end if

      ! rows of node own
      do k = 1, NbNodeOwn_Ncpus(commRank + 1)

         rowk = k ! row of k
         colk = rowk ! col of diag element in rowk

         call Jacobian_Alignment_diag_row(rowk, colk)
      end do

      ! rows of frac own
      do k = 1, NbFracOwn_Ncpus(commRank + 1)

         rowk = k + NbNodeOwn_Ncpus(commRank + 1)   ! row of k
         colk = k + NbNodeLocal_Ncpus(commRank + 1) ! col of diag element in rowk

         call Jacobian_Alignment_diag_row(rowk, colk)
      end do

   end subroutine Jacobian_Alignment_diag

   ! sub subroutine of Jacobian_Alignment: diag method
   ! used for Alignment for node/frac
   subroutine Jacobian_Alignment_diag_row(rowk, colk)

      integer, intent(in) :: rowk ! local row
      integer, intent(in) :: colk ! global row

      integer, parameter :: n = NbCompThermique
      ! optimal size of the WORK array, pre-queried for machine
      integer, parameter :: lwork = 320 ! Mac
      double precision, dimension(lwork) :: work
      integer :: ipival(n), info
      integer :: i, nz, errcode, Ierr

      double precision, dimension(n, n) :: AA, BB
      double precision, dimension(n) :: Smk

      nz = JacA%Pt(rowk) + findloc(JacA%Num((JacA%Pt(rowk) + 1):JacA%Pt(rowk + 1)), colk, 1)
      BB = JacA%Val(:, :, nz)

      ! BB = inv(JacA%Val(:,:,nz))
      ! ps. the index order of JacA%Val(:,:,nz) is (col, row)
      ! so the index order of BB is also (col, row)
      call dgetrf(n, n, BB, n, ipival, info)
      if (info /= 0) then
         print *, "dgetrf error", info, "in Alignment, rowk/colk = ", rowk, colk
         print *, "shape of BB", shape(BB)
         print *, BB
         print *, "Try to locate error..."
         call Jacobian_JacBigA_locate_frac_row(rowk)
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      call dgetri(n, BB, n, ipival, work, lwork, info)
      if (info /= 0) then
         print *, "dgetri error", info, "in Alignment, rowk/colk = ", rowk, colk
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      ! JacA%Val(:,:,i) = JacA%Val(:,:,i)*BB, row k
      do i = JacA%Pt(rowk) + 1, JacA%Pt(rowk + 1)
         AA(:, :) = JacA%Val(:, :, i)
         call dgemm('N', 'N', n, n, n, 1.d0, AA, n, BB, n, 0.d0, JacA%Val(:, :, i), n)
      end do

      ! Sm(:,rowk) = BB * Sm(:,rowk), rowk
      ! transpose of BB is necessary since the index of BB is (col, row)
      Smk(:) = Sm(:, rowk)
      call dgemv('T', n, n, 1.d0, BB, n, Smk, 1, 0.d0, Sm(:, rowk), 1)

   end subroutine Jacobian_Alignment_diag_row

   ! Alignment of Jacobian: manually method
   ! operation on JacA, Sm
   subroutine Jacobian_Alignment_man() &
      bind(C, name="Jacobian_Alignment_man")

      integer :: k, rowk
      integer :: glb_unk_idx_s, unk_idx_s, s, num_s, row_s

      ! rows of node own
      do k = 1, NbNodeOwn_Ncpus(commRank + 1)

         ! WARNING: l'alignement avec alignemat ne s'applique qu'aux eqs de conservation
         ! donc pas aux noeuds DIR DIR
         ! TODO: le cas DIR Neu ou Neu Dir reste a faire

         if ((IdNodeLocal(k)%P .ne. "d") .and. (IdNodeLocal(k)%T .ne. "d")) then
            call Jacobian_Alignment_man_row(k, IncNode(k)%ic)
         else if ((IdNodeLocal(k)%P .eq. "d") .and. (IdNodeLocal(k)%T .ne. "d")) then
            call CommonMPI_abort('in manual alignment Jacobian mix Dir/Neu node are not implemented')
         else if ((IdNodeLocal(k)%T .eq. "d") .and. (IdNodeLocal(k)%P .ne. "d")) then
            call CommonMPI_abort('in manual alignment Jacobian mix Dir/Neu node are not implemented')
         endif
      end do

      ! rows of frac own
      do k = 1, NbFracOwn_Ncpus(commRank + 1)

         rowk = k + NbNodeOwn_Ncpus(commRank + 1) ! row of k
         call Jacobian_Alignment_man_row(rowk, IncFrac(k)%ic)
      end do

      do k = 1, NbMSWellLocal
         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s)
            !Unkown indices of node s
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s) !row of s @ MSWells_JacA
            glb_unk_idx_s = unk_idx_s + mswells_JacA_row_off  !unkown idx @ JacA
            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            row_s = glb_unk_idx_s! row of s @ JacA
            call Jacobian_Alignment_man_row(row_s, IncMSWell(s)%coats%ic) !k is not used
         end do
      end do

   end subroutine Jacobian_Alignment_man

   ! sub subroutine of Jacobian_Alignment: manually method
   ! used for Alignment for node/frac
   subroutine Jacobian_Alignment_man_row(rowk, ic)

      integer, intent(in) :: rowk, ic
      integer :: i

      double precision, dimension(NbCompThermique, NbCompThermique) :: &
         AA, BB
      double precision, dimension(NbCompThermique) :: &
         Smk

      ! the index order of JacA%Val(:,:,nz) is (col, row)
      ! the index order of aligmethod(:,:,ic) is also (col, row)

      BB(:, :) = aligmat(:, :, ic)

      ! JacA%Val(:,:,i) = JacA%Val(:,:,i) * aligmethod(:,:,ic)
      ! since all the matrix are transpose
      do i = JacA%Pt(rowk) + 1, JacA%Pt(rowk + 1)

         AA(:, :) = JacA%Val(:, :, i)
         call dgemm('N', 'N', NbCompThermique, NbCompThermique, NbCompThermique, &
                    1.d0, AA, NbCompThermique, BB, NbCompThermique, 0.d0, &
                    JacA%Val(:, :, i), NbCompThermique)
      end do

      ! Sm(:,rowk) = BB * Sm(:,rowk), rowk
      ! transpose of aligmethod(:,:,ic) is necessary
      ! since the index order of BB is (col, row)
      Smk(:) = Sm(:, rowk)
      call dgemv('T', NbCompThermique, NbCompThermique, 1.d0, &
                 BB, NbCompThermique, Smk, 1, &
                 0.d0, Sm(:, rowk), 1)

   end subroutine Jacobian_Alignment_man_row

   !---------------------------------------------------------------------------
   !> @brief
   !> Compute the in place LU factorization of \f$J_{KK}\f$ diagonal blocks
   !---------------------------------------------------------------------------
   subroutine Jacobian_Schur_Jkk_LU_factorization(Mat)
      type(CSRArray2dble), intent(out) :: Mat

      integer, parameter :: n = NbCompThermique ! The size of squared block
      integer :: k, nb_cells, row_k, row_cells_offset
      integer :: info
      real(c_double), pointer, dimension(:, :) :: JkkT
#ifndef NDEBUG
      integer :: col_cells_offset, col_cells_end
      integer, pointer, dimension(:) :: cols
#endif

      nb_cells = NbCellLocal_Ncpus(commRank + 1)

#ifndef NDEBUG
      col_cells_offset = NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1)
      col_cells_end = col_cells_offset + nb_cells
#endif

      if (.not. allocated(pivot)) then
         allocate (pivot(n, nb_cells))
      end if

      row_cells_offset = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)

      do k = 1, nb_cells

         row_k = row_cells_offset + k

#ifndef NDEBUG
         cols => Mat%Num(Mat%Pt(row_k) + 1:Mat%Pt(row_k + 1) - 1)
         if (any(cols > col_cells_offset)) &
            call CommonMPI_abort("Jacobian_Schur_Jkk_LU_factorization: inconsistent column in Mat")
         ! Check J_{KW}=0
         if ( &
            Mat%Num(Mat%Pt(row_k + 1)) <= col_cells_offset .or. &
            Mat%Num(Mat%Pt(row_k + 1)) > col_cells_end &
            ) call CommonMPI_abort("Jacobian_Schur_Jkk_LU_factorization: inconsistent cell column in Mat")
#endif

         ! J_{KK} is block diagonal and J_{KW}=0
         ! so Jkk is the last block of Mat%Val(:, :, row_begin:row_end) with:
         ! row_begin = Mat%Pt(cells_row_offset + k) + 1
         ! row_end = Mat%Pt(cells_row_offset + k + 1)
         JkkT => Mat%Val(:, :, Mat%Pt(row_k + 1))
         call dgetrf(n, n, JkkT, n, pivot(:, k), info)

         if (info /= 0) then
            write (*, *) "Local block:", JkkT
            call CommonMPI_abort("dgetrf error in Jacobian_Schur_Jkk_LU_factorization")
         end if

      end do

   end subroutine Jacobian_Schur_Jkk_LU_factorization

   !---------------------------------------------------------------------------
   !> @brief
   !> Extract \f$J_{\nu\nu}\f$ \f$J_{\nuW}\f$ \f$S_{\nu}\f$...
   !>
   !> Extract \f[
   !> \left[\begin{array}{cc}
   !>    J_{\nu\nu} & J_{\nu W}\\
   !>    J_{W\nu} & J_{WW}
   !> \end{array}\right]
   !> \f]
   !---------------------------------------------------------------------------
   subroutine Jacobian_Schur_extract_submatrix(BigMat, Mat)
      type(CSRArray2dble), intent(in) :: BigMat
      type(CSRArray2dble), intent(out) :: Mat

      integer :: nu, mu, nb_rows
      integer :: nc, nb_cells, row_cells_offset, row_cells_end
      integer :: j, colj, col_cells_offset, col_cells_end
      integer, pointer, dimension(:) :: col

      nb_cells = NbCellLocal_Ncpus(commRank + 1)
      row_cells_offset = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)
      row_cells_end = row_cells_offset + nb_cells
      col_cells_offset = NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1)
      col_cells_end = col_cells_offset + nb_cells
      nb_rows = size(BigSm, 2)

#ifndef NDEBUG
      if (nb_rows /= (row_cells_end + NbWellInjOwn_Ncpus(commRank + 1) + NbWellProdOwn_Ncpus(commRank + 1) + NbMSWellNodeOwn)) &
         call CommonMPI_abort("Unconsistent number of rows in Jacobian_Schur_extract_submatrix!")
#endif

      ! We extract J_{\V\V} J_{\V\W} skipping J_{\V\K}
      do nu = 1, nb_rows
         mu = nu
         if (mu > row_cells_offset) then
            if (mu <= row_cells_end) cycle
            mu = mu - nb_cells
         end if
         nc = 0 ! number of cells on row nu
         col => BigMat%Num(BigMat%Pt(nu) + 1:BigMat%Pt(nu + 1))
         do j = 1, size(col)
            colj = col(j)
            if (colj <= col_cells_offset) then
               Mat%Val(:, :, Mat%Pt(mu) + j) = BigMat%Val(:, :, BigMat%Pt(nu) + j)
            else if (colj <= col_cells_end) then
               nc = nc + 1
            else
               Mat%Val(:, :, Mat%Pt(mu) + j - nc) = BigMat%Val(:, :, BigMat%Pt(nu) + j)
            end if
         end do
      end do

   end subroutine Jacobian_Schur_extract_submatrix

   !---------------------------------------------------------------------------
   !> @brief
   !> Extract \f$J_{\nu\nu}\f$ \f$J_{\nuW}\f$ \f$S_{\nu}\f$...
   !>
   !> Extract \f[
   !> \left[\begin{array}{c}
   !>    S_{\nu} \\ S_{W}
   !> \end{array}\right]
   !> \f]
   !---------------------------------------------------------------------------
   subroutine Jacobian_Schur_extract_subvector(BigVec, Vec)
      real(c_double), dimension(:, :), intent(in) :: BigVec
      real(c_double), dimension(:, :), intent(out) :: Vec
      integer :: nb_rows, row_cells_offset, row_cells_end

      row_cells_offset = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)
      row_cells_end = row_cells_offset + NbCellLocal_Ncpus(commRank + 1)
      nb_rows = size(BigVec, 2)

#ifndef NDEBUG
      if (nb_rows /= (row_cells_offset + NbCellLocal_Ncpus(commRank + 1) &
                      + NbWellInjOwn_Ncpus(commRank + 1) + NbWellProdOwn_Ncpus(commRank + 1) + NbMSWellNodeOwn)) &
         call CommonMPI_abort("Unconsistent number of rows in Jacobian_Schur_init_Sm!")
#endif

      Vec(:, 1:row_cells_offset) = BigVec(:, 1:row_cells_offset)
      if (size(Vec, 2) > row_cells_offset) &
         Vec(:, (row_cells_offset + 1):size(Vec, 2)) = BigVec(:, (row_cells_end + 1):nb_rows)

   end subroutine Jacobian_Schur_extract_subvector

   ! Schur complement
   ! cf. latex file in docs/tech/schur.tex
   subroutine Jacobian_Schur_substitution(BigMat, Mat, BigVec, Vec, RHS_only)
      type(CSRArray2dble), intent(in) :: BigMat
      type(CSRArray2dble), intent(inout) :: Mat
      real(c_double), dimension(:, :), intent(in) :: BigVec
      real(c_double), dimension(:, :), intent(inout) :: Vec
      logical, intent(in) :: RHS_only

      integer, parameter :: n = NbCompThermique ! The size of squared block
      integer :: i, j, k, l, col_k, row_k, nu, nup, info
      integer :: row_nodes_end, row_cells_offset, col_cells_offset, col_cells_end
      real(c_double) :: B(n, n)
      real(c_double), pointer, dimension(:, :) :: JkkT_LU
      integer, pointer, dimension(:) :: colsJnu, colsJk

      col_cells_offset = NbNodeLocal_Ncpus(commRank + 1) + NbFracLocal_Ncpus(commRank + 1)
      col_cells_end = col_cells_offset + NbCellLocal_Ncpus(commRank + 1)
      row_nodes_end = NbNodeOwn_Ncpus(commRank + 1)
      row_cells_offset = row_nodes_end + NbFracOwn_Ncpus(commRank + 1)

      do nu = 1, row_cells_offset

         ! skip Dirichlet nodes
         if (nu <= row_nodes_end) then
            if ( &
#ifdef _THERMIQUE_
               (IdNodeLocal(nu)%P == "d") .and. (IdNodeLocal(nu)%T == "d") &
#else
               IdNodeLocal(nu)%P == "d" &
#endif
               ) cycle
         end if

         colsJnu => BigMat%Num(BigMat%Pt(nu) + 1:BigMat%Pt(nu + 1))
         do i = 1, size(colsJnu)

            col_k = colsJnu(i)
            ! process only columns that correspond to cells
            if (col_k <= col_cells_offset) cycle
            if (col_k > col_cells_end) exit
            k = col_k - col_cells_offset

            row_k = row_cells_offset + k

            ! B <- J_{\nu k}^T
            B = BigMat%Val(:, :, BigMat%Pt(nu) + i)
            ! J_{KK} is block diagonal and J_{KW}=0
            ! so Jkk is the last block of BigMat%Val(:, :, row_begin:row_end) with:
            ! row_begin = BigMat%Pt(row_k) + 1
            ! row_end = BigMat%Pt(row_k + 1)
            JkkT_LU => BigMat%Val(:, :, BigMat%Pt(row_k + 1))
            ! B  <- (J_{kk}^T)^{-1} B = (J_{kk}^T)^{-1} J_{\nu k}^T = B_{\nu k}
            call dgetrs('N', n, n, JkkT_LU, n, pivot(:, k), B, n, info)
            if (info /= 0) &
               call CommonMPI_abort("dgetrs error in Jacobian_Schur_substitution")

            ! S_{\nu} <- S_{\nu} - B^{T} S_{k}
            call dgemv('T', n, n, -1.d0, B, n, BigVec(:, row_k), 1, 1.d0, Vec(:, nu), 1)

            if (RHS_only) cycle

            colsJk => BigMat%Num(BigMat%Pt(row_k) + 1:BigMat%Pt(row_k + 1))
#ifndef NDEBUG
            ! Given the present problem we know that cells are not connected to wells
            ! and that cell block is (block-) diagonal
            if (size(colsJk) < 2) &
               call CommonMPI_abort("Jacobian_Schur_substitution: wrong connectivities!")
            if (colsJk(size(colsJk) - 1) > col_cells_offset) &
               call CommonMPI_abort("Jacobian_Schur_substitution: column should be node or fracture!")
            if (colsJk(size(colsJk)) /= col_cells_offset + k) &
               call CommonMPI_abort("Jacobian_Schur_substitution: column should be cell (on the diagonal)!")
#endif
            do j = 1, size(colsJk) - 1
               nup = colsJk(j)
               ! VAG connects nodes and fractures to all cell nodes so that findloc will not fail
               l = findloc(colsJnu, nup, 1)
#ifndef NDEBUG
               if (l == 0) &
                  call CommonMPI_abort( &
                  "Jacobian_Schur_substitution: could not find a column of BigJacA in JacA: check connectivities")
               if (Mat%Num(Mat%Pt(nu) + l) /= nup) &
                  call CommonMPI_abort("Jacobian_Schur_substitution: unconsistencies between Mat and BigMat")
#endif
               ! because of C like storage of the blocks we compute the transpose of
               ! J_{\nu\nu'}-\sum_{k\in\K}B_{\nu k}^{T}J_{k\nu'}
               call dgemm('N', 'N', n, n, n, -1.d0, &
                          BigMat%Val(:, :, BigMat%Pt(row_k) + j), n, B, n, 1.d0, &
                          Mat%Val(:, :, Mat%Pt(nu) + l), n)
            end do ! end of colsJk

         end do ! end of colsJnu

      end do ! end of \nu \in \V

   end subroutine Jacobian_Schur_substitution

   ! Schur complement
   ! documented in docs/tech/schur.tex
   subroutine Jacobian_Schur() &
      bind(C, name="Jacobian_Schur")

      ! in place LU factorization of JKK diagonal blocks
      call Jacobian_Schur_Jkk_LU_factorization(JacBigA)

      call Jacobian_Schur_extract_submatrix(JacBigA, JacA)
      call Jacobian_Schur_extract_subvector(BigSm, Sm)

      call Jacobian_Schur_substitution(JacBigA, JacA, BigSm, Sm, .false.)

   end subroutine Jacobian_Schur

   ! Non-zero sturcture of Jacobian
   subroutine Jacobian_StrucJacBigA

      !            node, frac, cell, wellinj, wellprod, mswell_nodes
      !           | a11, a12, a13, a14, a15 a16 | node own
      ! JacBigA = | a21, a22, a23, 0,   0   0   | frac own
      !           | a31, a32, a33, 0,   0   0   | cell (own and ghost)
      !           ! a41, 0,   0,   a44, 0   0   | wellinj own
      !           ! a51, 0,   0,   0,   a55 0   | wellprod own
      !           ! a61, 0,   0,   0,   0   a66 | mswell_nodes own

      !            node, frac, wellinj, wellprod, mswell_nodes
      ! JacA =    | A11, A12, A13, A14  A15     | node own
      !           | A21, A22, 0,   0    0       | frac own
      !           | A31, 0,   A33, 0    0       | wellinj own,
      !           | A41, 0,   0,   A44  0       | wellprod own
      !           ! A51, 0,   0,   0,   0   A55 | mswell_nodes own
      ! the nonzero structure of Aij is same as aij, i,j=1,2 in JacBigA
      ! A31=a31, A33=a33, A41=a41, A44=a44, A55=a55

      ! The non zero structure of aij is based on the connectivities
      ! a11: NodebyNodeOwn if not dirichlet
      ! a12: FracbyNodeOwn if not dir
      ! a13: CellbyNodeOwn if not dir
      ! a14: WellInjbyNodeOwn if not dir
      ! a15: WellProdbyNodeOwn if not dir
      ! Rq: if node own is Dir then only one non zero/row for a1*

      ! a21: NodebyFracOwn
      ! a22: FracbyFracOwn
      ! a23: CellbyFracOwn

      ! a31: NodebyCellLocal
      ! a32: FracbyCellLocal
      ! a33: diag block

      ! a41: NodebyWellInjLocal
      ! a44: diag
      ! a51: NodebyWellProdLocal
      ! a55: diag

      ! Four steps in this subroutine
      !   Number of non zeros each line: NbNnzbyline
      !   bigA%Nb, bigA%Pt using NbNnzbyline
      !   bigA%Num, non zero structure of bigA
      !   arrange bigA%Num such that in each row, the cols is in inscreasing order

      integer, dimension(:), allocatable :: &
         nbNnzbyline ! number of non zeros each line

      integer :: i, j, k, Nz, start, s, num_s, unk_idx_s, glb_unk_col_idx_s
      integer :: tmp
      integer :: &
         nbNodeOwn, nbFracOwn, nbWellInjOwn, nbWellProdOwn, &
         nbNodeLocal, nbFracLocal, nbCellLocal, nbWellInjLocal, nbWellProdLocal

      nbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)
      nbFracOwn = NbFracOwn_Ncpus(commRank + 1)
      nbWellInjOwn = NbWellInjOwn_Ncpus(commRank + 1)
      nbWellProdOwn = NbWellProdOwn_Ncpus(commRank + 1)

      nbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)
      nbFracLocal = NbFracLocal_Ncpus(commRank + 1)
      nbCellLocal = NbCellLocal_Ncpus(commRank + 1)
      nbWellInjLocal = NbWellInjLocal_Ncpus(commRank + 1)
      nbWellProdLocal = NbWellProdLocal_Ncpus(commRank + 1)

      JacBigA%Nb = nbNodeOwn + nbFracOwn + nbCellLocal + nbWellInjOwn + nbWellProdOwn
      JacBigA%Nb = JacBigA%Nb + NbMSWellNodeOwn
      mswells_JacBigA_unk_lcl_off = nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal + nbWellProdLocal

      allocate (nbNnzbyLine(JacBigA%Nb))
      allocate (JacBigA%Pt(JacBigA%Nb + 1))
      allocate (BigSm(NbCompThermique, JacBigA%Nb))

      ! a1* in JacbigA
      do i = 1, nbNodeOwn
         nbNnzbyLine(i) = NodebyNodeOwn%Pt(i + 1) - NodebyNodeOwn%Pt(i) &
                          + FracbyNodeOwn%Pt(i + 1) - FracbyNodeOwn%Pt(i) &
                          + CellbyNodeOwn%Pt(i + 1) - CellbyNodeOwn%Pt(i) &
                          + WellInjbyNodeOwn%Pt(i + 1) - WellInjbyNodeOwn%Pt(i) &
                          + WellProdbyNodeOwn%Pt(i + 1) - WellProdbyNodeOwn%Pt(i)
      end do

      ! a2* in JacbigA, rq: frac is not dir face
      start = nbNodeOwn
      do i = 1, nbFracOwn
         nbNnzbyLine(start + i) = NodebyFracOwn%Pt(i + 1) - NodebyFracOwn%Pt(i) &
                                  + FracbyFracOwn%Pt(i + 1) - FracbyFracOwn%Pt(i) &
                                  + CellbyFracOwn%Pt(i + 1) - CellbyFracOwn%Pt(i)
      end do

      ! a3* in JacbigA
      start = start + nbFracOwn
      do i = 1, nbCellLocal
         nbNnzbyLine(start + i) = &
            NodebyCellLocal%Pt(i + 1) - NodebyCellLocal%Pt(i) &
            + FracbyCellLocal%Pt(i + 1) - FracbyCellLocal%Pt(i) + 1
      end do

      ! a4* in JacbigA
      start = start + nbCellLocal
      do i = 1, nbWellInjOwn
         nbNnzbyLine(start + i) = &
            NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i) + 1 ! a14 and a44, (+1 sicne a44 is diag)
      end do

      ! a5* in JacbigA
      start = start + nbWellInjOwn
      ! print*, 'DEBUG HERE', NbWellProdOwn_Ncpus(commRank+1)
      do i = 1, nbWellProdOwn
         nbNnzbyLine(start + i) = &
            NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i) + 1 ! a55 is diag
      end do

      ! a6* in JacbigA
      start = start + nbWellProdOwn
      mswells_JacBigA_row_off = start
      ! print*, 'DEBUG HERE', NbMSWellOwn_Ncpus(commRank+1)
      do i = 1, NbMSWellNodeOwn! this is equal to MSWells_JacA%Nb
         nbNnzbyLine(start + i) = &
            MSWells_JacA%Pt(i + 1) - MSWells_JacA%Pt(i) + 1 !copy the nbcols of MSWells_JacA and sum the reservoir paired node
      end do

      ! JacBigA%Pt
      Nz = 0
      JacBigA%Pt(1) = 0
      do i = 1, JacBigA%Nb
         Nz = Nz + nbNnzbyLine(i)
         JacBigA%Pt(i + 1) = Nz
      end do

      deallocate (nbNnzbyLine)

      allocate (JacBigA%Num(Nz))
      allocate (JacBigA%Val(NbCompThermique, NbCompThermique, Nz))

      start = 0

      do i = 1, nbNodeOwn
         ! a11(i,:)
         do j = 1, NodebyNodeOwn%Pt(i + 1) - NodebyNodeOwn%Pt(i)
            JacBigA%Num(start + j) = NodebyNodeOwn%Num(j + NodebyNodeOwn%Pt(i))
         end do
         start = start + NodebyNodeOwn%Pt(i + 1) - NodebyNodeOwn%Pt(i)

         ! a12(i,:)
         do j = 1, FracbyNodeOwn%Pt(i + 1) - FracbyNodeOwn%Pt(i)
            JacBigA%Num(start + j) = FracbyNodeOwn%Num(j + FracbyNodeOwn%Pt(i)) + nbNodeLocal ! col
         end do
         start = start + FracbyNodeOwn%Pt(i + 1) - FracbyNodeOwn%Pt(i)

         ! a13(i,:)
         do j = 1, CellbyNodeOwn%Pt(i + 1) - CellbyNodeOwn%Pt(i)
            JacBigA%Num(start + j) = CellbyNodeOwn%Num(j + CellbyNodeOwn%Pt(i)) + nbNodeLocal + nbFracLocal ! col
         end do
         start = start + CellbyNodeOwn%Pt(i + 1) - CellbyNodeOwn%Pt(i)

         ! a14(i,:)
         do j = 1, WellInjbyNodeOwn%Pt(i + 1) - WellInjbyNodeOwn%Pt(i)
            JacBigA%Num(start + j) = WellInjbyNodeOwn%Num(j + WellInjbyNodeOwn%Pt(i)) &
                                     + nbNodeLocal + nbFracLocal + nbCellLocal
         end do
         start = start + WellInjbyNodeOwn%Pt(i + 1) - WellInjbyNodeOwn%Pt(i)

         ! a15(i;:)
         do j = 1, WellProdbyNodeOwn%Pt(i + 1) - WellProdbyNodeOwn%Pt(i)
            JacBigA%Num(start + j) = WellProdbyNodeOwn%Num(j + WellProdbyNodeOwn%Pt(i)) &
                                     + nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal
         end do
         start = start + WellProdbyNodeOwn%Pt(i + 1) - WellProdbyNodeOwn%Pt(i)
      end do

      do i = 1, nbFracOwn

         ! a21(i,:)
         do j = 1, NodebyFracOwn%Pt(i + 1) - NodebyFracOwn%Pt(i)
            JacBigA%Num(start + j) = NodebyFracOwn%Num(j + NodebyFracOwn%Pt(i))
         end do
         start = start + NodebyFracOwn%Pt(i + 1) - NodebyFracOwn%Pt(i)

         ! a22(i,:)
         do j = 1, FracbyFracOwn%Pt(i + 1) - FracbyFracOwn%Pt(i)
            JacBigA%Num(start + j) = FracbyFracOwn%Num(j + FracbyFracOwn%Pt(i)) + nbNodeLocal ! col
         end do
         start = start + FracbyFracOwn%Pt(i + 1) - FracbyFracOwn%Pt(i)

         ! a23(i,:)
         do j = 1, CellbyFracOwn%Pt(i + 1) - CellbyFracOwn%Pt(i)
            JacBigA%Num(start + j) = CellbyFracOwn%Num(j + CellbyFracOwn%Pt(i)) + nbNodeLocal + nbFracLocal ! col
         end do
         start = start + CellbyFracOwn%Pt(i + 1) - CellbyFracOwn%Pt(i)

         ! a24 = a25 = 0
      end do

      do i = 1, nbCellLocal

         ! a31(i,:)
         do j = 1, NodebyCellLocal%Pt(i + 1) - NodebyCellLocal%Pt(i)
            JacBigA%Num(start + j) = NodebyCellLocal%Num(j + NodebyCellLocal%Pt(i))
         end do
         start = start + NodebyCellLocal%Pt(i + 1) - NodebyCellLocal%Pt(i)

         ! a32(i,:)
         do j = 1, FracbyCellLocal%Pt(i + 1) - FracbyCellLocal%Pt(i)
            JacBigA%Num(start + j) = FracbyCellLocal%Num(j + FracbyCellLocal%Pt(i)) + nbNodeLocal ! col
         end do
         start = start + FracbyCellLocal%Pt(i + 1) - FracbyCellLocal%Pt(i)

         ! a33(i,:)
         JacBigA%Num(start + 1) = i + nbNodeLocal + nbFracLocal ! col
         start = start + 1

         ! a34 = a35 = 0
      end do

      do i = 1, nbWellInjOwn

         ! a41(i,:)
         do j = 1, NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i)
            JacBigA%Num(start + j) = NodebyWellInjLocal%Num(j + NodebyWellInjLocal%Pt(i))
         end do
         start = start + NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i)

         ! a44(i,:)
         JacBigA%Num(start + 1) = i + nbNodeLocal + nbFracLocal + nbCellLocal
         start = start + 1
      end do

      do i = 1, nbWellProdOwn

         ! a51(i,:)
         do j = 1, NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i)
            JacBigA%Num(start + j) = NodebyWellProdLocal%Num(j + NodebyWellProdLocal%Pt(i))
         end do
         start = start + NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i)

         ! a55(i,:)
         JacBigA%Num(start + 1) = i + nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal
         start = start + 1
      end do

      ! a61 in JacBigA, pairing mswell_node with reservoir_node
      mswells_JacBigA_num_off = start
      do k = 1, NbMSWellLocal
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) ! index of s in the reservoir (local mesh)

            if (IdNodeLocal(num_s)%Proc /= "o") cycle ! node not own

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx @ MSWells_JacA
            JacBigA%Num(start + 1) = num_s
            start = start + 1 + MSWells_JacA%Pt(unk_idx_s + 1) - MSWells_JacA%Pt(unk_idx_s) !cols of a66 to be filled
         end do
      end do

      start = mswells_JacBigA_num_off + 1 !sum the offset of col a61
      do i = 1, NbMSWellNodeOwn! this is equal to MSWells_JacA%Nb
         ! a66 in JacBigA, copy cols from mswells_JacA using JacBigA offset
         do j = 1, MSWells_JacA%Pt(i + 1) - MSWells_JacA%Pt(i) !cols at row i
            JacBigA%Num(start + 1) = MSWells_JacA%Num(j + MSWells_JacA%Pt(i)) + mswells_JacBigA_unk_lcl_off
            start = start + 1
         end do
         start = start + 1 !sum the offset of col  a61
      end do

      ! Sort JacbigA%Num until mswellnodes
      do i = 1, mswells_JacBigA_row_off !JacBigA%Nb

         do k = 1, JacBigA%Pt(i + 1) - JacBigA%Pt(i)
            do j = JacBigA%Pt(i) + 1, JacBigA%Pt(i + 1) - k

               if (JacBigA%Num(j) > JacBigA%Num(j + 1)) then
                  tmp = JacBigA%Num(j)
                  JacBigA%Num(j) = JacBigA%Num(j + 1)
                  JacBigA%Num(j + 1) = tmp
               end if

            end do
         end do
      end do

   end subroutine Jacobian_StrucJacBigA

   ! We speak in terms in blocks - that means one non-zero means one block that is non zero
   !         node, frac, wellinj, wellprod, mswell_nodes
   ! JacA = | A11, A12, A13, A14 A15 | node own
   !        | A21, A22, 0,   0   0   | frac own
   !        | A31, 0,   A33, 0   0   | wellinj own
   !        | A41, 0,   0,   A44 0   | wellprod own
   !        | A51, 0,   0,   0   A55 | mswell_nodes own
   ! the nonzero structure of Aij is same as aij in JacBigA, i,j=1,2

   ! The non zero structure of aij is based on the connectivities
   ! A11: NodebyNodeOwn if not dir
   ! A12: FracbyNodeOwn if not dir
   ! Rq: if node own is Dir then only one non zero per row for A1*

   ! A21: NodebyFracOwn
   ! A22: FracbyFracOwn
   ! A13: WellInjbyNodeOwn
   ! A14: WellProdbyNodeOwn

   ! A31: NodebyWellInjLocal
   ! A41: NodebyWellProdLocal
   ! A33, A44: diag

   ! Four steps in this subroutine
   !   Number of non zeros each line: NbNnzbyline
   !   bigA%Nb, bigA%Pt using NbNnzbyline
   !   bigA%Num, non zero structure of bigA
   !   arrange bigA%Num such that in each row, the cols is in inscreasing order
   subroutine Jacobian_StrucJacA_fill_Pt(Pt)
      integer, dimension(:), intent(out) :: Pt

      integer, dimension(:), allocatable :: &
         nbNnzbyline ! number of non zeros each line
      integer :: i, nb_rows, nnz, start
      integer :: &
         nbNodeOwn, nbFracOwn, nbWellInjOwn, nbWellProdOwn

      nbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)
      nbFracOwn = NbFracOwn_Ncpus(commRank + 1)
      nbWellInjOwn = NbWellInjOwn_Ncpus(commRank + 1)
      nbWellProdOwn = NbWellProdOwn_Ncpus(commRank + 1)

      nb_rows = nbNodeOwn + nbFracOwn + nbWellInjOwn + nbWellProdOwn + nbMSWellNodeOwn
      allocate (nbNnzbyLine(nb_rows))

#ifndef NDEBUG
      if (size(Pt) /= nb_rows + 1) &
         call CommonMPI_abort("Inconsitent Pt size in Jacobian_StrucJacA_fill_Pt")
#endif

      ! A1* in JacA
      do i = 1, nbNodeOwn

         nbNnzbyLine(i) = NodebyNodeOwn%Pt(i + 1) - NodebyNodeOwn%Pt(i) &
                          + FracbyNodeOwn%Pt(i + 1) - FracbyNodeOwn%Pt(i) &
                          + WellInjbyNodeOwn%Pt(i + 1) - WellInjbyNodeOwn%Pt(i) &
                          + WellProdbyNodeOwn%Pt(i + 1) - WellProdbyNodeOwn%Pt(i) &
                          + MSWellbyNodeOwn%Pt(i + 1) - MSWellbyNodeOwn%Pt(i)
      end do

      ! A2* in JacA, rq: frac is not dir face
      start = nbNodeOwn
      do i = 1, nbFracOwn
         nbNnzbyLine(start + i) = NodebyFracOwn%Pt(i + 1) - NodebyFracOwn%Pt(i) &
                                  + FracbyFracOwn%Pt(i + 1) - FracbyFracOwn%Pt(i)
      end do

      ! A3* in JacA
      start = start + nbFracOwn
      do i = 1, nbWellInjOwn
         nbNnzbyLine(start + i) = &
            NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i) + 1 ! A44 is diagonal
      end do

      ! A4* in JacA
      start = start + nbWellInjOwn
      do i = 1, nbWellProdOwn
         nbNnzbyLine(start + i) = &
            NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i) + 1 ! A55 is diagonal
      end do

      ! A5* in JacA
      start = start + nbWellProdOwn
      mswells_JacA_row_off = start
      do i = 1, NbMSWellNodeOwn! this is equal to MSWells_JacA%Nb
         nbNnzbyLine(start + i) = &
            MSWells_JacA%Pt(i + 1) - MSWells_JacA%Pt(i) + 1 ! copy the nbcols of MSWells_Jac
         ! and sum the contribution of the reservoir paired node
      end do

      ! JacA%Pt
      nnz = 0
      Pt(1) = 0
      do i = 1, nb_rows
         nnz = nnz + nbNnzbyLine(i)
         Pt(i + 1) = nnz
      end do

      deallocate (nbNnzbyLine)

   end subroutine Jacobian_StrucJacA_fill_Pt

   ! FIXME: This a convenience temporary function that is to be removed later
   !> Fills the column number of a sparse matrix (SM)
   !> with dof stored in an adjacency CSRInt whose main attributes are Pt and Num (cf. type(CSR))
   !> j is the rank to the node of the adjacency graph currently considered
   !> coff is the column offset due to the face that there is several adjacency graph/tables to be considered
   !> col is the array storing the columns number in a contiguous fashion
   !> roff is the row offset in cols to find the first columns corresponding to the adjacency relations
   !>    it is incremented so that ad the end of the subroutine it can be reused (the same way as an output iterator)
   pure subroutine fill_SM_column_number(adj, j, coff, col, roff)
      type(CSR), intent(in) :: adj
      integer, intent(in) :: j
      integer, intent(in) :: coff
      integer, dimension(:), intent(inout) :: col
      integer, intent(inout) :: roff

      integer :: n

      n = adj%Pt(j + 1) - adj%Pt(j)
      col(roff + 1:roff + n) = adj%Num(adj%Pt(j) + 1:adj%Pt(j + 1)) + coff
      roff = roff + n

   end subroutine fill_SM_column_number

   subroutine Jacobian_StrucJacA_fill_columns(M)
      type(CSRArray2dble), intent(inout) :: M

      integer :: i, j, start, k, s, num_s, unk_idx_s, glb_unk_col_idx_s
      integer :: nbNodeOwn, nbFracOwn, nbWellInjOwn, nbWellProdOwn
      integer :: column_offset(4)

      nbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)
      nbFracOwn = NbFracOwn_Ncpus(commRank + 1)
      nbWellInjOwn = NbWellInjOwn_Ncpus(commRank + 1)
      nbWellProdOwn = NbWellProdOwn_Ncpus(commRank + 1)

      column_offset(1) = 0
      column_offset(2) = column_offset(1) + NbNodeLocal_Ncpus(commRank + 1)
      column_offset(3) = column_offset(2) + NbFracLocal_Ncpus(commRank + 1)
      column_offset(4) = column_offset(3) + NbWellInjLocal_Ncpus(commRank + 1)

      M%Num = 0
      start = 0

#ifndef NDEBUG
      if (.not. associated(M%Pt)) &
         call CommonMPI_abort("JacA%Pt should be a valid pointer.")
      if (M%Pt(1) /= 0) &
         call CommonMPI_abort("JacA%Pt(0) should be 0.")
#endif

      ! First set of (block) rows conservation law at reversoir nodes
      do i = 1, nbNodeOwn
         call fill_SM_column_number(NodebyNodeOwn, i, column_offset(1), M%Num, start)
         call fill_SM_column_number(FracbyNodeOwn, i, column_offset(2), M%Num, start)
         call fill_SM_column_number(WellInjbyNodeOwn, i, column_offset(3), M%Num, start)
         call fill_SM_column_number(WellProdbyNodeOwn, i, column_offset(4), M%Num, start)
         ! A15(i;:)
         do j = 1, MSWellbyNodeOwn%Pt(i + 1) - MSWellbyNodeOwn%Pt(i)
            k = MSWellbyNodeOwn%Num(j + MSWellbyNodeOwn%Pt(i)) !number of the mswell which is linked with the reservoir node

            !Find the reservoir node
            do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

               if (IncIdxNodebyMSWellLocal%Val(s) == i) then

                  unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)           !mswell_node unkown idx @ MSWells_JacA
                  glb_unk_col_idx_s = unk_idx_s + mswells_JacA_unk_lcl_off  !unkown col idx @ JacA
                  JacA%Num(start + j) = glb_unk_col_idx_s
                  start = start + 1
                  exit

               end if
            end do
         end do
#ifndef NDEBUG
         if (M%Pt(i + 1) /= start) &
            call CommonMPI_abort("Wrong offset in filling JacA node rows.")
#endif
      end do

      ! Second set of (block) rows conservation law at fracture faces
      ! wells are not connected to fracture faces
      do i = 1, nbFracOwn
         call fill_SM_column_number(NodebyFracOwn, i, column_offset(1), M%Num, start)
         call fill_SM_column_number(FracbyFracOwn, i, column_offset(2), M%Num, start)
#ifndef NDEBUG
         if (M%Pt(NbNodeOwn + i + 1) /= start) &
            call CommonMPI_abort("Wrong offset in filling JacA fracture rows.")
#endif
      end do

      ! Third set of (block) rows conservation law at injection wells
      do i = 1, nbWellInjOwn
         call fill_SM_column_number(NodebyWellInjLocal, i, column_offset(1), M%Num, start)
         M%Num(start + 1) = i + column_offset(3)
         start = start + 1
#ifndef NDEBUG
         if (M%Pt(nbNodeOwn + nbFracOwn + i + 1) /= start) &
            call CommonMPI_abort("Wrong offset in filling JacA injection wells rows.")
#endif
      end do

      ! Fourth set of (block) rows conservation law at production wells
      do i = 1, nbWellProdOwn
         call fill_SM_column_number(NodebyWellProdLocal, i, column_offset(1), M%Num, start)
         M%Num(start + 1) = i + column_offset(4)
         start = start + 1
#ifndef NDEBUG
         if (M%Pt(nbNodeOwn + nbFracOwn + nbWellInjOwn + i + 1) /= start) &
            call CommonMPI_abort("Wrong offset in filling JacA production wells rows.")
#endif
      end do

      ! A51 in JacBigA, pairing mswell_node with reservoir_node
      mswells_JacA_num_off = start
      do k = 1, NbMSWellLocal
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) ! index of s in the reservoir (local mesh)

            if (IdNodeLocal(num_s)%Proc /= "o") cycle ! node not own

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx @ MSWells_JacA
            JacA%Num(start + 1) = num_s
            start = start + 1 + MSWells_JacA%Pt(unk_idx_s + 1) - MSWells_JacA%Pt(unk_idx_s) !cols of A55 to be filled
         end do
      end do

      ! A55 in JacBigA, copy cols from mswells_JacA using JacA offset
      start = mswells_JacA_num_off + 1 !sum the offset of col  A51
      do i = 1, NbMSWellNodeOwn! this is equal to MSWells_JacA%Nb
         do j = 1, MSWells_JacA%Pt(i + 1) - MSWells_JacA%Pt(i) !cols at row i
            JacA%Num(start + 1) = MSWells_JacA%Num(j + MSWells_JacA%Pt(i)) + mswells_JacA_unk_lcl_off
            start = start + 1
         end do
         start = start + 1 !sum the offset of col  A51
      end do

   end subroutine Jacobian_StrucJacA_fill_columns

   ! Non zero sturcture of JacA: JacBigA after Schur
   subroutine Jacobian_StrucJacA
      ! We speak in terms in blocks - that means one non-zero means one block that is non zero
      !         node, frac, wellinj, wellprod
      ! JacA = | A11, A12, A13, A14 | node own
      !        | A21, A22, 0,   0   | frac own
      !        | A31, 0,   A33, 0   | wellinj own
      !        | A41, 0,   0,   A44 | wellprod own
      ! the nonzero structure of Aij is same as aij in JacBigA, i,j=1,2

      ! The non zero structure of aij is based on the connectivities
      ! A11: NodebyNodeOwn if not dir
      ! A12: FracbyNodeOwn if not dir
      ! Rq: if node own is Dir then only one non zero per row for A1*

      ! A21: NodebyFracOwn
      ! A22: FracbyFracOwn
      ! A13: WellInjbyNodeOwn
      ! A14: WellProdbyNodeOwn

      ! A31: NodebyWellInjLocal
      ! A41: NodebyWellProdLocal
      ! A33, A44: diag

      ! Four steps in this subroutine
      !   Number of non zeros each line: NbNnzbyline
      !   bigA%Nb, bigA%Pt using NbNnzbyline
      !   bigA%Num, non zero structure of bigA
      !   arrange bigA%Num such that in each row, the cols is in inscreasing order

      integer :: i, nb_rows, nb_rows_no_mswells, nb_cols, nnz

      nb_rows_no_mswells = NbNodeOwn_Ncpus(commRank + 1) &
                           + NbFracOwn_Ncpus(commRank + 1) &
                           + NbWellInjOwn_Ncpus(commRank + 1) &
                           + NbWellProdOwn_Ncpus(commRank + 1)

      nb_rows = nb_rows_no_mswells + NbMSWellNodeOwn

      nb_cols = NbNodeLocal_Ncpus(commRank + 1) &
                + NbFracLocal_Ncpus(commRank + 1) &
                + NbCellLocal_Ncpus(commRank + 1) &
                + NbWellInjLocal_Ncpus(commRank + 1) &
                + NbWellProdLocal_Ncpus(commRank + 1) &
                + NbMSWellNodeLocal

      mswells_JacA_unk_lcl_off = NbNodeLocal_Ncpus(commRank + 1) &
                                 + NbFracLocal_Ncpus(commRank + 1) &
                                 + NbWellInjLocal_Ncpus(commRank + 1) &
                                 + NbWellProdLocal_Ncpus(commRank + 1)

      JacA%Nb = nb_rows
      allocate (JacA%Pt(nb_rows + 1))

      call Jacobian_StrucJacA_fill_Pt(JacA%Pt)

#ifndef NDEBUG
      if (nb_rows_no_mswells /= mswells_JacA_row_off) &
         call CommonMPI_abort("Wrong offset in filling JacA mswells rows.")
#endif
      nnz = JacA%Pt(nb_rows + 1)

      allocate (JacA%Num(nnz))

      call Jacobian_StrucJacA_fill_columns(JacA)

      ! Sort JacA%Num so that column indices are in ascending order
      !               but rows of  mswells_nodes
      do i = 1, nb_rows_no_mswells
         call quicksort(JacA%Num, JacA%Pt(i) + 1, JacA%Pt(i + 1))
      end do

      allocate (JacA%Val(NbCompThermique, NbCompThermique, nnz)) ! number of non zero
      allocate (Sm(NbCompThermique, nb_rows))
      allocate (csrK(nb_cols))
      csrK = 0
      allocate (csrSR(nb_cols))
      csrSR = 0

   end subroutine Jacobian_StrucJacA

   ! TODO: to be moved elsewhere - is there already a quicksort routine for wells?
   ! Adapted from https://gist.github.com/t-nissie/479f0f16966925fa29ea
   ! Author: t-nissie
   ! License: GPLv3
   pure recursive subroutine quicksort(a, first, last)
      integer, dimension(:), intent(inout) ::  a
      integer, intent(in) :: first, last

      integer i, j, n, tmp

      n = a((first + last)/2)
      i = first
      j = last
      do
         do while (a(i) < n)
            i = i + 1
         end do
         do while (n < a(j))
            j = j - 1
         end do
         if (i >= j) exit
         tmp = a(i); a(i) = a(j); a(j) = tmp
         i = i + 1
         j = j - 1
      end do
      if (first < i - 1) call quicksort(a, first, i - 1)
      if (j + 1 < last) call quicksort(a, j + 1, last)
   end subroutine quicksort

   subroutine Jacobian_free

      deallocate (JacBigA%Pt)
      deallocate (JacBigA%Num)
      deallocate (JacBigA%Val)

      deallocate (JacA%Pt)
      deallocate (JacA%Num)
      deallocate (JacA%Val)

      deallocate (bigSm)
      deallocate (Sm)

      deallocate (csrK)
      deallocate (csrSR)

   end subroutine Jacobian_free

   ! num transformation, used for TkLocal(k)%(i,j), where i, j is node or frac
   ! from i, j to its num (local)
   ! k: loop of cell
   subroutine Jacobian_RowCol_KSR(k, &
                                  nbNodeCell, nbFracCell, &
                                  rowk, colk, &
                                  rowSR, colSR)

      integer, intent(in) :: k, &
                             nbNodeCell, nbFracCell

      integer, intent(out) :: &
         rowk, colk, &
         rowSR(NbNodeCellMax + NbFracCellMax), &
         colSR(NbNodeCellMax + NbFracCellMax)

      integer :: i, in

      ! rowk
      rowk = k &
             + NbNodeOwn_Ncpus(commRank + 1) &
             + NbFracOwn_Ncpus(commRank + 1)

      ! colk
      colk = k &
             + NbNodeLocal_Ncpus(commRank + 1) &
             + NbFracLocal_Ncpus(commRank + 1)

      ! colSR, nodes
      do i = 1, nbNodeCell
         colSR(i) = NodebyCellLocal%Num(i + NodebyCellLocal%Pt(k))
      end do

      ! colSR, frac
      do i = 1, nbFracCell
         colSR(nbNodeCell + i) = FracbyCellLocal%Num(i + FracbyCellLocal%Pt(k)) &
                                 + NbNodeLocal_Ncpus(commRank + 1)
      end do

      rowSR(:) = 0

      ! rowSR, nodes
      do i = 1, nbNodeCell

         ! node
         ! Rq. if i is dir Darcy, rowSR(i) will not be used for Darcy
         !     if i is dir Fourier, rowSR(i) will not be used for Fourier
         rowSR(i) = NodebyCellLocal%Num(i + NodebyCellLocal%Pt(k))
      end do

      ! rowSR, frac
      do i = 1, nbFracCell

         in = FracbyCellLocal%Num(i + FracbyCellLocal%Pt(k))
         if (FracToFaceLocal(in) <= NbFaceOwn_Ncpus(commRank + 1)) then ! own, rq: in is frac number
            rowSR(nbNodeCell + i) = in + NbNodeOwn_Ncpus(commRank + 1)
         end if
      end do

   end subroutine Jacobian_RowCol_KSR

   ! num transformation, used for TkFracLocal(k)%(i,j), where i, j is node
   ! from i, j to its num (local)
   ! k: loop of frac
   subroutine Jacobian_RowCol_FR(k, &
                                 nbNodeFrac, &
                                 rowk, colk, &
                                 rowSR, colSR)

      integer, intent(in) :: k, nbNodeFrac

      integer, intent(out) :: &
         rowk, colk, &
         rowSR(NbNodeFaceMax), &
         colSR(NbNodeFaceMax)

      integer :: i, fk

      fk = FracToFaceLocal(k) ! this frac is which face

      ! rowk
      rowk = k + NbNodeOwn_Ncpus(commRank + 1)

      ! colk
      colk = k + NbNodeLocal_Ncpus(commRank + 1)

      ! rowR
      rowSR(:) = 0
      do i = 1, nbNodeFrac
         rowSR(i) = NodebyFaceLocal%Num(i + NodebyFaceLocal%Pt(fk))
      end do

      ! colSR = rowR
      colSR(:) = rowSR(:)

   end subroutine Jacobian_RowCol_FR

   subroutine Jacobian_GetSolCell_C(increments_pointers) &
      bind(C, name="Jacobian_GetSolCell")

      type(Newton_increments_pointers), intent(in), value :: increments_pointers
      type(Newton_increments) :: increments

      call Newton_pointers_to_values(increments_pointers, increments)
      call Jacobian_GetSolCell( &
         increments%nodes, increments%fractures, increments%cells &
         )

   end subroutine Jacobian_GetSolCell_C

   ! A3 * IncrementNodeFrac + A4 * IncrementCell = bigSm (cell part)
   ! A4 is a block diag matrix,
   ! its block element has been done LU factorization
   subroutine Jacobian_GetSolCell( &
      NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell &
      )

      real(c_double), dimension(:, :), intent(in) :: &
         NewtonIncreNode, &
         NewtonIncreFrac

      real(c_double), dimension(:, :), intent(out) :: &
         NewtonIncreCell

      integer :: k, rowk, s, cols, nums, i, j, nz
      integer :: info, errcode, Ierr
      integer :: NbNodeFracOwn, NbNodeLocal
      double precision :: tmpInc(NbCompThermique)

      NbNodeFracOwn = NbNodeOwn_Ncpus(commRank + 1) + NbFracOwn_Ncpus(commRank + 1)
      NbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)

      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         rowk = k + NbNodeFracOwn

         tmpInc(:) = 0.d0

         do s = JacbigA%Pt(rowk) + 1, JacBigA%Pt(rowk + 1) - 1  ! A3
            cols = JacbigA%Num(s)

            nz = JacBigA%Pt(rowk + 1) ! A4

            ! cols is in part node
            if (cols <= NbNodeLocal) then

               do i = 1, NbCompThermique
                  do j = 1, NbCompThermique

                     tmpInc(i) = tmpInc(i) &
                                 + JacBigA%Val(j, i, s)*NewtonIncreNode(j, cols)
                  end do
               end do

            else
               ! cols is in part frac
               ! this col corresponts to frac nums
               nums = cols - NbNodeLocal

               do i = 1, NbCompThermique
                  do j = 1, NbCompThermique

                     tmpInc(i) = tmpInc(i) &
                                 + JacBigA%Val(j, i, s)*NewtonIncreFrac(j, nums)
                  end do
               end do
            end if
         end do

         NewtonIncreCell(1:NbCompThermique, k) = &
            bigSm(1:NbCompThermique, rowk) - tmpInc(1:NbCompThermique)

         call dgetrs('T', NbCompThermique, 1, &
                     JacBigA%Val(:, :, nz), NbCompThermique, pivot(:, k), &
                     NewtonIncreCell(1:NbCompThermique, k), NbCompThermique, info)

         if (info /= 0) then
            write (0, '(A,I5)', advance='no') "dgetrs error in inverse Schur complement k =", k
            write (0, '(A,I5)') ",  info =", info
            call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
         end if

      end do

   end subroutine Jacobian_GetSolCell

end module Jacobian
