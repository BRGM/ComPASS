!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecdMSWells

#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use mpi
   use CommonMPI
   use IncCVReservoir
   use DefModel
   use NumbyContext
   use MeshSchema
   use Thermodynamics
   use IncPrimSecdTypes
   use IncPrimSecd
   use IncCVMSWells
   use IncPrimSecd
   use MeshSchemaMSWells
   use Newton
#else
   use iso_c_binding, only: c_double
   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir

   use DefModel, only: &
      NbComp, NbPhase, MCP, &
      NbIncTotalPrimMax, NbEqFermetureMax, NbIncTotalMax, &
      pschoice, psprim, pssecd, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IndThermique, NbCompThermique

   use NumbyContext, only: &
      NbEqFermeture_ctx, NbIncTotalPrim_ctx, &
      NbIncPTC_ctx, NbIncTotal_ctx, NbCompCtilde_ctx

   use MeshSchema, only: &
      NbNodeLocal_Ncpus, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus

   use Thermodynamics, only: f_Fugacity
   use IncPrimSecdTypes, only: ControlVolumeInfo, IncPrimSecdTypes_collect_cv_info

   use IncPrimSecd, only: &
      IncPrimSecdTypes_collect_cv_info, &
      IncPrimSecd_dFsurdX_cv, &
      IncPrimSecd_ps_cv, &
      IncPrimSecd_dXssurdXp_cv

   use IncCVMSWells, only: &
      TYPE_IncCVReservoir, &
      IncMSWell

   use MeshSchemaMSWells, only: &
      NbMSWellNodeOwn, &
      NbMSWellNodeLocal, &
      IncIdxNodebyMSWellLocal

   use Newton, only: &
      Newton_increments_pointers, Newton_increments, Newton_pointers_to_values
#endif
   implicit none

   type :: TYPE_IncPrimSecd
      double precision  :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax), &
         SmF(NbEqFermetureMax)

      !num inc prim secd
      integer :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)
   end type TYPE_IncPrimSecd

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   TYPE(TYPE_IncPrimSecd), allocatable, dimension(:), target, protected  :: &
      IncPrimSecdMSWell !<  Injector and Producer MSWells-IncPrimSecd

   public :: &
      IncPrimSecdMSWells_allocate, &
      IncPrimSecdMSWells_free, &
      IncPrimSecdMSWells_compute

   private :: &
      IncPrimSecdMSWells_compute_cv

contains

   subroutine IncPrimSecdMSWells_allocate()

      integer :: Nb, Nnz

      Nb = NodebyMSWellLocal%Nb
      Nnz = NodebyMSWellLocal%Pt(Nb + 1)
      allocate (IncPrimSecdMSWell(Nnz))

   end subroutine IncPrimSecdMSWells_allocate

   subroutine IncPrimSecdMSWells_free()

      deallocate (IncPrimSecdMSWell)

   end subroutine IncPrimSecdMSWells_free

   !> \brief Main subroutine of this module,
   !! compute dXssurdXp, SmdXs for each mswell
   subroutine IncPrimSecdMSWells_compute() &
      bind(C, name="IncPrimSecdMSWells_compute")

      integer :: s, k, nbwells

      !For all mswell
      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do k = 1, nbwells

         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
            call IncPrimSecdMSWells_compute_cv( &
               IncMSWell(s)%coats, &
               IncPrimSecdMSWell(s))

         end do
      end do

   end subroutine IncPrimSecdMSWells_compute

   subroutine IncPrimSecdMSWells_compute_cv(inc, inc_primsecd)
      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc
      type(TYPE_IncPrimSecd), intent(out)::inc_primsecd

      ! tmp
      double precision :: &
         dFsurdX(NbIncTotalMax, NbEqFermetureMax) ! (col,row) index order
      type(ControlVolumeInfo) :: cv_info
      real(c_double)          :: dpadS(NbPhase)

      dpadS = 0.d0 ! derivative of phase/capillar pressure with respect saturations
      ! init tmp values for each cv
      call IncPrimSecdTypes_collect_cv_info(inc%ic, cv_info)

      !< compute dF/dX
      !< dFsurdX: (col, row) index order
      inc_primsecd%SmF(:) = 0.d0
      call IncPrimSecd_dFsurdX_cv(cv_info, inc, dpadS, dFsurdX, inc_primsecd%SmF)

      !< choose inconnues prim and secd
      call IncPrimSecd_ps_cv(cv_info, inc, dFsurdX, pschoice, &
                             inc_primsecd%NumIncTotalPrimCV, inc_primsecd%NumIncTotalSecondCV)

      !< compute dXssurdxp and SmdXs
      inc_primsecd%dXssurdXp(:, :) = 0.d0
      inc_primsecd%SmdXs(:) = 0.d0
      call IncPrimSecd_dXssurdXp_cv(cv_info, inc, dFsurdX, inc_primsecd%SmF, &
                                    inc_primsecd%NumIncTotalPrimCV, inc_primsecd%NumIncTotalSecondCV, &
                                    inc_primsecd%dXssurdXp, inc_primsecd%SmdXs)

   end subroutine IncPrimSecdMSWells_compute_cv

   !> \brief  Link with the C code
   subroutine IncPrimSecdMSWells_PrimToSecd_C(increments_pointers) &
      bind(C, name="IncPrimSecdMSWells_PrimToSecd")

      type(Newton_increments_pointers), intent(in), value :: increments_pointers
      type(Newton_increments) :: increments

      call Newton_pointers_to_values(increments_pointers, increments)
      call IncPrimSecdMSWells_PrimToSecd(increments%mswell_nodes)

   end subroutine IncPrimSecdMSWells_PrimToSecd_C

   !> \brief Compute secd values using prim values
   !! secd = SmdX - dXssurdXp * prim
   !! for node and frac and cell
   ! v for variation ?
   subroutine IncPrimSecdMSWells_PrimToSecd(var_inc)

      real(c_double), dimension(:, :), intent(inout) :: &
         var_inc

      integer :: s, l, nbwells, unk_idx_s
      integer ::  ic, i, iph
      integer :: &
         NbPhasePresente, NbEqFermeture, NbNodeLocal, &
         NbIncPTC, NbIncTotal, NbIncPTCPrim, NbIncTotalPrim

      double precision :: &
         xp(NbCompThermique), &
         xs(NbEqFermetureMax)

      !For all mswell
      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do l = 1, nbwells
         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(l) + 1, NodebyMSWellLocal%Pt(l + 1)

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown mswell-idx. To be used to access var_inc

            ic = IncMSWell(s)%coats%ic
            NbPhasePresente = NbPhasePresente_ctx(ic)
            NbEqFermeture = NbEqFermeture_ctx(ic)
            NbIncPTC = NbIncPTC_ctx(ic)
            NbIncTotal = NbIncTotal_ctx(ic)
            NbIncTotalPrim = NbIncTotalPrim_ctx(ic)
            NbIncPTCPrim = NbIncPTC - NbEqFermeture

            xp(1:NbCompThermique) = var_inc(1:NbCompThermique, unk_idx_s)
            xs(1:NbEqFermeture) = IncPrimSecdMSWell(s)%SmdXs(1:NbEqFermeture)

            ! http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
            ! y := alpha*A*x + beta*y
            ! Compute xs = xs - dXssurdXp * xp, where xs=SmdXs and xp=var_prim
            call dgemv('T', NbIncTotalPrim, NbEqFermeture, &
                       -1.d0, IncPrimSecdMSWell(s)%dXssurdXp(:, :), NbIncTotalPrimMax, &
                       xp(:), 1, -1.d0, xs(:), 1)

            !-----------------------------------------------------
            ! update var_inc with the primary/secondary increments
            !-----------------------------------------------------

            var_inc(:, unk_idx_s) = 0.d0

            ! copy prim P,T,C,S
            do i = 1, NbIncTotalPrim
               var_inc(IncPrimSecdMSWell(s)%NumIncTotalPrimCV(i), unk_idx_s) = xp(i)
            end do

            ! copy secd P,T,C
            do i = 1, NbEqFermeture
               var_inc(IncPrimSecdMSWell(s)%NumIncTotalSecondCV(i), unk_idx_s) = xs(i)
            end do

            ! fill secd S
            ! if NbPhasePresente=1,
            !    then this phase is eliminated
            ! else last saturation is eliminated, the others are prim (in reservoir dof)
            !    eliminated S = - sum_{S^alpha is prim} S^alpha
            ! iph is last present phase in vector (P,T,C,S,n_i)
            iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente, ic)
            var_inc(iph, unk_idx_s) = 0.d0 !> \todo FIXME: is useful, because wrong numerotation of FreeFlow sat
            do i = 1, NbPhasePresente - 1
               var_inc(iph, unk_idx_s) = var_inc(iph, unk_idx_s) - xp(NbIncPTCPrim + i)
            end do

            ! term prim n_k(X_j^n), components which are present only in absent phase(s)
            ! n_k(X_j^n) are not part of NbIncTotal
            ! copy prim n_i
            do i = 1, NbCompCtilde_ctx(ic) ! =NbCompThermique-NbIncTotalPrim
               var_inc(NbIncTotal + i, unk_idx_s) = xp(NbIncTotalPrim + i)
            end do

         end do
      end do

   end subroutine IncPrimSecdMSWells_PrimToSecd

end module IncPrimSecdMSWells
