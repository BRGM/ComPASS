!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Newton

   use, intrinsic :: iso_c_binding, only: &
      c_double, c_ptr, c_f_pointer

   use mpi ! FIXME: when using 'use only' MPI_Allreduce is not found on some platform

   use CommonMPI, only: &
      commRank, ComPASS_COMM_WORLD, CommonMPI_abort
   use DefModel, only: &
      NbIncTotalMax, NbPhase, NbComp, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IndThermique, MCP
#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      DIPHASIC_CONTEXT, GAS_PHASE, LIQUID_PHASE
#endif

   use MeshSchema, only: &
      IdNodeLocal, &
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbCellLocal_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus, &
      NbCellOwn_Ncpus, NbFracOwn_Ncpus, NbNodeOwn_Ncpus, &
      NbMSWellNodeLocal_Ncpus, NbMSWellNodeOwn_Ncpus, &
      NodebyMSWellLocal, NbMSWellLocal_Ncpus

   use NumbyContext, only: &
      NbIncPTC_ctx, NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   use SchemeParameters, only: &
      NewtonIncreObj_C, NewtonIncreObj_P, NewtonIncreObj_S, NewtonIncreObj_T, eps

   use IncCVReservoir, only: &
      IncNode, IncCell, IncFrac

   use MeshSchemaMSWells, only: &
      NbMSWellNodeOwn, &
      NbMSWellNodeLocal, &
      IncIdxNodebyMSWellLocal

   use IncCVMSWells, only: &
      IncMSWell

   implicit none

   type, bind(C) :: Newton_increments_pointers
      type(c_ptr) :: nodes
      type(c_ptr) :: fractures
      type(c_ptr) :: cells
      type(c_ptr) :: injectors
      type(c_ptr) :: producers
      type(c_ptr) :: mswell_nodes
   end type Newton_increments_pointers

   type Newton_increments
      real(c_double), pointer :: nodes(:, :)
      real(c_double), pointer :: fractures(:, :)
      real(c_double), pointer :: cells(:, :)
      real(c_double), pointer :: injectors(:)
      real(c_double), pointer :: producers(:)
      real(c_double), pointer :: mswell_nodes(:, :)
   end type Newton_increments

   real(c_double) last_max_inc_mswells

   public :: &
      Newton_pointers_to_values

contains

   subroutine Newton_pointers_to_values(increment_pointers, increment_values)

      type(Newton_increments_pointers), intent(in), value :: increment_pointers
      type(Newton_increments), intent(out) :: increment_values

      call c_f_pointer( &
         increment_pointers%nodes, increment_values%nodes, &
         shape=[NbIncTotalMax, NbNodeLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%fractures, increment_values%fractures, &
         shape=[NbIncTotalMax, NbFracLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%cells, increment_values%cells, &
         shape=[NbIncTotalMax, NbCellLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%injectors, increment_values%injectors, &
         shape=[NbWellInjLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%producers, increment_values%producers, &
         shape=[NbWellProdLocal_Ncpus(commRank + 1)] &
         )

      call c_f_pointer( &
         increment_pointers%mswell_nodes, increment_values%mswell_nodes, &
         shape=[NbIncTotalMax, NbMSWellNodeLocal_Ncpus(commRank + 1)] &
         )

   end subroutine Newton_pointers_to_values

   !This function does not compute any relaxation for MSWellNodes
   function Newton_compute_relaxation_C(increments_pointers) &
      result(relaxation) &
      bind(C, name="Newton_compute_relaxation")

      type(Newton_increments_pointers), intent(in), value :: increments_pointers
      real(c_double) :: relaxation
      type(Newton_increments) :: increments

      call Newton_pointers_to_values(increments_pointers, increments)
      call Newton_compute_relaxation(increments, relaxation)

   end function Newton_compute_relaxation_C

   !> \brief Compute relaxation in Newton.
   !!
   !! relax = min(1, IncreObj/NewtonIncreObjMax) <br>
   !! where IncreObj is set by the user in SchemeParameters.F90 <br>
   !! and NewtonIncreObjMax is the maximum of the Nemton increment
   !! in current iteration
   subroutine Newton_compute_relaxation(increments, relax)
      type(Newton_increments), intent(in) :: increments
      real(c_double), intent(out) :: relax

      double precision :: &
         incremaxlocal_P, &
         incremaxlocal_T, &
         incremaxlocal_C(NbComp, NbPhase), &
         incremaxlocal_S(NbPhase), &
         relaxlocal

      integer :: k, i, ic, iph, icp, j, Ierr
      integer :: NbIncPTC

      incremaxlocal_P = 0.d0

#ifdef _THERMIQUE_
      incremaxlocal_T = 0.d0
#endif
      incremaxlocal_C(:, :) = 0.d0
      incremaxlocal_S(:) = 0.d0

      ! FIXME: factorize nodes/fractures/cells (3x the same code !)
      ! max Newton increment node
      do k = 1, NbNodeOwn_Ncpus(commRank + 1)

         if (IdNodeLocal(k)%P /= "d") then

            ic = IncNode(k)%ic
            NbIncPTC = NbIncPTC_ctx(ic)

            incremaxlocal_P = max(incremaxlocal_P, abs(increments%nodes(1, k)))

            do i = 2 + IndThermique, NbIncPTC
               icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
               iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

               incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(increments%nodes(i, k)))
            enddo

            do i = 1, NbPhasePresente_ctx(ic)
               iph = NumPhasePresente_ctx(i, ic)

               incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(increments%nodes(iph + NbIncPTC, k)))
            end do
         end if

#ifdef _THERMIQUE_
         if (IdNodeLocal(k)%T /= "d") then
            incremaxlocal_T = max(incremaxlocal_T, abs(increments%nodes(2, k)))
         end if
#endif

      end do

      ! max Newton increment fracture face
      do k = 1, NbFracOwn_Ncpus(commRank + 1)

         ic = IncFrac(k)%ic
         NbIncPTC = NbIncPTC_ctx(ic)

         incremaxlocal_P = max(incremaxlocal_P, abs(increments%fractures(1, k)))

#ifdef _THERMIQUE_
         incremaxlocal_T = max(incremaxlocal_T, abs(increments%fractures(2, k)))
#endif

         do i = 2 + IndThermique, NbIncPTC
            icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
            iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

            incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(increments%fractures(i, k)))
         enddo

         do i = 1, NbPhasePresente_ctx(ic)
            iph = NumPhasePresente_ctx(i, ic)

            incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(increments%fractures(iph + NbIncPTC, k)))
         end do
      end do

      ! max Newton increment cell
      do k = 1, NbCellOwn_Ncpus(commRank + 1)

         ic = IncCell(k)%ic
         NbIncPTC = NbIncPTC_ctx(ic)

         incremaxlocal_P = max(incremaxlocal_P, abs(increments%cells(1, k)))

#ifdef _THERMIQUE_
         incremaxlocal_T = max(incremaxlocal_T, abs(increments%cells(2, k)))
#endif

         do i = 2 + IndThermique, NbIncPTC
            icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
            iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

            incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(increments%cells(i, k)))
         enddo

         do i = 1, NbPhasePresente_ctx(ic)
            iph = NumPhasePresente_ctx(i, ic)

            incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(increments%cells(iph + NbIncPTC, k)))
         end do
      end do

      ! Consider well pressure head variation for relaxation
      ! incremaxlocal_P = max(incremaxlocal_P, maxval(abs(increments%injectors)))
      ! incremaxlocal_P = max(incremaxlocal_P, maxval(abs(increments%producers)))

      ! -- Compute relaxation value --

      ! relax local = min(1, increobj/incremax)
      relaxlocal = min(1.d0, NewtonIncreObj_P/incremaxlocal_P) ! P
      ! if(NewtonIncreObj_P<incremaxlocal_P) &
      !    write(*,*) "Local pressure relaxation on", commRank, ":", incremaxlocal_P, "->", NewtonIncreObj_P, "relax=", relaxlocal

#ifdef _THERMIQUE_
      relaxlocal = min(relaxlocal, NewtonIncreObj_T/incremaxlocal_T) ! T
      ! if(NewtonIncreObj_T<incremaxlocal_T) &
      !    write(*,*) "Local pressure relaxation on", commRank, ":", incremaxlocal_T, "->", &
      !    NewtonIncreObj_T, "relax=", NewtonIncreObj_T/incremaxlocal_T
#endif

      do i = 1, NbPhase ! C_i^alpha
         do j = 1, NbComp

            if (MCP(j, i) == 1 .and. abs(incremaxlocal_C(j, i)) > eps) then
               relaxlocal = min(relaxlocal, NewtonIncreObj_C/incremaxlocal_C(j, i))
            end if
         end do
      end do

      do i = 1, NbPhase ! S^alpha
         if (abs(incremaxlocal_S(i)) > eps) then
            ! if(relaxlocal>NewtonIncreObj_S/incremaxlocal_S(i)) &
            !     write(*,*) "Relaxation induced by delta S max = ", incremaxlocal_S(i), " > ", NewtonIncreObj_S, " for phase ", i
            relaxlocal = min(relaxlocal, NewtonIncreObj_S/incremaxlocal_S(i))
         end if
      end do

      ! relax global
      call MPI_Allreduce(relaxlocal, relax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

      ! pmaxlocal = increobj_P/incremaxlocal_P
      ! !tmaxlocal = increobj_T/incremaxlocal_T
      ! do i=1, NbPhase
      !    smaxlocal(i) = increobj_S/incremaxlocal_S(i)
      ! end do

      ! call MPI_Allreduce(pmaxlocal, pmax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Allreduce(tmaxlocal, tmer": "compass",
      ! call MPI_Allreduce(smaxlocal, smax, NbPhase, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

      ! if(commRank==0) then
      !    write(*,'(F16.5,F16.5,F16.5)',advance="no"), pmax, smax(:)
      ! end if

   end subroutine Newton_compute_relaxation

   !This function  computes the relaxation only for MSWellNodes
   function Newton_compute_relaxation_mswells_C(increments_pointers) &
      result(relaxation) &
      bind(C, name="Newton_compute_relaxation_mswells")
      type(Newton_increments_pointers), intent(in), value :: increments_pointers
      real(c_double) :: relaxation
      type(Newton_increments) :: increments

      call Newton_pointers_to_values(increments_pointers, increments)
      call Newton_compute_relaxation_mswells(increments, relaxation)

   end function Newton_compute_relaxation_mswells_C

   !This function  get the last max increment only for MSWellNodes
   function Newton_get_last_max_inc_mswells_C() &
      result(dxmax) &
      bind(C, name="Newton_get_last_max_inc_mswells")
      real(c_double) :: dxmax

      dxmax = last_max_inc_mswells

   end function Newton_get_last_max_inc_mswells_C

   !> \brief Compute relaxation in Newton only for mswells.
   !!
   subroutine Newton_compute_max_increment_mswells(increments, dxmax)
      type(Newton_increments), intent(in) :: increments
      real(c_double), intent(out) :: dxmax

      double precision :: &
         incremaxlocal_P, &
         incremaxlocal_T, &
         incremaxlocal_C(NbComp, NbPhase), &
         incremaxlocal_S(NbPhase), &
         incremax_local(3), incremax_global(3), &
         relax_global, incremax_S

      integer :: k, i, ic, iph, icp, j, Ierr
      integer :: nbwells, s, num_s, unk_idx_s
      integer :: NbIncPTC

#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort("Multi-segmented wells are only implemented for diphasic physics!")

#else

      incremaxlocal_P = 0.d0

#ifdef _THERMIQUE_
      incremaxlocal_T = 0.d0
#endif
      incremaxlocal_C(:, :) = 0.d0
      incremaxlocal_S(:) = 0.d0

      nbwells = NbMSWellLocal_Ncpus(commRank + 1) !local mswells

      do k = 1, nbwells
         ! looping from  queue to the  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s)
            !Unkown indices of node s
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)

            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            ic = IncMSWell(s)%coats%ic
            NbIncPTC = NbIncPTC_ctx(ic)

            incremaxlocal_P = max(incremaxlocal_P, abs(increments%mswell_nodes(1, unk_idx_s)))

            do i = 2 + IndThermique, NbIncPTC
               icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
               iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

               incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(increments%mswell_nodes(i, unk_idx_s)))
            enddo

            do i = 1, NbPhasePresente_ctx(ic)
               iph = NumPhasePresente_ctx(i, ic)

               if (IncMSWell(s)%coats%ic .eq. DIPHASIC_CONTEXT) then
                  incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(increments%mswell_nodes(iph + NbIncPTC, unk_idx_s)))
               endif
            end do

#ifdef _THERMIQUE_
            incremaxlocal_T = max(incremaxlocal_T, abs(increments%mswell_nodes(2, unk_idx_s)))

#endif

         end do
      end do

      incremax_S = 0.d0
      do i = 1, NbPhase ! S^alpha
         incremax_S = max(incremax_S, incremaxlocal_S(i))
      end do

      incremax_local(1) = incremax_S
      incremax_local(2) = incremaxlocal_P
      incremax_local(3) = incremaxlocal_T

      ! Get all incremax
      call MPI_Allreduce(incremax_local, incremax_global, 3, MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

      !TODO: Should we implement relaxation-routine  as it is done for the Reservoir (previous function)?
      dxmax = incremax_global(1) + incremax_global(2)/1.d+5 + incremax_global(3)/1.d+2
      last_max_inc_mswells = dxmax !Save last maximum increment

#endif
   end subroutine Newton_compute_max_increment_mswells

   subroutine Newton_compute_relaxation_mswells(increments, relaxation)
      type(Newton_increments), intent(in) :: increments
      real(c_double) :: relaxation
      double precision ::   dxmax

      call Newton_compute_max_increment_mswells(increments, dxmax)
      relaxation = min(1.d0, 0.2d0/dxmax)

   end subroutine Newton_compute_relaxation_mswells

end module Newton
