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
      commRank, ComPASS_COMM_WORLD

   use DefModel, only: &
      NbIncTotalMax, NbPhase, NbComp, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IndThermique, MCP

   use MeshSchema, only: &
      IdNodeLocal, &
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbCellLocal_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus, &
      NbCellOwn_Ncpus, NbFracOwn_Ncpus, NbNodeOwn_Ncpus

   use NumbyContext, only: &
      NbIncPTC_ctx, NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   use SchemeParameters, only: &
      NewtonIncreObj_C, NewtonIncreObj_P, NewtonIncreObj_S, NewtonIncreObj_T, eps

   use IncCVReservoir, only: &
      IncNode, IncCell, IncFrac

   implicit none

   type, bind(C) :: Newton_increments_pointers
      type(c_ptr) :: nodes
      type(c_ptr) :: fractures
      type(c_ptr) :: cells
      type(c_ptr) :: injectors
      type(c_ptr) :: producers
   end type Newton_increments_pointers

   type Newton_increments
      real(c_double), pointer :: nodes(:, :)
      real(c_double), pointer :: fractures(:, :)
      real(c_double), pointer :: cells(:, :)
      real(c_double), pointer :: injectors(:)
      real(c_double), pointer :: producers(:)
   end type Newton_increments

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

   end subroutine Newton_pointers_to_values

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
      ! call MPI_Allreduce(tmaxlocal, tmax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Allreduce(smaxlocal, smax, NbPhase, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

      ! if(commRank==0) then
      !    write(*,'(F16.5,F16.5,F16.5)',advance="no"), pmax, smax(:)
      ! end if

   end subroutine Newton_compute_relaxation

end module Newton
