module IncCVMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use DefModel
   use mpi
   use CommonMPI
   use IncCVReservoir
   use MeshSchema
   use MeshSchemaMSWells
#else
   use iso_c_binding, only: c_double, c_bool

   use DefModel, only: &
      NbPhase, NbComp, NbContexte
#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      LIQUID_PHASE, GAS_PHASE, LIQUID_CONTEXT, GAS_CONTEXT, DIPHASIC_CONTEXT
#endif

   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, IncNode, NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IncCVReservoir_NewtonIncrement_reservoir

   use MeshSchema, only: &
      XNodeLocal, &
      NodebyMSWellLocal, &
      NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus

   use MeshSchemaMSWells, only: &
      IncIdxNodebyMSWellLocal, NbMSWellNodeLocal_Ncpus
#endif
   implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set this to one to print mswells info like pressure, temperature, etc
! See function IncCVMSWells_print_info_to_file below
#define _DEBUG_INCCVM_MSWELLS_ 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   type, bind(C) :: TYPE_IncCVMSWells
      !Unkowns per Node
      type(TYPE_IncCVReservoir) :: coats
   end type TYPE_IncCVMSWells

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   TYPE(TYPE_IncCVMSWells), allocatable, dimension(:), target, public :: &
      IncMSWell   !< Producer & Injector   MSWell unkowns

   ! Inc for previous time step: mswell_nodes
   type(TYPE_IncCVReservoir), allocatable, dimension(:), target, public :: IncCoatsMSWellPreviousTimeStep

   !Experimental:
   !    Flag to compute coupling terms between the reservoir and mswells.
   !    This is used in the Jacobian, JacobianMSWells, and Residu  modules
   !    By default has value equel to .false.
   logical(c_bool), target, private :: compute_coupling_mswells

   public :: &
      IncCVMSWells_allocate, &
      IncCVMSWells_free, &
      IncCVMSWells_compute_coupling, &
      IncCVMSWells_set_compute_coupling, &
      IncCVMSWells_NewtonIncrement, &
      IncCVMSWells_print_info_to_file
contains

   !> \brief Allocate well unknowns vectors

   subroutine IncCVMSWells_allocate()

      integer :: Nb, Nnz

      Nb = NodebyMSWellLocal%Nb
      Nnz = NodebyMSWellLocal%Pt(Nb + 1)
      allocate (IncMSWell(Nnz))
      allocate (IncCoatsMSWellPreviousTimeStep(Nnz))
      compute_coupling_mswells = .false.
#if _DEBUG_INCCVM_MSWELLS_ == 1
      if (commRank == 0) then !Master proc
         open (unit=110, file='mswell.val', status='unknown', action='write')
         write (110, *) "# x-coord, z-coord, w_gas_sat,  w_pression,  w_T,  w_context"
         close (110)

         open (unit=310, file='mswell_exd.val', status='unknown', action='write')
         write (310, *) "# w_total_gas_sat"
         close (310)
      endif
#endif
   end subroutine IncCVMSWells_allocate

   !> \brief Deallocate well unknowns vectors
   subroutine IncCVMSWells_free()

      deallocate (IncCoatsMSWellPreviousTimeStep)
      deallocate (IncMSWell)

   end subroutine IncCVMSWells_free

   !> \brief Save current status if it is necessary to start again current time iteration.
   !!
   !! Copy IncObj to IncObjPreviousTimeStep
   subroutine IncCMSWells_SaveIncPreviousTimeStep

      integer :: k

      ! save current status
      do k = 1, NbMSWellNodeLocal_Ncpus(commRank + 1)
         IncCoatsMSWellPreviousTimeStep(k) = IncMSWell(k)%coats
      end do

   end subroutine IncCMSWells_SaveIncPreviousTimeStep

   !! Load IncObj to IncObjPreviousTimeStep
   subroutine IncCMSWells_LoadIncPreviousTimeStep

      integer :: k

      ! save current status
      do k = 1, NbMSWellNodeLocal_Ncpus(commRank + 1)
         IncMSWell(k)%coats = IncCoatsMSWellPreviousTimeStep(k)
      end do

   end subroutine IncCMSWells_LoadIncPreviousTimeStep

   !Copy states from reservoir and init some data at well nodes for all mswells
   subroutine IncCVMSWells_copy_states_from_reservoir() &
      bind(C, name="IncCVMSWells_copy_states_from_reservoir")

      integer :: s, k, nums, nbwells
#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort( &
         "In function IncCVMSWells_copy_states_from_reservoir: Multi-segmented wells are only implemented for diphasic physics!")

#else

      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do k = 1, nbwells

         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
            nums = NodebyMSWellLocal%Num(s)

            IncMSWell(s)%coats = IncNode(nums)
            !IncMSWell(s)%coats%Pression = IncMSWell(s)%coats%Pression
            ! FIXME: this a hard coded trick to make the producer produce
            !        at the beginning of the simulation by lowering slightly
            !        the well pressure below the reservoir pressure
            IncMSWell(s)%coats%Pression = IncMSWell(s)%coats%Pression - 1
            !Set phase pressure for all phases
            IncMSWell(s)%coats%phase_pressure(:) = IncMSWell(s)%coats%Pression

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!Init to diphasic
            !!!IncMSWell(s)%coats%ic = 3 !Diphasic
            !!!IncMSWell(s)%coats%Saturation(1) = 0.8d0!Sg
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IncMSWell(s)%coats%Saturation(2) = 1.d0 - IncMSWell(s)%coats%Saturation(1) !Sl

            !Init state
            IncMSWell(s)%coats%AccVol = 0.d0
         end do

      end do

#endif

   end subroutine IncCVMSWells_copy_states_from_reservoir

!!   !> \brief Loops over mswell_nodes to increment the unknowns
   subroutine IncCVMSWells_NewtonIncrement( &
      NewtonIncreMSWellNode, &
      relax)

      double precision, dimension(:, :), intent(in) :: &
         NewtonIncreMSWellNode

      double precision, intent(in) :: relax
      integer :: k, nbwells, s, unk_idx_s

      nbwells = NbMSWellLocal_Ncpus(commRank + 1) !local mswells

      do k = 1, nbwells
         ! looping from  queue to the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            !Unkown indices of node s
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)
            call IncCVReservoir_NewtonIncrement_reservoir(IncMSWell(s)%coats, NewtonIncreMSWellNode(:, unk_idx_s), relax, .true.)

         end do
      end do

   end subroutine IncCVMSWells_NewtonIncrement

   function IncCVMSWells_compute_coupling() result(flag)

      logical(c_bool) :: flag

      flag = compute_coupling_mswells

   end function IncCVMSWells_compute_coupling

   subroutine IncCVMSWells_set_compute_coupling(flag) &
      bind(C, name="IncCVMSWells_set_compute_coupling")

      logical(c_bool), intent(in), value:: flag

      compute_coupling_mswells = flag

   end subroutine IncCVMSWells_set_compute_coupling

   subroutine IncCVMSWells_print_info_to_file() &
      bind(C, name="IncCVMSWells_print_info_to_file")

      integer :: s, k, nbwells, num_s, num_ss
      double precision :: dz, total_gas_sat

#if _DEBUG_INCCVM_MSWELLS_ == 1

#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort( &
         "In function IncCVMSWells_print_info_to_file: Multi-segmented wells are only implemented for diphasic physics!")

#else
      if (commRank == 0) then !Master proc, watch it: this only works if the proc 0 has a local mswell

         open (unit=110, file='mswell.val', status='old', position='append', action='write')
         open (unit=310, file='mswell_exd.val', status='old', position='append', action='write')

         nbwells = NbMSWellLocal_Ncpus(commRank + 1)
         do k = 1, nbwells

            ! Check if the well is closed
            if (DataMSWellLocal(k)%IndWell == 'c') cycle

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!Compute total gas sat
            total_gas_sat = 0.d0

            ! looping  over the nodes stricly between the  queue to  head
            do s = NodebyMSWellLocal%Pt(k) + 2, NodebyMSWellLocal%Pt(k + 1) - 1
               num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
               num_ss = NodebyMSWellLocal%Num(s - 1)
               dz = XNodeLocal(3, num_s) - XNodeLocal(3, num_ss)
               total_gas_sat = total_gas_sat + dz*IncMSWell(s)%coats%Saturation(GAS_PHASE)

            end do
            !add bottom
            s = NodebyMSWellLocal%Pt(k) + 1
            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_ss = NodebyMSWellLocal%Num(s + 1)  ! next well node
            dz = XNodeLocal(3, num_ss) - XNodeLocal(3, num_s)
            total_gas_sat = total_gas_sat + 0.5*dz*IncMSWell(s)%coats%Saturation(GAS_PHASE)
            !add head
            s = NodebyMSWellLocal%Pt(k + 1)
            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_ss = NodebyMSWellLocal%Num(s - 1)  ! previous well node
            dz = XNodeLocal(3, num_s) - XNodeLocal(3, num_ss)
            total_gas_sat = total_gas_sat + 0.5*dz*IncMSWell(s)%coats%Saturation(GAS_PHASE)

            write (310, *) total_gas_sat

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!Print info

            ! looping from  queue to  head
            do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
               num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)

               write (110, *) XNodeLocal(1, num_s), XNodeLocal(3, num_s), &
                  IncMSWell(s)%coats%Saturation(GAS_PHASE), &
                  IncMSWell(s)%coats%Pression, &
                  IncMSWell(s)%coats%Temperature, &
                  IncMSWell(s)%coats%ic

            end do

         end do
         write (110, *) " "
         close (110)
         write (310, *) " "
         close (310)

      end if
#endif
#endif
   end subroutine IncCVMSWells_print_info_to_file

end module IncCVMSWells
