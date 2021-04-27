!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncCVWells

#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use mpi
   use CommonType
   use CommonMPI
   use DefModel
   use Physics
   use Thermodynamics
   use IncCVReservoir
   use MeshSchema
#else
! use, intrinsic :: iso_c_binding
   use iso_c_binding
   use ieee_arithmetic
   use mpi, only: MPI_Abort
   use CommonType, only: CSR
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use DefModel, only: NbComp, NbPhase, LIQUID_PHASE
   use Physics, only: gravity

   use Thermodynamics, only: f_DensiteMassique

   use IncCVReservoir, only: IncNode, NumPhasePresente_ctx, NbPhasePresente_ctx, IncCVReservoir_compute_density

   use MeshSchema, only: &
      XNodeLocal, &
      DataWellProdLocal, DataWellInjLocal, &
      NodebyWellProdLocal, NodebyWellInjLocal, &
      NodeDatabyWellProdLocal, NodeDatabyWellInjLocal, DataWellInjLocal, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus
#endif

   implicit none

   !> Type for the perforations, stores informations which are not constant in the well
   type, bind(C) :: WellPerforationState_type
      real(c_double) :: &
         Pression, & !< Pressure at the perforation
         Temperature, & !< Temperature at the perforation
         Density, & !< Density at the perforation: constant per edge, stored at node parent
         Saturation(NbPhase), & !< Phases saturation
         PressureDrop, & !< Pressure drop at the perforation, used to construct Pressure from the head pressure
         MolarFlowrate(NbComp), & !< Molar flux at the perforation, q_{w,s,i}
         EnergyFlowrate !< Energy flux at the perforation, q_{w,s,e}
   end type WellPerforationState_type

   type, bind(C) :: WellPerforations_type
      type(c_ptr) :: perforations_begin
      integer(c_size_t) :: nb_perforations
   end type WellPerforations_type

   ! well pressure for current time step
   real(c_double), allocatable, dimension(:), target, public :: &
      IncPressionWellInj, & !< Injection Well unknown: head pressure for current time step
      IncPressionWellProd !< Production Well unknown: head pressure for current time step

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   TYPE(WellPerforationState_type), allocatable, dimension(:), target, public :: &
      PerfoWellInj, & !< Injection Well informations at each perforation for current time step
      PerfoWellProd !< Production Well informations at each perforation for current time step

   ! well pressure from previous time step
   double precision, allocatable, dimension(:), target, public :: &
      IncPressionWellInjPreviousTimeStep, & !< Injection Well unknown: head pressure for previous time step
      IncPressionWellProdPreviousTimeStep !< Production Well unknown: head pressure for previous time step

   public :: &
      IncCVWells_allocate, &
      IncCVWells_free, &
      IncCVWells_PressureDropWellProd, &
      IncCVWells_PressureDropWellInj, &
      IncCVWells_estimate_producers_density, &
      IncCVWells_UpdatePressureDrop, &
      IncCVWells_UpdateWellPressures, &
      IncCVWells_NewtonIncrement, &
      IncCVWells_SaveIncPreviousTimeStep, &
      IncCVWells_LoadIncPreviousTimeStep, &
      get_producing_perforations, &
      get_injecting_perforations

contains

   function get_producing_perforations(well) result(perforations) &
      bind(C, name="get_producing_perforations")

      integer(c_size_t), intent(in), value :: well
      type(WellPerforations_type) :: perforations

      integer :: w

      w = well + 1 ! Fortran index starts at 1
      if (allocated(PerfoWellProd)) then
         perforations%perforations_begin = c_loc(PerfoWellProd(NodebyWellProdLocal%Pt(w) + 1))
         perforations%nb_perforations = NodebyWellProdLocal%Pt(w + 1) - NodebyWellProdLocal%Pt(w)
      else
         perforations%perforations_begin = c_null_ptr
         perforations%nb_perforations = 0
      end if

   end function get_producing_perforations

   function get_injecting_perforations(well) result(perforations) &
      bind(C, name="get_injecting_perforations")

      integer(c_size_t), intent(in), value :: well
      type(WellPerforations_type) :: perforations

      integer :: w

      w = well + 1 ! Fortran index starts at 1
      if (allocated(PerfoWellInj)) then
         perforations%perforations_begin = c_loc(PerfoWellInj(NodebyWellInjLocal%Pt(w) + 1))
         perforations%nb_perforations = NodebyWellInjLocal%Pt(w + 1) - NodebyWellInjLocal%Pt(w)
      else
         perforations%perforations_begin = c_null_ptr
         perforations%nb_perforations = 0
      end if

   end function get_injecting_perforations

   subroutine IncCVWells_set_density_from_reservoir(producer)
      integer, intent(in) :: producer

      integer :: s, nums
      do s = NodebyWellProdLocal%Pt(producer) + 1, NodebyWellProdLocal%Pt(producer + 1)
         nums = NodebyWellProdLocal%Num(s)
         PerfoWellProd(s)%Density = IncCVReservoir_compute_density(IncNode(nums))
      end do

   end subroutine IncCVWells_set_density_from_reservoir

   function IncCVWells_minimum_density(producer) result(rhomin)
      integer, intent(in) :: producer
      real(c_double) :: rhomin

      integer :: s

      rhomin = ieee_value(rhomin, ieee_positive_inf)
      do s = NodebyWellProdLocal%Pt(producer) + 1, NodebyWellProdLocal%Pt(producer + 1)
         rhomin = min(rhomin, PerfoWellProd(s)%Density)
      end do

   end function IncCVWells_minimum_density

   subroutine IncCVWells_use_minimum_density(producer)
      integer, intent(in) :: producer

      integer :: s
      real(c_double) :: rhomin

      rhomin = IncCVWells_minimum_density(producer)
      do s = NodebyWellProdLocal%Pt(producer) + 1, NodebyWellProdLocal%Pt(producer + 1)
         if (abs(NodeDatabyWellProdLocal%Val(s)%WID) < 1d-20) PerfoWellProd(s)%Density = rhomin
      end do

   end subroutine IncCVWells_use_minimum_density

   subroutine IncCVWells_estimate_producers_density(use_minimum_density) &
      bind(C, name="IncCVWells_estimate_producers_density")
      logical(c_bool), value, intent(in) :: use_minimum_density

      integer :: k

#ifndef NDEBUG
      if (NbWellProdLocal_Ncpus(commRank + 1) /= NodebyWellProdLocal%Nb) &
         call CommonMPI_abort("Well numbers are inconsistent")
#endif

      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)

         call IncCVWells_set_density_from_reservoir(k)
         if (use_minimum_density) call IncCVWells_use_minimum_density(k)

      end do

   end subroutine IncCVWells_estimate_producers_density

   !> \brief Compute well pressure drops and P_{w,s} using Pw (pressure head) and density for Well Producers
   subroutine IncCVWells_PressureDropWellProd
      integer :: k, s, nums, sparent
      double precision :: zp, zs

#ifndef NDEBUG
      if (NbWellProdLocal_Ncpus(commRank + 1) /= NodebyWellProdLocal%Nb) &
         call CommonMPI_abort("Well numbers are inconsistent")
#endif

      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)

         ! Check if the well is closed
         if (DataWellProdLocal(k)%IndWell == 'c') cycle

         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(k + 1), NodebyWellProdLocal%Pt(k) + 1, -1 !Reverse order, recall the numbering of parents & sons
            nums = NodebyWellProdLocal%Num(s)

            if (s == NodebyWellProdLocal%Pt(k + 1)) then ! head node, P = Pw

               PerfoWellProd(s)%PressureDrop = 0.d0

            else ! Pws = P_{w,parent} + \Delta P_{w,parent}

               zs = XNodeLocal(3, nums) ! z-cordinate of node s
               zp = XNodeLocal(3, NodeDatabyWellProdLocal%Val(s)%Parent) ! z-cordinate of parent of s
               sparent = NodeDatabyWellProdLocal%Val(s)%PtParent ! parent pointer
               PerfoWellProd(s)%PressureDrop = PerfoWellProd(sparent)%PressureDrop &
                                               + PerfoWellProd(sparent)%Density*gravity*(zp - zs)

            end if

         end do

      end do

   end subroutine IncCVWells_PressureDropWellProd

   !> \brief Compute well pressure drops and P_{w,s} using Pw (pressure head) and density for Well Injectors
   !! integration from node head (w) to node (s)
   subroutine IncCVWells_PressureDropWellInj

      integer :: s, sp, n, k, nbwells, nums, nump
      real(c_double) :: Pw_head, Pws, zp, zs, T, C(NbComp)

      ! nb pieces for discrete integration
      ! FIXME: call quad or something similar
      integer, parameter :: Npiece = 100
      real(c_double) :: dz, Rhotmp, dPf, dTf, dCf(NbComp)
      ! FIXME: this is temporary a work array, its size could be much smaller
      real(c_double), dimension(NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb + 1)) :: parent_pressure

      nbwells = NbWellInjLocal_Ncpus(commRank + 1)

#ifndef NDEBUG
      if (NbWellInjLocal_Ncpus(commRank + 1) /= NodebyWellInjLocal%Nb) &
         call CommonMPI_abort("Well numbers are inconsistent")
#endif

      do k = 1, nbwells

         ! Check if the well is closed
         if (DataWellInjLocal(k)%IndWell == 'c') cycle

         T = DataWellInjLocal(k)%InjectionTemperature
         C = DataWellInjLocal(k)%CompTotal

         ! integrate over head to queue
         s = NodebyWellInjLocal%Pt(k + 1)

#ifndef NDEBUG
         if (s <= NodebyWellInjLocal%Pt(k)) &
            call CommonMPI_abort("Well has no nodes!")
         if (NodeDatabyWellInjLocal%Val(s)%Parent /= -1) &
            call CommonMPI_abort("Inconsistent well head node!")
#endif

         parent_pressure(s) = IncPressionWellInj(k)
         PerfoWellInj(s)%PressureDrop = 0.d0

         do while (s > NodebyWellInjLocal%Pt(k) + 1)
            s = s - 1
            nums = NodebyWellInjLocal%Num(s)
            nump = NodeDatabyWellInjLocal%Val(s)%Parent
            zs = XNodeLocal(3, nums) ! z-cordinate of node s
            zp = XNodeLocal(3, nump) ! z-cordinate of parent of s
            dz = (zp - zs)/Npiece
            sp = NodeDatabyWellInjLocal%Val(s)%PtParent ! parent pointer
            Pws = parent_pressure(sp)
            ! integrate from zp to zs
            do n = 1, Npiece
               call f_DensiteMassique(LIQUID_PHASE, Pws, T, C, Rhotmp, dPf, dTf, dCf)
               Pws = Pws + gravity*Rhotmp*dz
            end do
            parent_pressure(s) = Pws
            PerfoWellInj(s)%PressureDrop = PerfoWellInj(sp)%PressureDrop + Pws - parent_pressure(sp)
         end do

      end do

   end subroutine IncCVWells_PressureDropWellInj

   subroutine IncCVWells_UpdatePressureDrop() &
      bind(C, name="IncCVWells_UpdatePressureDrop")

      call IncCVWells_PressureDropWellInj
      call IncCVWells_PressureDropWellProd

   end subroutine IncCVWells_UpdatePressureDrop

   !> \brief Update  only the  pressures of all wells, while the well pressures drops are kept constant.
   subroutine IncCVWells_UpdateWellPressures() &
      bind(C, name="IncCVWells_UpdateWellPressures")

      call IncCVWells_UpdateWellPressures_loop(NodebyWellInjLocal, IncPressionWellInj, PerfoWellInj)
      call IncCVWells_UpdateWellPressures_loop(NodebyWellProdLocal, IncPressionWellProd, PerfoWellProd)

   end subroutine IncCVWells_UpdateWellPressures

   !> \brief Update well pressures actual loop implementation
   subroutine IncCVWells_UpdateWellPressures_loop(wellnodes, reference_pressure, perforations)
      type(CSR), intent(in) :: wellnodes
      real(c_double), dimension(:), intent(in) :: reference_pressure
      type(WellPerforationState_type), dimension(:), intent(inout) :: perforations

      integer :: k, nbwells, s, s_head
      double precision :: Pw_head

      nbwells = wellnodes%Nb
#ifndef NDEBUG
      if (size(reference_pressure) /= nbwells) &
         call CommonMPI_abort("Inconsistent well data.")
#endif
      do k = 1, nbwells
         s_head = wellnodes%Pt(k + 1)
         if (wellnodes%Pt(k) < s_head) then
#ifndef NDEBUG
            if (abs(perforations(s_head)%PressureDrop) > 0.d0) &
               call CommonMPI_abort("There should be no pressure drop at well head.")
#endif
            Pw_head = reference_pressure(k)
            do s = s_head, wellnodes%Pt(k) + 1, -1 !Reverse order, recall the numbering of parents & sons
               perforations(s)%Pression = Pw_head + perforations(s)%PressureDrop
            end do
#ifndef NDEBUG
         else
            call CommonMPI_abort("Well without nodes.")
#endif
         endif
      end do

   end subroutine IncCVWells_UpdateWellPressures_loop

   !> \brief Allocate well unknowns vectors
   subroutine IncCVWells_allocate() &
      bind(C, name="IncCVWells_allocate")

      integer :: Nb, Nnz

      allocate (IncPressionWellInj(NbWellInjLocal_Ncpus(commRank + 1)))
      allocate (IncPressionWellProd(NbWellProdLocal_Ncpus(commRank + 1)))

      Nb = NodebyWellInjLocal%Nb
      Nnz = NodebyWellInjLocal%Pt(Nb + 1)
      allocate (PerfoWellInj(Nnz))

      Nb = NodebyWellProdLocal%Nb
      Nnz = NodebyWellProdLocal%Pt(Nb + 1)
      allocate (PerfoWellProd(Nnz))

      allocate (IncPressionWellInjPreviousTimeStep(NbWellInjLocal_Ncpus(commRank + 1)))
      allocate (IncPressionWellProdPreviousTimeStep(NbWellProdLocal_Ncpus(commRank + 1)))

   end subroutine IncCVWells_allocate

   !> \brief Deallocate well unknowns vectors
   subroutine IncCVWells_free() &
      bind(C, name="IncCVWells_free")

      deallocate (IncPressionWellInj)
      deallocate (IncPressionWellProd)

      deallocate (PerfoWellInj)
      deallocate (PerfoWellProd)

      deallocate (IncPressionWellInjPreviousTimeStep)
      deallocate (IncPressionWellProdPreviousTimeStep)

   end subroutine IncCVWells_free

   subroutine IncCVWells_NewtonIncrement( &
      NewtonIncreWellInj, NewtonIncreWellProd, relax)

      double precision, dimension(:), intent(in) :: &
         NewtonIncreWellInj, &
         NewtonIncreWellProd

      double precision, intent(in) :: relax

      integer :: k

      ! injection wells (head Pressure)
      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)
         IncPressionWellInj(k) = IncPressionWellInj(k) + relax*NewtonIncreWellInj(k)
      end do

      ! production wells
      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)
         IncPressionWellProd(k) = IncPressionWellProd(k) + relax*NewtonIncreWellProd(k)
      end do

   end subroutine IncCVWells_NewtonIncrement
   !> \brief Load previous status to start again current time iteration.
   !!
   !! Copy IncObjPreviousTimeStep to IncObj
   subroutine IncCVWells_LoadIncPreviousTimeStep

      integer :: k

      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)
         IncPressionWellInj(k) = IncPressionWellInjPreviousTimeStep(k)
      end do
      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)
         IncPressionWellProd(k) = IncPressionWellProdPreviousTimeStep(k)
      end do

   end subroutine IncCVWells_LoadIncPreviousTimeStep

   subroutine IncCVWells_SaveIncPreviousTimeStep

      integer :: k

      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)
         IncPressionWellInjPreviousTimeStep(k) = IncPressionWellInj(k)
      end do
      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)
         IncPressionWellProdPreviousTimeStep(k) = IncPressionWellProd(k)
      end do

   end subroutine IncCVWells_SaveIncPreviousTimeStep

end module IncCVWells
