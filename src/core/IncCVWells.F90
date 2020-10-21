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
   use CommonMPI
   use DefModel
   use Physics
   use Thermodynamics
   use IncCVReservoir
   use MeshSchema
#else
! use, intrinsic :: iso_c_binding
   use iso_c_binding
   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use DefModel, only: NbComp, NbPhase, LIQUID_PHASE
   use Physics, only: gravity

   use Thermodynamics, only: f_DensiteMassique

   use IncCVReservoir, only: IncNode, NumPhasePresente_ctx, NbPhasePresente_ctx

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
      IncCVWells_InitPressureDrop, &
      IncCVWells_UpdatePressureDrop, &
      IncCVWells_UpdateProdWellPressures, &
      IncCVWells_UpdateInjWellPressures, &
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

   !> \brief Compute well pressure drops and P_{w,s} using Pw (pressure head) and density for Well Producers
   subroutine IncCVWells_PressureDropWellProd(use_avg_dens)
      !Flag to use the an avg_density  from the reservoir or the density computed from  the function DefFlashWells_TimeFlash  to compute the pressure drops
      !We want to use the avg_density usually to init the pressure drops
      logical, intent(in) :: use_avg_dens

      integer :: k, s, m, mph, nums, sparent
      double precision :: Pws, zp, zs, Pdrop, Rhotmp
      double precision :: dPf, dTf, dCf(NbComp)

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

            if (use_avg_dens) then
               ! average density
               PerfoWellProd(s)%Density = 0.d0
               do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic)
                  mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

                  call f_DensiteMassique(mph, IncNode(nums)%Pression, IncNode(nums)%Temperature, &
                                         IncNode(nums)%Comp(:, mph), Rhotmp, dPf, dTf, dCf)
                  PerfoWellProd(s)%Density = PerfoWellProd(s)%Density + Rhotmp*IncNode(nums)%Saturation(mph)
               end do
            end if

            if (s == NodebyWellProdLocal%Pt(k + 1)) then ! head node, P = Pw

               Pws = IncPressionWellProd(k) ! P_{w,s} = Pw
               PerfoWellProd(s)%Pression = Pws
               PerfoWellProd(s)%PressureDrop = 0.d0

            else ! Pws = P_{w,parent} + \Delta P_{w,parent}

               zs = XNodeLocal(3, nums) ! z-cordinate of node s
               zp = XNodeLocal(3, NodeDatabyWellProdLocal%Val(s)%Parent) ! z-cordinate of parent of s

               sparent = NodeDatabyWellProdLocal%Val(s)%PtParent ! parent pointer

               Pdrop = PerfoWellProd(sparent)%Density*gravity*(zp - zs)
               Pws = PerfoWellProd(sparent)%Pression + Pdrop ! Pws

               PerfoWellProd(s)%Pression = Pws
               PerfoWellProd(s)%PressureDrop = PerfoWellProd(sparent)%PressureDrop + Pdrop
            end if

         end do
      end do

   end subroutine IncCVWells_PressureDropWellProd

   !> \brief Compute well pressure drops and P_{w,s} using Pw (pressure head) and density for Well Injectors
   !! integration from node head (w) to node (s)
   subroutine IncCVWells_PressureDropWellInj

      integer :: s, sp, n, k, nbwells, nums, nump
      double precision :: Pw_head, Pws, zp, zs, Pdrop, T, C(NbComp)

      ! nb pieces for discrete integration
      ! FIXME: call quad or something similar
      integer, parameter :: Npiece = 100
      double precision :: Ptmp, z1, z2, dz, Rhotmp, dPf, dTf, dCf(NbComp)

      nbwells = NbWellInjLocal_Ncpus(commRank + 1)

#ifndef NDEBUG
      if (NbWellInjLocal_Ncpus(commRank + 1) /= NodebyWellInjLocal%Nb) &
         call CommonMPI_abort("Well numbers are inconsistent")
#endif

      do k = 1, nbwells

         ! Check if the well is closed
         if (DataWellInjLocal(k)%IndWell == 'c') cycle

         Pw_head = IncPressionWellInj(k)
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

         PerfoWellInj(s)%Pression = Pw_head
         PerfoWellInj(s)%PressureDrop = 0.d0

         do while (s > NodebyWellInjLocal%Pt(k) + 1)
            s = s - 1
            nums = NodebyWellInjLocal%Num(s)
            nump = NodeDatabyWellInjLocal%Val(s)%Parent
            zs = XNodeLocal(3, nums) ! z-cordinate of node s
            zp = XNodeLocal(3, nump) ! z-cordinate of parent of s
            dz = (zp - zs)/Npiece
            sp = NodeDatabyWellInjLocal%Val(s)%PtParent ! parent pointer
            Pws = PerfoWellInj(sp)%Pression
            ! integrate from zp to zs
            do n = 1, Npiece
               call f_DensiteMassique(LIQUID_PHASE, Pws, T, C, Rhotmp, dPf, dTf, dCf)
               Pws = Pws + gravity*Rhotmp*dz
            end do
            PerfoWellInj(s)%Pression = Pws
            PerfoWellInj(s)%PressureDrop = Pws - Pw_head
         end do

      end do

   end subroutine IncCVWells_PressureDropWellInj

   subroutine IncCVWells_InitPressureDrop() &
      bind(C, name="IncCVWells_InitPressureDrop")

      call IncCVWells_PressureDropWellInj
      call IncCVWells_PressureDropWellProd(.true.)

   end subroutine IncCVWells_InitPressureDrop

   subroutine IncCVWells_UpdatePressureDrop() &
      bind(C, name="IncCVWells_UpdatePressureDrop")

      call IncCVWells_PressureDropWellInj
      call IncCVWells_PressureDropWellProd(.false.)

   end subroutine IncCVWells_UpdatePressureDrop

   !   !> \brief Update well Pressures of injection well,
   subroutine IncCVWells_UpdateInjWellPressures
      integer :: s, n, k, nbwells
      double precision :: Pw_head

      nbwells = NbWellInjLocal_Ncpus(commRank + 1)
      do k = 1, nbwells
         Pw_head = IncPressionWellInj(k)
         ! looping from head to queue
         do s = NodebyWellInjLocal%Pt(k + 1), NodebyWellInjLocal%Pt(k) + 1, -1 !Reverse order, recall the numbering of parents & sons
            if (s == NodebyWellInjLocal%Pt(k + 1)) then ! head node, P = Pw
               PerfoWellInj(s)%Pression = Pw_head
            else ! explicit computation
               PerfoWellInj(s)%Pression = Pw_head + PerfoWellInj(s)%PressureDrop

            end if
         end do
      end do
   end subroutine IncCVWells_UpdateInjWellPressures

   !> \brief Update  well pressures of producer well
   subroutine IncCVWells_UpdateProdWellPressures

      integer :: s, n, k, nbwells
      double precision :: Pw_head

      nbwells = NbWellProdLocal_Ncpus(commRank + 1)
      do k = 1, nbwells
         Pw_head = IncPressionWellProd(k)
         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(k + 1), NodebyWellProdLocal%Pt(k) + 1, -1 !Reverse order, recall the numbering of parents & sons
            if (s == NodebyWellProdLocal%Pt(k + 1)) then ! head node, P = Pw
               PerfoWellProd(s)%Pression = Pw_head
            else ! explicit computation
               PerfoWellProd(s)%Pression = Pw_head + PerfoWellProd(s)%PressureDrop

            end if
         end do
      end do
   end subroutine IncCVWells_UpdateProdWellPressures

   !> \brief Update  only the  pressures of all wells, while the well pressures drops are kept constant.
   subroutine IncCVWells_UpdateWellPressures() &
      bind(C, name="IncCVWells_UpdateWellPressures")

      call IncCVWells_UpdateInjWellPressures
      call IncCVWells_UpdateProdWellPressures

   end subroutine IncCVWells_UpdateWellPressures

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
