!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncCVWells

   use MeshSchema
   use DefModel
   use Thermodynamics

   use NumbyContext
   use CommonMPI
   use Physics
   use SchemeParameters
   use IncCVReservoir

   use iso_c_binding

   implicit none

   !> Type for the perforations, stores informations which are not constant in the well
   TYPE TYPE_PhysPerfoWell
      double precision :: &
         Pression, & !< Pressure at the perforation
         Temperature, & !< Temperature at the perforation
         Density, & !< Density at the perforation: constant per edge, stored at node parent
         PressureDrop !< Pressure drop at the perforation, used to construct Pressure from the head pressure
      ! FluxMolar(NbComp), & !< Molar flux at the perforation, q_{w,s,i}
      ! FluxEnergy           !< Energy flux at the perforation, q_{w,s,e}
   end TYPE TYPE_PhysPerfoWell

   ! well pressure for current time step
   real(c_double), allocatable, dimension(:), target, public :: &
      IncPressionWellInj, & !< Injection Well unknown: head pressure for current time step
      IncPressionWellProd !< Production Well unknown: head pressure for current time step

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   TYPE(TYPE_PhysPerfoWell), allocatable, dimension(:), target, public :: &
      PerfoWellInj, & !< Injection Well informations at each perforation for current time step
      PerfoWellProd !< Production Well informations at each perforation for current time step

   ! well pressure from previous time step
   double precision, allocatable, dimension(:), target, public :: &
      IncPressionWellInjPreviousTimeStep, & !< Injection Well unknown: head pressure for previous time step
      IncPressionWellProdPreviousTimeStep !< Production Well unknown: head pressure for previous time step

   ! Sorted NodebyWellInjLocal%Num and %Val according to z-coordinate for each well
   integer, allocatable, dimension(:), protected :: ZSortedInj_Znum
   double precision, allocatable, dimension(:), protected :: ZSortedInj_Zval

   public :: &
      IncCVWells_allocate, &
      IncCVWells_free, &
      IncCVWells_SortHeightWellInj, &
      IncCVWells_PressureDropWellProd, &
      IncCVWells_PressureDropWellInj, &
      IncCVWells_PressureDrop, &
      IncCVWells_PressureDropWellInj_integrate, &
      IncCVWells_NewtonIncrement, &
      IncCVWells_SaveIncPreviousTimeStep, &
      IncCVWells_LoadIncPreviousTimeStep

contains

   ! sort the nodes of wells by z-cordinate from the smallest to the largest
   ! the results are stored in ZSortedInj_Znum (num) and in ZSortedinj_Zval (z-cordinate)
   subroutine IncCVWells_SortHeightWellInj

      integer :: s, k, j, Nnz, nums
      integer :: tmp_num
      double precision :: tmp_val

      Nnz = NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb + 1)
      do s = 1, Nnz
         nums = NodebyWellInjLocal%Num(s)
         ZSortedInj_Znum(s) = s
         ZSortedInj_Zval(s) = XNodeLocal(3, nums)
      end do

      do k = 1, NodebyWellInjLocal%Nb ! = NbWellInjLocal(commRank+1)

         do s = 1, NodebyWellInjLocal%Pt(k + 1) - NodebyWellInjLocal%Pt(k)
            do j = NodebyWellInjLocal%Pt(k) + 1, NodebyWellInjLocal%Pt(k + 1) - s

               if (ZSortedInj_Zval(j) > ZSortedInj_Zval(j + 1)) then

                  tmp_num = ZSortedInj_Znum(j + 1)
                  tmp_val = ZSortedInj_Zval(j + 1)

                  ZSortedInj_Znum(j + 1) = ZSortedInj_Znum(j)
                  ZSortedInj_Zval(j + 1) = ZSortedInj_Zval(j)

                  ZSortedInj_Znum(j) = tmp_num
                  ZSortedInj_Zval(j) = tmp_val
               end if
            end do
         end do
      end do

   end subroutine IncCVWells_SortHeightWellInj

   ! compute P_{w,s} using Pw (pressure head) and density
   subroutine IncCVWells_PressureDropWellProd

      integer :: k, s, m, mph, nums, sparent
      double precision :: Pws, zp, zs, Pdrop, Rhotmp
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)

      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)

         ! ! Init PressureDrop as zero
         ! do s=NodebyWellProdLocal%Pt(k)+1, NodebyWellProdLocal%Pt(k+1)
         !    PerfoWellProd(s)%PressureDrop = 0.d0
         ! end do

         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(k + 1), NodebyWellProdLocal%Pt(k) + 1, -1
            nums = NodebyWellProdLocal%Num(s)

            ! average density
            PerfoWellProd(s)%Density = 0.d0
            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic)
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

               call f_DensiteMolaire(NbPhase, IncNode(nums)%Pression, IncNode(nums)%Temperature, &
                                     IncNode(nums)%Comp(:, mph), IncNode(nums)%Saturation, Rhotmp, dPf, dTf, dCf, dSf)
               PerfoWellProd(s)%Density = PerfoWellProd(s)%Density + Rhotmp*IncNode(nums)%Saturation(m)
            end do

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

   subroutine IncCVWells_PressureDropWellInj_integrate(sfirst, slast, direction, Pfirst, T, C)

      integer, intent(in) :: sfirst, slast
      ! direction = 1 upwards / -1 downwards
      integer, intent(in) :: direction
      ! T and C are constant !!!
      double precision, intent(in) :: Pfirst, T, C(NbComp)

      integer :: n, s
      ! nb pieces for discrete integration
      ! FIXME: call quad or something similar
      integer, parameter :: Npiece = 100
      double precision :: Ptmp, Stmp(NbPhase), &
         z1, z2, dz, Pdrop, Rhotmp, dPf, dTf, dCf(NbComp), dSf(NbPhase)
#ifndef NDEBUG
      integer :: Ierr, errcode ! used for MPI_Abort
#endif

      Ptmp = Pfirst
      ! Saturation is fixed to liquid
      Stmp(:) = 0.d0
      Stmp(LIQUID_PHASE) = 1.d0
      PerfoWellInj(sfirst)%Pression = Pfirst
      do s = sfirst, slast - direction, direction
         z1 = ZSortedInj_Zval(s)
         z2 = ZSortedInj_Zval(s + direction)
#ifndef NDEBUG
         if (direction*(z2 - z1) < 0) then
            write (*, *) 'Nodes are badly sorted.'
            call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
         end if
#endif
         dz = (z2 - z1)/Npiece
         do n = 1, Npiece
            call f_DensiteMolaire(LIQUID_PHASE, Ptmp, T, C, Stmp, &
                                  Rhotmp, dPf, dTf, dCf, dSf)
            Ptmp = Ptmp - direction*gravity*Rhotmp*dz
         end do
         PerfoWellInj(s + direction)%Pression = Ptmp
      end do

   end subroutine IncCVWells_PressureDropWellInj_integrate

   ! compute pressure drop of injection well
   ! integration from node head (w) to node (s)
   subroutine IncCVWells_PressureDropWellInj

      integer :: s, wk, nbwells, bottom, head, headsortedpos
      double precision :: Phead, T, C(NbComp)

      nbwells = NbWellInjLocal_Ncpus(commRank + 1)
      do s = NodebyWellInjLocal%Pt(1) + 1, NodebyWellInjLocal%Pt(nbwells + 1)
         PerfoWellInj(s)%PressureDrop = 0.d0
      end do
      do wk = 1, nbwells
         ! Locate head
         bottom = NodebyWellInjLocal%Pt(wk) + 1
         head = NodebyWellInjLocal%Pt(wk + 1)
         do headsortedpos = bottom, head
            if (ZsortedInj_Znum(headsortedpos) == head) exit
         end do
         Phead = IncPressionWellInj(wk)
         T = DataWellInjLocal(wk)%Temperature
         C = DataWellInjLocal(wk)%CompTotal
         call IncCVWells_PressureDropWellInj_integrate(headsortedpos, head, 1, Phead, T, C)
         call IncCVWells_PressureDropWellInj_integrate(headsortedpos, bottom, -1, Phead, T, C)
         do s = bottom, head
            PerfoWellInj(s)%PressureDrop = Phead - PerfoWellInj(s)%Pression
         end do
      end do

   end subroutine IncCVWells_PressureDropWellInj

   subroutine IncCVWells_PressureDrop() &
      bind(C, name="IncCVWells_PressureDrop")

      call IncCVWells_PressureDropWellInj
      call IncCVWells_PressureDropWellProd

   end subroutine IncCVWells_PressureDrop

   subroutine IncCVWells_allocate

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

      Nnz = NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb + 1)
      allocate (ZSortedInj_Znum(Nnz))
      allocate (ZSortedInj_Zval(Nnz))

      ! Sort injection well
      ! and save the results in vector ZsortedInj_Znum and ZsortedInj_Zval
      call IncCVWells_SortHeightWellInj

   end subroutine IncCVWells_allocate

   !> \brief Deallocate unknowns vectors
   subroutine IncCVWells_free

      deallocate (IncPressionWellInj)
      deallocate (IncPressionWellProd)

      deallocate (PerfoWellInj)
      deallocate (PerfoWellProd)

      deallocate (IncPressionWellInjPreviousTimeStep)
      deallocate (IncPressionWellProdPreviousTimeStep)

      deallocate (ZSortedInj_Znum)
      deallocate (ZSortedInj_Zval)

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
