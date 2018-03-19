!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp, context switch

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use Thermodynamics
   use IncCV
   use VAGFrac ! for rocktypes
   use DefFlashWells

   implicit none

   public :: &
      DefFlash_Flash ! Flash after each Newton iteration

   private :: &
      DefFlash_Flash_cv

contains

   !> \brief Main surboutine, after each Newton iteration
   !! execute the flash to determine the phases
   !! which are actualy present, and
   !! the mode of the well (flowrate or pressure).
   subroutine DefFlash_Flash

      integer :: k
      double precision :: Psat, dTsat, Tsat

      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncNode(k), NodeRocktypeLocal(:, k), PoroVolDarcyNode(k))
      end do

      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncFrac(k), FracRocktypeLocal(:, k), PoroVolDarcyFrac(k))
      end do

      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         call DefFlash_Flash_cv(IncCell(k), CellRocktypeLocal(:, k), PoroVolDarcyCell(k))
      end do

      ! choose between linear or non-linear update of the Newton unknown Pw
      ! The next subroutines also compute the mode of the wells ('pressure' or 'flowrate')
      call DefFlashWells_NewtonFlashLinWells

   end subroutine DefFlash_Flash

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume ?????
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

      type(Type_IncCV), intent(inout) :: inc
      INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      integer :: i, iph, j, icp, m, mph, ic
      double precision :: DensiteMolaire(NbComp), acc1, acc2, &
         dPf, dTf, dCf(NbComp), dSf(NbPhase)

      double precision :: Tsat, dTsatdP, Psat, dPsatdT

      ic = inc%ic

      if (ic == 1) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression > Psat) then
            inc%ic = 3
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(1) = 1.d0
            inc%Saturation(2) = 0.d0
         end if

      else if (ic == 2) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression < Psat) then
            inc%ic = 3
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(1) = 0.d0
            inc%Saturation(2) = 1.d0
         end if

      else if (ic == 3) then

         call FluidThermodynamics_Tsat(inc%Pression, Tsat, dTsatdP)
         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         inc%Temperature = Tsat
         inc%Pression = Psat

         if (inc%Saturation(1) < 0.d0) then
            inc%ic = 2
            inc%Saturation(1) = 0.d0
            inc%Saturation(2) = 1.d0
         else if (inc%Saturation(2) < 0.d0) then
            inc%ic = 1
            inc%Saturation(1) = 1.d0
            inc%Saturation(2) = 0.d0
         end if

      else
         print *, "Error in Flash: no such context"
      end if

   end subroutine DefFlash_Flash_cv

end module DefFlash
