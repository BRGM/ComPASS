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
   use IncCVReservoir
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

   subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

      type(TYPE_IncCVReservoir), intent(inout) :: inc
      INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      inc%ic = 1

   end subroutine DefFlash_Flash_cv

end module DefFlash
