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

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume ?????
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

      type(TYPE_IncCVReservoir), intent(inout) :: inc
      INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      integer :: i, iph, j, icp, m, mph, ic
      double precision :: DensiteMolaire(NbComp), acc1, acc2, &
         dPf, dTf, dCf(NbComp), dSf(NbPhase)

      double precision :: Tsat, dTsatdP, Psat, dPsatdT

      if (ic /= 1) then
         print *, "Error in Flash: no such context"
      end if

      inc%ic = 1

   end subroutine DefFlash_Flash_cv

end module DefFlash
