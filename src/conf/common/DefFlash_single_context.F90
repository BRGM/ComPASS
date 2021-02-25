!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! No flash, typically when a single phase is present

module DefFlash

   use iso_c_binding, only: c_double
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir
   use DefModel, only: IndThermique, NbPhase

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
   subroutine DefFlash_Flash_cv(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      inc%ic = ComPASS_SINGLE_CONTEXT

   end subroutine DefFlash_Flash_cv

end module DefFlash
