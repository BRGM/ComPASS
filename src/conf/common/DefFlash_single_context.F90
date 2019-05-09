!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!
    
! No flash, typically when a single phase is present
    
module DefFlash

   use IncCVReservoir, only: TYPE_IncCVReservoir
   use DefModel, only: IndThermique

   implicit none

   integer, parameter :: SINGLE_CONTEXT = 1
   
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
      integer, intent(in) :: rocktype(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      inc%ic = SINGLE_CONTEXT

   end subroutine DefFlash_Flash_cv

end module DefFlash
