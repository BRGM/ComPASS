!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp, thermal well

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use iso_c_binding, only: c_int, c_double
   use CommonMPI, only: CommonMPI_abort
   use IncCVReservoirTypes, only: Type_IncCVReservoir
   use DefModel, only: &
      NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP
   use Thermodynamics, only: f_Fugacity

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   include "../common/DiphasicFlash.F90"

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume ?????
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

#ifndef NDEBUG
      if (NbPhase /= 2) then
         call CommonMPI_abort("Wrong number of phases.")
      endif
#endif

      call DiphasicFlash_Flash_cv(inc)
      call DiphasicFlash_enforce_consistent_molar_fractions(inc)

   end subroutine DefFlash_Flash_cv

end module DefFlash
