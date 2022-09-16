!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>   This file contains an example on how to deal with change of phase when the physic contains
!!   components which are not contained in every phase (MCP not full). It has been commented
!!   because in this physic each component is present in only one phase. Thus it is not necessary
!!   to deal with change of phase, the simulation can be run only with "diphasic" context because
!!   the equilibrium is never set between the components.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Model: 2 phase 2 comp, context switch

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use iso_c_binding, only: c_double
   use IncCVReservoirTypes, only: Type_IncCVReservoir
   use DefModel, only: &
      IndThermique, NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP
#ifndef NDEBUG
   use CommonMPI, only: CommonMPI_abort
#endif

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      integer :: context
      double precision :: eps, pure_phase_molar_fraction(2, 2)

      eps = 1.d-18 ! 0 reference
      pure_phase_molar_fraction(:, :) = 0.d0
      pure_phase_molar_fraction(AIR_COMP, GAS_PHASE) = 1.d0
      pure_phase_molar_fraction(WATER_COMP, LIQUID_PHASE) = 1.d0
      context = inc%ic

      ! In this physic we are considering immiscible gas-air // liquid-water
      ! so the flash is not necessary (put always diphasic context)
      if (context == GAS_CONTEXT) then ! if bad initialization, move to diphasic

         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 1.d0
         inc%Saturation(LIQUID_PHASE) = 0.d0
         inc%Comp = pure_phase_molar_fraction

      else if (context == LIQUID_CONTEXT) then ! if bad initialization, move to diphasic

         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
         inc%Comp = pure_phase_molar_fraction

      else if (context == DIPHASIC_CONTEXT) then

         inc%Saturation(GAS_PHASE) = MIN(MAX(inc%Saturation(GAS_PHASE), 0.d0), 1.d0)
         inc%Saturation(LIQUID_PHASE) = MIN(MAX(inc%Saturation(LIQUID_PHASE), 0.d0), 1.d0)

#ifndef NDEBUG
         ! needs that MCP = 0 implies C = 0 : refer to issue #552
         if (dabs(inc%Comp(WATER_COMP, GAS_PHASE)) .gt. 1e-7) &
            call CommonMPI_abort("C(WATER_COMP, GAS_PHASE) != 0")
         if (dabs(inc%Comp(AIR_COMP, LIQUID_PHASE)) .gt. 1e-7) &
            call CommonMPI_abort("C(AIR_COMP, LIQUID_PHASE) != 0")
#endif
      else
         print *, "Error in Flash: no such context"
      end if

   end subroutine DefFlash_Flash_cv

end module DefFlash
