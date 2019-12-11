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

#ifndef _THERMIQUE_
   use CommonMPI, only: CommonMPI_abort
#endif
   use IncCVReservoir, only: Type_IncCVReservoir
   use DefModel, only: &
      IndThermique, NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   ! !> \brief Determine the phases
   ! !! which are actualy present.
   ! !!
   ! !! Applied to IncNode, IncFrac and IncCell.
   ! !! \param[in]      porovol   porous Volume ?????
   ! !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   ! subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

   !    type(TYPE_IncCVReservoir), intent(inout) :: inc
   !    INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
   !    double precision, intent(in) :: porovol ! porovol

   !    integer :: context, iph
   !    double precision :: eps, pure_phase_molar_fraction(2,2)

   !    eps = 1.d-18 ! 0 reference
   !    pure_phase_molar_fraction(:,:) = 0.d0
   !    pure_phase_molar_fraction(AIR_COMP, GAS_PHASE) = 1.d0
   !    pure_phase_molar_fraction(WATER_COMP, LIQUID_PHASE) = 1.d0
   !    context = inc%ic

   !    ! In this physic we are considering immiscible gas-air // liquid-water
   !    ! so the flash is done over the component present 
   !    ! only in the missing phase (Ctilde)
   !    if (context == GAS_CONTEXT) then

   !       if(dabs(inc%AccVol(WATER_COMP))<eps) &
   !          inc%AccVol(WATER_COMP) = 0.d0
   !       ! inc%%AccVol contains the increment of Ctilde
   !       if (inc%AccVol(WATER_COMP) > eps) then
   !          print*,"liquid appearance",inc%AccVol(WATER_COMP)
   !          inc%ic = DIPHASIC_CONTEXT
   !          ! inc%Temperature is the saturation temperature (by construction)
   !          inc%Saturation(GAS_PHASE) = 1.d0
   !          inc%Saturation(LIQUID_PHASE) = 0.d0
   !          inc%Comp = pure_phase_molar_fraction
   !       end if

   !    else if (context == LIQUID_CONTEXT) then

   !       if(dabs(inc%AccVol(AIR_COMP))<eps) &
   !          inc%AccVol(AIR_COMP) = 0.d0
   !       ! inc%%AccVol contains the increment of Ctilde
   !       if (inc%AccVol(AIR_COMP) > eps) then
   !          print*,"gas appearance",inc%AccVol(AIR_COMP)
   !          inc%ic = DIPHASIC_CONTEXT
   !          ! inc%Temperature is the saturation temperature (by construction)
   !          inc%Saturation(GAS_PHASE) = 0.d0
   !          inc%Saturation(LIQUID_PHASE) = 1.d0
   !          inc%Comp = pure_phase_molar_fraction
   !       end if

   !    else if (context == DIPHASIC_CONTEXT) then

   !       do iph = 1, NbPhase
   !          if(dabs(inc%Saturation(iph))<eps) &
   !             inc%Saturation(iph) = 0.d0
   !       enddo
   !       if (inc%Saturation(GAS_PHASE) < -eps) then
   !          print*,"gas disappearance",inc%Saturation(GAS_PHASE)
   !          inc%ic = LIQUID_CONTEXT
   !          inc%Saturation(GAS_PHASE) = 0.d0
   !          inc%Saturation(LIQUID_PHASE) = 1.d0
   !          inc%Comp = pure_phase_molar_fraction
   !          inc%AccVol(AIR_COMP) = 0.d0

   !       else if (inc%Saturation(LIQUID_PHASE) < -eps) then
   !          print*,"liquid disappearance",inc%Saturation(LIQUID_PHASE)
   !          inc%ic = GAS_CONTEXT
   !          inc%Saturation(GAS_PHASE) = 1.d0
   !          inc%Saturation(LIQUID_PHASE) = 0.d0
   !          inc%Comp = pure_phase_molar_fraction
   !          inc%AccVol(WATER_COMP) = 0.d0

   !       end if

   !    else
   !       print *, "Error in Flash: no such context"
   !    end if

   ! end subroutine DefFlash_Flash_cv

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

   type(TYPE_IncCVReservoir), intent(inout) :: inc
   INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
   double precision, intent(in) :: porovol ! porovol

   integer :: context, iph
   double precision :: eps, pure_phase_molar_fraction(2,2)

   eps = 1.d-18 ! 0 reference
   pure_phase_molar_fraction(:,:) = 0.d0
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

   else
      print *, "Error in Flash: no such context"
   end if

end subroutine DefFlash_Flash_cv

end module DefFlash
