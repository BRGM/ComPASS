!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!
! Model: 2 phase 2 comp, thermal well
! water can be present only in the liquid phase

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).

module DefFlash

   use iso_c_binding, only: c_int, c_double
   use CommonMPI, only: CommonMPI_abort
   use IncCVReservoirTypes, only: Type_IncCVReservoir
   use DefModel, only: &
      NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      LIQUID2DIPHASIC_CONTEXT, DIPHASIC2LIQUID_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, CO2_COMP, WATER_COMP
   use Thermodynamics, only: f_Fugacity, f_Fugacity_coeff_CO2_gas

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   subroutine DiphasicFlash_liquid_to_liq2diph(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      real(c_double) :: fa, phi
      real(c_double) :: dPf, dTf, dCf(NbComp) ! dummy values

      ! compute liquid fugacity of CO2 in liquid
      call f_Fugacity(CO2_COMP, LIQUID_PHASE, inc%phase_pressure(LIQUID_PHASE), &
                      inc%Temperature, inc%Comp(:, LIQUID_PHASE), fa, dPf, dTf, dCf)
      call f_Fugacity_coeff_CO2_gas(inc%phase_pressure(GAS_PHASE), inc%Temperature, inc%Comp(:, GAS_PHASE), phi, dPf, dTf, dCf)
      if (fa > phi*inc%phase_pressure(GAS_PHASE)) then
         ! gas should appear, change into intermediate liq2diph
         ! write(*, *) "liquid to liq2diph"
         inc%ic = LIQUID2DIPHASIC_CONTEXT
      endif

   end subroutine DiphasicFlash_liquid_to_liq2diph

   subroutine DiphasicFlash_liq2diph_switches(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      real(c_double) :: fa, phi
      real(c_double) :: dPf, dTf, dCf(NbComp) ! dummy values

      ! compute liquid fugacity of CO2 in liquid
      call f_Fugacity(CO2_COMP, LIQUID_PHASE, inc%phase_pressure(LIQUID_PHASE), &
                      inc%Temperature, inc%Comp(:, LIQUID_PHASE), fa, dPf, dTf, dCf)
      call f_Fugacity_coeff_CO2_gas(inc%phase_pressure(GAS_PHASE), inc%Temperature, inc%Comp(:, GAS_PHASE), phi, dPf, dTf, dCf)
      if (fa > phi*inc%phase_pressure(GAS_PHASE)) then
         ! write(*, *) "liq2diph to diphasic"
         ! FIXME: is it ok to compare to Pg when gas is not ideal? (cf. issue #159)
         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
      else
         ! write(*, *) "leave liq2diph (go back to liquid) "
         ! leave this intermediate context, go back to liquid_ctx
         inc%ic = LIQUID_CONTEXT
      endif

   end subroutine DiphasicFlash_liq2diph_switches

   pure subroutine DiphasicFlash_diphasic_switches(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      if (inc%Saturation(GAS_PHASE) < 0.d0) then ! gas vanishes
         ! write(*, *) "diphasic to diph2liq"
         inc%ic = DIPHASIC2LIQUID_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
      else if (inc%Saturation(LIQUID_PHASE) < 0.d0) then ! liquid vanishes
         ! write(*, *) "diphasic to gas"
         inc%ic = GAS_CONTEXT
         inc%Saturation(GAS_PHASE) = 1.d0
         inc%Saturation(LIQUID_PHASE) = 0.d0
      endif

   end subroutine DiphasicFlash_diphasic_switches

   pure subroutine DiphasicFlash_diph2liq_switches(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      if (inc%Saturation(GAS_PHASE) < 0.d0) then ! gas vanishes
         ! write(*, *) "diph2liq to liquid"
         inc%ic = LIQUID_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
      else
         ! write(*, *) "leave diph2liq (go back to diphasic) "
         ! leave this intermediate context, go back to diphasic_ctx
         inc%ic = DIPHASIC_CONTEXT
      endif

   end subroutine DiphasicFlash_diph2liq_switches

   pure subroutine DiphasicFlash_gas_to_diphasic(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc
      real(c_double), parameter :: eps = 1.0d-20

      ! Here we test the existence of liquid with the existence of water
      ! inc%AccVol contains ni when the component is not present in the context
      ! There is no water in the gas phase
      ! eps=0 creates oscillations with values such as ni = 1e-90...
      if (inc%AccVol(WATER_COMP) > eps) then
         ! write(*, *) "gas to diphasic"
         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 1.d0
         inc%Saturation(LIQUID_PHASE) = 0.d0
      endif

   end subroutine DiphasicFlash_gas_to_diphasic

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume ?????
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DiphasicFlash_Flash_cv(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      integer(c_int) :: context

      context = inc%ic

      if (context == LIQUID_CONTEXT) then

         call DiphasicFlash_liquid_to_liq2diph(inc)

      elseif (context == LIQUID2DIPHASIC_CONTEXT) then

         call DiphasicFlash_liq2diph_switches(inc)

      elseif (context == DIPHASIC_CONTEXT) then

         call DiphasicFlash_diphasic_switches(inc)

      elseif (context == DIPHASIC2LIQUID_CONTEXT) then

         call DiphasicFlash_diph2liq_switches(inc)

      elseif (context == GAS_CONTEXT) then

         call DiphasicFlash_gas_to_diphasic(inc)

      endif

   end subroutine DiphasicFlash_Flash_cv

   !< enforce gas molar fractions and Cl in [0,1] and sum equal to 1
   pure subroutine DiphasicFlash_enforce_consistent_molar_fractions(inc) &
      bind(C, name="DiphasicFlash_enforce_consistent_molar_fractions")
      type(Type_IncCVReservoir), intent(inout) :: inc

      integer :: alpha
      real(c_double) :: Cal

      ! there is no gas water
      inc%Comp(CO2_COMP, GAS_PHASE) = 1.d0
      inc%Comp(WATER_COMP, GAS_PHASE) = 0.d0

      ! enforce Cl in [0,1] and sum equal to 1
      Cal = min(max(inc%Comp(CO2_COMP, LIQUID_PHASE), 0.d0), 1.d0)
      inc%Comp(CO2_COMP, LIQUID_PHASE) = Cal
      inc%Comp(WATER_COMP, LIQUID_PHASE) = 1.d0 - Cal

   end subroutine DiphasicFlash_enforce_consistent_molar_fractions

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      integer :: iph, context, i
      double precision :: T, f(NbPhase)
      double precision :: S(NbPhase)
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: Cal, Cwl
      double precision :: PgCag, PgCwg
      double precision :: Pref, Pg

#ifndef NDEBUG
      if (NbPhase /= 2) then
         call CommonMPI_abort("Wrong number of phases.")
      endif
#endif

      call DiphasicFlash_Flash_cv(inc)

      call DiphasicFlash_enforce_consistent_molar_fractions(inc)

   end subroutine DefFlash_Flash_cv

end module DefFlash
