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
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      GAS_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP
   use Thermodynamics, only: f_Fugacity, f_Fugacity_coefficient

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   subroutine DiphasicFlash_fugacity_coeffficients(icp, inc, fg, fl)
      integer(c_int), intent(in) :: icp
      type(Type_IncCVReservoir), intent(inout) :: inc
      real(c_double), intent(out) :: fg, fl

      real(c_double) :: dPf, dTf, dCf(NbComp) ! dummy values

      call f_Fugacity_coefficient(icp, GAS_PHASE, inc%phase_pressure(GAS_PHASE), &
                                  inc%Temperature, inc%Comp(:, GAS_PHASE), fg, dPf, dTf, dCf)
      call f_Fugacity_coefficient(icp, LIQUID_PHASE, inc%phase_pressure(LIQUID_PHASE), &
                                  inc%Temperature, inc%Comp(:, LIQUID_PHASE), fl, dPf, dTf, dCf)

   end subroutine DiphasicFlash_fugacity_coeffficients

   subroutine DiphasicFlash_liquid_fugacities(inc, fa, fw)
      type(Type_IncCVReservoir), intent(inout) :: inc
      real(c_double), intent(out) :: fa, fw

      real(c_double) :: dPf, dTf, dCf(NbComp) ! dummy values

      call f_Fugacity(AIR_COMP, LIQUID_PHASE, inc%phase_pressure(LIQUID_PHASE), &
                      inc%Temperature, inc%Comp(:, LIQUID_PHASE), fa, dPf, dTf, dCf)
      call f_Fugacity(WATER_COMP, LIQUID_PHASE, inc%phase_pressure(LIQUID_PHASE), &
                      inc%Temperature, inc%Comp(:, LIQUID_PHASE), fw, dPf, dTf, dCf)

   end subroutine DiphasicFlash_liquid_fugacities

   subroutine DiphasicFlash_liquid_to_diphasic(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      real(c_double) :: fa, fw

      ! compute liquid fugacity of both components
      call DiphasicFlash_liquid_fugacities(inc, fa, fw)

      ! WARNING: the following relies on the ideal gas assumption
      !          and assumes fa = PaCag and fw = PgCwg
      ! WARNING: don't divide inequality by Pg (migth be negative during Newton iteration)
      if (fa + fw > inc%phase_pressure(GAS_PHASE)) then
         ! FIXME: is it ok to compare to Pg when gas is not ideal? (cf. issue #159)
         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
      endif

   end subroutine DiphasicFlash_liquid_to_diphasic

   subroutine DiphasicFlash_gas_to_diphasic(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      real(c_double) :: fg, fl
      real(c_double) :: Cla, Clw

      call DiphasicFlash_fugacity_coeffficients(AIR_COMP, inc, fg, fl)
      Cla = (fg/fl)*inc%Comp(AIR_COMP, GAS_PHASE)
      call DiphasicFlash_fugacity_coeffficients(WATER_COMP, inc, fg, fl)
      Clw = (fg/fl)*inc%Comp(WATER_COMP, GAS_PHASE)
      if (Cla + Clw > 1.d0) then ! Liquid appears
         inc%ic = DIPHASIC_CONTEXT
         inc%Saturation(GAS_PHASE) = 1.d0
         inc%Saturation(LIQUID_PHASE) = 0.d0
         inc%Comp(AIR_COMP, LIQUID_PHASE) = Cla
      endif

   end subroutine DiphasicFlash_gas_to_diphasic

   pure subroutine DiphasicFlash_diphasic_switches(inc)
      type(Type_IncCVReservoir), intent(inout) :: inc

      if (inc%Saturation(GAS_PHASE) < 0.d0) then ! gas vanishes
         inc%ic = LIQUID_CONTEXT
         inc%Saturation(GAS_PHASE) = 0.d0
         inc%Saturation(LIQUID_PHASE) = 1.d0
      else if (inc%Saturation(LIQUID_PHASE) < 0.d0) then ! liquid vanishes
         inc%ic = GAS_CONTEXT
         inc%Saturation(GAS_PHASE) = 1.d0
         inc%Saturation(LIQUID_PHASE) = 0.d0
      endif

   end subroutine DiphasicFlash_diphasic_switches

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

         call DiphasicFlash_liquid_to_diphasic(inc)

      elseif (context == DIPHASIC_CONTEXT) then

         call DiphasicFlash_diphasic_switches(inc)

      elseif (context == GAS_CONTEXT) then

         call DiphasicFlash_gas_to_diphasic(inc)

      endif

   end subroutine DiphasicFlash_Flash_cv

   !< enforce C in [0,1] and sum equal to 1
   pure subroutine DiphasicFlash_enforce_consistent_molar_fractions(inc) &
      bind(C, name="DiphasicFlash_enforce_consistent_molar_fractions")
      type(Type_IncCVReservoir), intent(inout) :: inc

      integer :: alpha
      real(c_double) :: Ca

      do alpha = 1, NbPhase ! NbPhase = 2
         Ca = min(max(inc%Comp(AIR_COMP, alpha), 0.d0), 1.d0)
         inc%Comp(AIR_COMP, alpha) = Ca
         inc%Comp(WATER_COMP, alpha) = 1.d0 - Ca
      enddo

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

      context = inc%ic

      if (context == GAS_FF_NO_LIQ_OUTFLOW_CONTEXT) then

         call DiphasicFlash_gas_to_diphasic(inc)
         if (inc%ic == DIPHASIC_CONTEXT) inc%ic = DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT

      elseif (context == DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT) then

         ! CHECKME: we assume no entry pressure in the atmosphere
         call DiphasicFlash_diphasic_switches(inc)
         if (inc%ic == LIQUID_CONTEXT) then
            inc%ic = DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
         elseif (inc%ic == GAS_CONTEXT) then
            inc%ic = GAS_FF_NO_LIQ_OUTFLOW_CONTEXT
         endif
         ! The following always holds even at liquid outflow apparition
         inc%FreeFlow_flowrate(LIQUID_PHASE) = 0.d0

      elseif (context == DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT) then

         if (inc%FreeFlow_flowrate(LIQUID_PHASE) < 0.d0) then ! liquid outflow vanishes
            inc%ic = DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT
            inc%FreeFlow_flowrate(LIQUID_PHASE) = 0.d0
         endif

      else

         call DiphasicFlash_Flash_cv(inc)

      endif

      call DiphasicFlash_enforce_consistent_molar_fractions(inc)

   end subroutine DefFlash_Flash_cv

end module DefFlash
