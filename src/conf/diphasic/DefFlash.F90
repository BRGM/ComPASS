!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefFlash

   use iso_c_binding, only: c_int, c_double
   use CommonMPI, only: CommonMPI_abort
   use IncCVReservoirTypes, only: Type_IncCVReservoir
   use DefModel, only: &
      NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      GAS_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP
   ! f_Fugacity_coefficient called in ../common/DiphasicFlash.F90
   use Thermodynamics, only: f_Fugacity, f_Fugacity_coefficient

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

#include "../common/DiphasicFlash.F90"

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
