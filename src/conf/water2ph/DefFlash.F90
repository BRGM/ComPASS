!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp, context switch

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use iso_c_binding, only: c_double
   use IncCVReservoirTypes, only: Type_IncCVReservoir
   use DefModel, only: &
      locked_context, &
      IndThermique, NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, GAS_PHASE, LIQUID_PHASE
   use Thermodynamics, only: FluidThermodynamics_Psat, FluidThermodynamics_Tsat

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

      integer :: context
      double precision :: Tsat, Psat, unused

      context = inc%ic

      if (locked_context(context)) return

      if (context == GAS_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, unused)

         if (inc%phase_pressure(GAS_PHASE) > Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat ! Pression = Pref is Pg in this physical model
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0
         end if

      else if (context == LIQUID_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, unused)

         if (inc%phase_pressure(GAS_PHASE) < Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat! Pression = Pref is Pg in this physical model
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0
         end if

      else if (context == DIPHASIC_CONTEXT) then

#ifndef _THERMIQUE_
         write (*, *) 'ERROR: Diphasic context is meaningless without energy transfer!'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
#endif

         ! inc%Pression = Pref = Pg in this physic
         call FluidThermodynamics_Tsat(inc%Pression, Tsat, unused)
         inc%Temperature = Tsat

         if (inc%Saturation(GAS_PHASE) < 0.d0) then
            inc%ic = LIQUID_CONTEXT
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0
         else if (inc%Saturation(LIQUID_PHASE) < 0.d0) then
            inc%ic = GAS_CONTEXT
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0
         end if

      else
         print *, "Error in Flash: no such context"
      end if

   end subroutine DefFlash_Flash_cv

end module DefFlash
