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

   use IncCVReservoir, only: Type_IncCVReservoir
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
   subroutine DefFlash_Flash_cv(inc, rocktype, porovol)

      type(TYPE_IncCVReservoir), intent(inout) :: inc
      INTEGER, INTENT(IN) :: rocktype(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      integer :: context

      double precision :: Tsat, dTsatdP, Psat, dPsatdT

      context = inc%ic

      if (locked_context(context)) return

      if (context == GAS_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression > Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0
         end if

      else if (context == LIQUID_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression < Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0
         end if

      else if (context == DIPHASIC_CONTEXT) then

#ifndef _THERMIQUE_
         write (*, *) 'ERROR: Diphasic context is meaningless without energy transfer!'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
#endif

         call FluidThermodynamics_Tsat(inc%Pression, Tsat, dTsatdP)
         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         inc%Temperature = Tsat
         inc%Pression = Psat

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
