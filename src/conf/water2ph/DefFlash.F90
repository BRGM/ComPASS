!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp, context switch

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlash

   use Thermodynamics
   use IncCVReservoir
   use VAGFrac ! for rocktypes

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

      integer :: i, iph, j, icp, m, mph, ic, errcode, Ierr
      double precision :: DensiteMolaire(NbComp), acc1, acc2, &
         dPf, dTf, dCf(NbComp), dSf(NbPhase)

      double precision :: Tsat, dTsatdP, Psat, dPsatdT

      ic = inc%ic

      if (ic == GAS_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression > Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(1) = 1.d0
            inc%Saturation(2) = 0.d0
         end if

      else if (ic == LIQUID_CONTEXT) then

         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         if (inc%Pression < Psat) then
            inc%ic = DIPHASIC_CONTEXT
            inc%Pression = Psat
            ! inc%Temperature is the saturation temperature (by construction)
            inc%Saturation(1) = 0.d0
            inc%Saturation(2) = 1.d0
         end if

      else if (ic == DIPHASIC_CONTEXT) then

#ifndef _THERMIQUE_
         write(*,*) 'ERROR: Diphasic context is meaningless without energy transfer!'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
#endif
          
         call FluidThermodynamics_Tsat(inc%Pression, Tsat, dTsatdP)
         call FluidThermodynamics_Psat(inc%Temperature, Psat, dPsatdT)

         inc%Temperature = Tsat
         inc%Pression = Psat

         if (inc%Saturation(1) < 0.d0) then
            inc%ic = LIQUID_CONTEXT
            inc%Saturation(1) = 0.d0
            inc%Saturation(2) = 1.d0
         else if (inc%Saturation(2) < 0.d0) then
            inc%ic = GAS_CONTEXT
            inc%Saturation(1) = 1.d0
            inc%Saturation(2) = 0.d0
         end if

      else
         print *, "Error in Flash: no such context"
      end if

   end subroutine DefFlash_Flash_cv

end module DefFlash
