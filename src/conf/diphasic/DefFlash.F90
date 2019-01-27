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

   use IncCVReservoir
   use Physics
   use VAGFrac ! to have rocktypes

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
   subroutine DefFlash_Flash_cv(inc, rt, porovol)

      type(Type_IncCVReservoir), intent(inout) :: inc
      integer, intent(in) :: rt(IndThermique + 1)
      double precision, intent(in) :: porovol ! porovol

      integer :: iph, icp, ic
      double precision :: T, f(NbPhase)
      double precision :: S(NbPhase), Pc, DSPc(NbPhase)
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: Cag, Cal, Cel
      double precision :: PgCag, PgCeg
      double precision :: Pref,P(NbPhase)

      ic = inc%ic
      T = inc%Temperature
      S = inc%Saturation
      Pref = inc%Pression
      do iph = 1,NbPhase
        call f_PressionCapillaire(rt, iph, S(iph), Pc, DSPc)
        P(iph) = Pref + Pc
      enddo

      if (ic == LIQUID_CONTEXT) then
         ! air liq fugacity
         iph = PHASE_WATER
         call f_Fugacity(rt,iph,1,P(iph),T,inc%Comp(:,iph),S(iph),f(iph),DPf,DTf,DCf,DSf)
         PgCag = inc%Comp(1, iph)*f(iph)

         ! water liq fugacity
         call f_Fugacity(rt,iph,2,P(iph),T,inc%Comp(:,iph),S(iph),f(iph),DPf,DTf,DCf,DSf)
         PgCeg = inc%Comp(2, iph)*f(iph)

         ! don't divide inequality by Pg (migth be negative during Newton iteration)
         if (PgCag + PgCeg > P(PHASE_GAS)) then

            ! write(*,*)' apparition gas ', P(PHASE_GAS), T

            inc%ic = DIPHASIC_CONTEXT
            inc%Saturation(PHASE_GAS) = 0
            inc%Saturation(PHASE_WATER) = 1
            inc%Comp(1,PHASE_GAS) = MIN(MAX(inc%Comp(1,PHASE_GAS), 0.d0), 1.d0)
            inc%Comp(2,PHASE_GAS) = 1.d0 - inc%Comp(1,PHASE_GAS)

         endif

      elseif (ic == DIPHASIC_CONTEXT) then

         if (S(PHASE_GAS) < 0.d0) then

            ! write(*,*)' disapparition gas ', P(PHASE_GAS), T

            inc%ic = LIQUID_CONTEXT
            inc%Saturation(PHASE_GAS) = 0.d0
            inc%Saturation(PHASE_WATER) = 1.d0

         elseif (S(PHASE_WATER) < 0.d0) then

            ! write (*, *) ' disapparition liquid ', P(PHASE_GAS), T

            inc%ic = GAS_CONTEXT
            inc%Saturation(PHASE_GAS) = 1.d0
            inc%Saturation(PHASE_WATER) = 0.d0

         endif

         ! force comp to be in [0,1] and sum equal to 1
         do iph = 1, NbPhase
            inc%Comp(1,iph) = MIN(MAX(inc%Comp(1,iph), 0.d0), 1.d0)
            inc%Comp(2,iph) = 1.d0 - inc%Comp(1,iph)
         enddo

      elseif (ic == GAS_CONTEXT) then
         ! air
         do iph = 1, NbPhase
            call f_Fugacity(rt,iph,1,P(iph),T,inc%Comp(:,iph),S(iph),f(iph),DPf,DTf,DCf,DSf)
         enddo
         Cal = inc%Comp(1,PHASE_GAS)*f(PHASE_GAS)/f(PHASE_WATER)
         ! water
         do iph = 1, NbPhase
            call f_Fugacity(rt,iph,2,P(iph),T,inc%Comp(:,iph),S(iph),f(iph),DPf,DTf,DCf,DSf)
         enddo
         Cel = inc%Comp(2,PHASE_GAS)*f(PHASE_GAS)/f(PHASE_WATER)

         if(Cal + Cel > 1.d0) then

            ! write(*,*)' apparition liquid ', P(PHASE_GAS), T

            inc%ic = DIPHASIC_CONTEXT

            inc%Saturation(PHASE_GAS) = 1.d0
            inc%Saturation(PHASE_WATER) = 0.d0
            Cal = MIN(MAX(Cal, 0.d0), 1.d0)
            inc%Comp(1,PHASE_WATER) = Cal
            inc%Comp(2,PHASE_WATER) = 1.d0 - Cal
         endif

      else
         print *, "Error in Flash: unknown context"
         stop
      endif

   end subroutine DefFlash_Flash_cv

end module DefFlash
