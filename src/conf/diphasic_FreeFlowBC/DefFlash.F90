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

   use CommonMPI, only: CommonMPI_abort
   use IncCVReservoir, only: Type_IncCVReservoir
   use DefModel, only: &
      IndThermique, NbPhase, NbComp, &
      DIPHASIC_CONTEXT, LIQUID_CONTEXT, GAS_CONTEXT, &
      GAS_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT, DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT, &
      GAS_PHASE, LIQUID_PHASE, AIR_COMP, WATER_COMP
   use Thermodynamics, only: f_Fugacity, f_PressionCapillaire

   implicit none

   public :: &
      DefFlash_Flash_cv

contains

   !> \brief Determine the phases
   !! which are actualy present.
   !!
   !! Applied to IncNode, IncFrac and IncCell.
   !! \param[in]      porovol   porous Volume
   !! \param[in]      rt        rocktype
   !! \param[inout]   inc       Unknown (IncNode, IncFrac or IncCell)
   subroutine DefFlash_Flash_cv(inc, rt, porovol)

      type(Type_IncCVReservoir), intent(inout) :: inc
      integer, intent(in) :: rt
      double precision, intent(in) :: porovol ! porous volume

      integer :: iph, ic, i
      double precision :: T, f(NbPhase)
      double precision :: S(NbPhase), Pc, DSPc(NbPhase)
      double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: Cal, Cwl
      double precision :: PgCag, PgCwg
      double precision :: Pref, Pg

      ic = inc%ic
      T = inc%Temperature
      S = inc%Saturation
      Pref = inc%Pression
      ! compute Pg
      call f_PressionCapillaire(rt, GAS_PHASE, S, Pc, DSPc)
      Pg = Pref + Pc

      ! RESRVOIR DOF
      if (ic == LIQUID_CONTEXT) then

         ! air liq fugacity
         iph = LIQUID_PHASE
         ! f_Fugacity is computed using Pref
         call f_Fugacity(rt, iph, AIR_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         PgCag = inc%Comp(AIR_COMP, iph)*f(iph)

         ! water liq fugacity, f_Fugacity is computed using Pref
         call f_Fugacity(rt, iph, WATER_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         PgCwg = inc%Comp(WATER_COMP, iph)*f(iph)

         ! don't divide inequality by Pg (migth be negative during Newton iteration)
         if (PgCag + PgCwg > Pg) then
            ! write(*,*)' appearance gas ', Pg, T
            inc%ic = DIPHASIC_CONTEXT
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0
            ! force comp to be in [0,1]
            do i = 1, NbComp
               inc%Comp(i, GAS_PHASE) = min(max(inc%Comp(i, GAS_PHASE), 0.d0), 1.d0)
            enddo

         endif

         ! force comp to be in [0,1]
         do i = 1, NbComp
            inc%Comp(i, LIQUID_PHASE) = min(max(inc%Comp(i, LIQUID_PHASE), 0.d0), 1.d0)
         enddo

      elseif (ic == DIPHASIC_CONTEXT) then

         if (S(GAS_PHASE) < 0.d0) then
            ! write(*,*)' disappearance gas ', Pg, T
            inc%ic = LIQUID_CONTEXT
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0

         elseif (S(LIQUID_PHASE) < 0.d0) then
            ! write (*, *) ' disappearance liquid ', Pg, T
            inc%ic = GAS_CONTEXT
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0

         endif

         ! force comp to be in [0,1]
         do iph = 1, NbPhase
            do i = 1, NbComp
               inc%Comp(i, iph) = min(max(inc%Comp(i, iph), 0.d0), 1.d0)
            enddo
         enddo

      elseif (ic == GAS_CONTEXT) then

         ! air
         do iph = 1, NbPhase
            call f_Fugacity(rt, iph, AIR_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         enddo
         Cal = inc%Comp(AIR_COMP, GAS_PHASE)*f(GAS_PHASE)/f(LIQUID_PHASE)
         ! water
         do iph = 1, NbPhase
            call f_Fugacity(rt, iph, WATER_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         enddo
         Cwl = inc%Comp(WATER_COMP, GAS_PHASE)*f(GAS_PHASE)/f(LIQUID_PHASE)

         if (Cal + Cwl > 1.d0) then
            ! write(*,*)' appearance liquid '
            inc%ic = DIPHASIC_CONTEXT
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0
            ! force comp to be in [0,1]
            do i = 1, NbComp
               inc%Comp(i, LIQUID_PHASE) = min(max(inc%Comp(i, LIQUID_PHASE), 0.d0), 1.d0)
            enddo
         endif

         ! force comp to be in [0,1]
         do i = 1, NbComp
            inc%Comp(i, GAS_PHASE) = min(max(inc%Comp(i, GAS_PHASE), 0.d0), 1.d0)
         enddo

         ! FREEFLOW BC DOF
      elseif (ic == GAS_FF_NO_LIQ_OUTFLOW_CONTEXT) then

         ! air
         do iph = 1, NbPhase
            call f_Fugacity(rt, iph, AIR_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         enddo
         Cal = inc%Comp(AIR_COMP, GAS_PHASE)*f(GAS_PHASE)/f(LIQUID_PHASE)
         ! water
         do iph = 1, NbPhase
            call f_Fugacity(rt, iph, WATER_COMP, Pref, T, inc%Comp(:, iph), S, f(iph), DPf, DTf, DCf, DSf)
         enddo
         Cwl = inc%Comp(WATER_COMP, GAS_PHASE)*f(GAS_PHASE)/f(LIQUID_PHASE)

         if (Cal + Cwl > 1.d0) then
            inc%ic = DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0
            ! force comp to be in [0,1]
            do i = 1, NbComp
               inc%Comp(i, LIQUID_PHASE) = min(max(inc%Comp(i, LIQUID_PHASE), 0.d0), 1.d0)
            enddo
         endif
         ! force comp to be in [0,1]
         do i = 1, NbComp
            inc%Comp(i, GAS_PHASE) = min(max(inc%Comp(i, GAS_PHASE), 0.d0), 1.d0)
         enddo

      elseif (ic == DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT) then

         if (S(GAS_PHASE) < 0.d0) then ! because there is no entry pressure in the atmosphere
            ! write(*,*)' appearance liquid outflow in Freeflow BC', Pg, T
            inc%ic = DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
            ! important to set the following saturations
            ! because they are not unknowns in this context
            inc%Saturation(GAS_PHASE) = 0.d0
            inc%Saturation(LIQUID_PHASE) = 1.d0
            inc%FreeFlow_flowrate(LIQUID_PHASE) = 0.d0

         elseif (S(LIQUID_PHASE) < 0.d0) then
            ! write (*, *) ' disappearance liquid phase in Freeflow BC', Pg, T
            inc%ic = GAS_FF_NO_LIQ_OUTFLOW_CONTEXT
            inc%Saturation(GAS_PHASE) = 1.d0
            inc%Saturation(LIQUID_PHASE) = 0.d0

         endif

         ! force comp to be in [0,1]
         do iph = 1, NbPhase
            do i = 1, NbComp
               inc%Comp(i, iph) = min(max(inc%Comp(i, iph), 0.d0), 1.d0)
            enddo
         enddo

      elseif (ic == DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT) then

         if (inc%FreeFlow_flowrate(LIQUID_PHASE) < 0.d0) then
            ! write(*,*)' disappearance liquid outflow in Freeflow BC'
            inc%ic = DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT
            inc%FreeFlow_flowrate(LIQUID_PHASE) = 0.d0
         endif

         ! force comp to be in [0,1]
         do iph = 1, NbPhase
            do i = 1, NbComp
               inc%Comp(i, iph) = min(max(inc%Comp(i, iph), 0.d0), 1.d0)
            enddo
         enddo

      else
         call CommonMPI_abort('Error in Flash: unknown context')

      endif

   end subroutine DefFlash_Flash_cv

end module DefFlash
