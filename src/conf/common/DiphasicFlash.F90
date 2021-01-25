
    pure subroutine DiphasicFlash_liquid_to_diphasic(inc, pa, dpadS)
       type(Type_IncCVReservoir), intent(inout) :: inc
       real(c_double), intent(in) :: pa(NbPhase) ! p^\alpha: phase pressure
       real(c_double), intent(in) :: dpadS(NbPhase)

       real(c_double) :: f
       real(c_double) :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! dummy values
       real(c_double) :: PgCag, PgCwg

       call f_Fugacity(LIQUID_PHASE, AIR_COMP, inc, pa, dpadS, f, DPf, DTf, DCf, DSf)
       PgCag = f*inc%Comp(AIR_COMP, LIQUID_PHASE)
       call f_Fugacity(LIQUID_PHASE, WATER_COMP, inc, pa, dpadS, f, DPf, DTf, DCf, DSf)
       PgCwg = f*inc%Comp(WATER_COMP, LIQUID_PHASE)
       ! WARNING: don't divide inequality by Pg (migth be negative during Newton iteration)
       if (PgCag + PgCwg > pa(GAS_PHASE)) then
          inc%ic = DIPHASIC_CONTEXT
          inc%Saturation(GAS_PHASE) = 0.d0
          inc%Saturation(LIQUID_PHASE) = 1.d0
       endif

    end subroutine DiphasicFlash_liquid_to_diphasic

    pure subroutine DiphasicFlash_gas_to_diphasic(inc, pa, dpadS)
       type(Type_IncCVReservoir), intent(inout) :: inc
       real(c_double), intent(in) :: pa(NbPhase) ! p^\alpha: phase pressure
       real(c_double), intent(in) :: dpadS(NbPhase)

       real(c_double) :: fg, fl
       real(c_double) :: dPf, dTf, dCf(NbComp), dSf(NbPhase) ! dummy values
       real(c_double) :: Cla, Clw

       call f_Fugacity(GAS_PHASE, AIR_COMP, inc, pa, dpadS, fg, DPf, DTf, DCf, DSf)
       call f_Fugacity(LIQUID_PHASE, AIR_COMP, inc, pa, dpadS, fl, DPf, DTf, DCf, DSf)
       Cla = (fg/fl)*inc%Comp(AIR_COMP, GAS_PHASE)
       call f_Fugacity(GAS_PHASE, WATER_COMP, inc, pa, dpadS, fg, DPf, DTf, DCf, DSf)
       call f_Fugacity(LIQUID_PHASE, WATER_COMP, inc, pa, dpadS, fl, DPf, DTf, DCf, DSf)
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
    pure subroutine DiphasicFlash_Flash_cv(inc, pa, dpadS)
       type(Type_IncCVReservoir), intent(inout) :: inc
       real(c_double), intent(in) :: pa(NbPhase) ! p^\alpha: phase pressure
       real(c_double), intent(in) :: dpadS(NbPhase)

       integer(c_int) :: context

       context = inc%ic

       if (context == LIQUID_CONTEXT) then

          call DiphasicFlash_liquid_to_diphasic(inc, pa, dpadS)

       elseif (context == DIPHASIC_CONTEXT) then

          call DiphasicFlash_diphasic_switches(inc)

       elseif (context == GAS_CONTEXT) then

          call DiphasicFlash_gas_to_diphasic(inc, pa, dpadS)

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
