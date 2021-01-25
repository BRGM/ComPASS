
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
