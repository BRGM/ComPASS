!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Thermodynamics

  use iso_c_binding, only: c_double, c_int
  use DefModel, only: NbPhase, NbComp, IndThermique, &
    GAS_PHASE, LIQUID_PHASE, WATER_COMP, AIR_COMP
  use CommonMPI, only: CommonMPI_abort

  implicit none

  public :: &
      f_Fugacity,           &  ! Fucacity
      f_DensiteMolaire,     &  ! \xi^alpha(P,T,C,S)
      f_DensiteMassique,    &  ! \rho^alpha(P,T,C,S)
      f_Viscosite,          &  ! \mu^alpha(P,T,C,S)
      f_PermRel,            &  ! k_{r_alpha}(S)
      f_PressionCapillaire, &  ! P_{c,alpha}(S)
      air_henry,            &
      air_henry_dT,         &
      liquid_pressure

#ifdef _THERMIQUE_
  public :: &
      f_EnergieInterne,     &
      f_Enthalpie, &
      f_SpecificEnthalpy, &
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat
#endif

contains

  function liquid_pressure(z_ref, p_ref, rho, g, z)
    real(c_double), intent(in) :: z_ref
    real(c_double), intent(in) :: p_ref
    real(c_double), intent(in) :: rho
    real(c_double), intent(in) :: g
    real(c_double), intent(in) :: z
    real(c_double) :: liquid_pressure

    liquid_pressure = p_ref - rho * g * (z - z_ref)
  end function liquid_pressure

  ! *** Physics *** !

  ! Fugacity coefficient
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  ! P is the pressure of the reference phase (P=Pg)
  subroutine f_Fugacity(rt,iph,icp,P,T,C,S,f,DPf,DTf,DCf,DSf)

    ! input
    integer(c_int), intent(in) :: rt(IndThermique+1)
    integer(c_int), intent(in) :: iph, icp
    real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp),DSf(NbPhase)

    real(c_double) :: PSat, dTSat, Pc, DSPc(NbPhase)
    real(c_double) :: RZetal

    RZetal = 8.314d0 * 1000.d0 / 0.018d0

    dPf = 0.d0
    dTf = 0.d0
    dCf = 0.d0
    dSf = 0.d0

    if(iph == GAS_PHASE)then
       f = P

       dPf = 1.d0
    else if(iph == LIQUID_PHASE)then
      if(icp==AIR_COMP)then
        call air_henry(T,f)
        call air_henry_dT(dTf)
      else if(icp==WATER_COMP)then
        call f_PressionCapillaire(rt,iph,S,Pc,DSPc)  ! Pl=Pref + Pc, then Pc = - capillary pressure ! 
        call FluidThermodynamics_Psat(T, Psat, dTSat) 

        f = Psat * dexp(Pc/(T*RZetal))

        dTf = (dTSat - Psat*Pc/RZetal/(T**2))*dexp(Pc/(T*RZetal))
        dSf = DSPc * f/(T*RZetal)
      endif
    endif
   end subroutine f_Fugacity


   subroutine air_henry(T,H)

     real(c_double), intent(in) :: T
     real(c_double), intent(out) :: H

     real(c_double) :: T1, T2
     real(c_double) :: H1, H2

     T1 = 293.d0
     T2 = 353.d0

     H1 = 6.d+9
     H2 = 10.d+9

     H = H1 + (H2-H1)*(T-T1)/(T2-T1)

   end subroutine


   subroutine air_henry_dT(H_dt)

     real(c_double), intent(out) :: H_dt

     real(c_double) :: T1, T2
     real(c_double) :: H1, H2

     T1 = 293.d0
     T2 = 353.d0

     H1 = 6.d+9
     H2 = 10.d+9

     H_dt = (H2-H1)/(T2-T1)

   end subroutine

  subroutine air_MasseMolaire(m)

    real(c_double), intent(out) :: m

    m = 29.d-3
  end subroutine air_MasseMolaire


  subroutine H2O_MasseMolaire(m)

    real(c_double), intent(out)  :: m

    m = 18.d-3
  end subroutine H2O_MasseMolaire

  ! Densite molaire
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_DensiteMolaire(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_density")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: Rgp, H2O_m

    Rgp = 8.314d0
    call H2O_MasseMolaire(H2O_m)

    if(iph==GAS_PHASE)then
      f = P/(Rgp*T)

      dPf = 1/(Rgp*T)
      dTf = -P/Rgp/T**2
      dCf = 0.d0
      dSf = 0.d0
    else if(iph == LIQUID_PHASE)then
      f = 1000.d0/H2O_m

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0
    endif
  end subroutine f_DensiteMolaire


  ! Densite Massique
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_DensiteMassique(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer(c_int), intent(in) :: iph
    real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    real(c_double) :: zeta, air_m, H2O_m, m

    call f_DensiteMolaire(iph,P,T,C,S,zeta,dPf,dTf,dCf,dSf)

    call air_MasseMolaire(air_m)
    call H2O_MasseMolaire(H2O_m)

    m = C(AIR_COMP)*air_m + C(WATER_COMP)*H2O_m

    f = zeta * m

    dPf = dPf * m
    dTf = dTf * m
    dCf(AIR_COMP) = dCf(AIR_COMP) * m + zeta * air_m
    dCf(WATER_COMP) = dCf(WATER_COMP) * m + zeta * H2O_m
    dSf = dSf * m

  end subroutine f_DensiteMassique


  ! Viscosities
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_Viscosite(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")

    ! input
    integer(c_int), value, intent(in) :: iph
    real(c_double), value, intent(in) :: P, T
    real(c_double), intent(in) :: C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: &
        f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    if(iph==GAS_PHASE)then
      f = 2.d-5
    else if(iph == LIQUID_PHASE)then
      f = 1.d-3
    endif

    dPf = 0.d0
    dTf = 0.d0
    dCf = 0.d0
    dSf = 0.d0

  end subroutine f_Viscosite


  ! Permeabilites = S**2
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_PermRel(rt,iph,S,f,DSf)

    ! input
    integer(c_int), intent(in) :: rt(IndThermique + 1)
    integer(c_int), intent(in) :: iph
    real(c_double), intent(in) :: S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, DSf(NbPhase)

    f = S(iph)**2
    dSf = 0.d0
    dSf(iph) = 2.d0*S(iph)
    if(S(iph)<0.d0) stop 2
    if(S(iph)>1.d0) stop 2


  end subroutine f_PermRel


  ! Pressions Capillaires des Phases et leurs derivees
  subroutine f_PressionCapillaire(rt,iph,S,f,DSf)

    ! input
    integer(c_int), intent(in) :: rt(IndThermique + 1)
    integer(c_int), intent(in) :: iph
    real(c_double), intent(in) :: S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, DSf(NbPhase)

    real(c_double) :: Pc_cst, Sg0, A
    real(c_double) :: Sg

    dSf = 0.d0
    if(iph==GAS_PHASE)then
      f = 0.d0
    else if(iph == LIQUID_PHASE)then
        Sg = 1.d0 - S(LIQUID_PHASE)
        Pc_cst = 2.d5
        Sg0 = 1.d0 - 1.d-2
        A = - Pc_cst * log(1.d0-Sg0) - Pc_cst/(1.d0-Sg0) * Sg0
        if( Sg < Sg0 )   then
            f = - Pc_cst * log(1.d0-Sg)
            dSf(iph) = Pc_cst/(1.d0-Sg) ! wrt Sg
        else
            f = Pc_cst * Sg / ( 1.d0 - Sg0 ) + A
            dSf(iph) = Pc_cst / ( 1.d0 - Sg0 )  ! wrt Sg
        endif
        dSf(iph) = -dSf(iph) ! wrt Sl = 1 - Sg
    endif


  end subroutine f_PressionCapillaire


  subroutine f_Sl(rt,Pc,Sl)

    integer(c_int), INTENT(IN) :: rt(IndThermique+1)
    real(c_double), INTENT(IN) :: Pc

    real(c_double), INTENT(OUT) :: Sl

    call CommonMPI_abort('entered in f_Sl, but Pc=0')

  end subroutine


#ifdef _THERMIQUE_

  ! EnergieInterne
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
subroutine f_EnergieInterne(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer(c_int), intent(in) :: iph
    real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    CALL f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

  end subroutine f_EnergieInterne


  ! Enthalpie
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")

    ! input
    integer(c_int), value, intent(in) :: iph
    real(c_double), value, intent(in) :: P, T
    real(c_double), intent(in) :: C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    real(c_double) :: H2O_m, air_m
    real(c_double) :: a,b,cc,d,Ts,T0,ss,dTss,cp

    if(iph==GAS_PHASE)then
      a = 1990.89d+3
      b = 190.16d+3
      cc = -1.91264d+3
      d = 0.2997d+3

      Ts = T/100.d0

      ss = a + b*Ts + cc*Ts**2.d0 + d*Ts**3.d0
      dTss = (b + 2.d0*cc*Ts + 3.d0*d*Ts**2.d0)/100.d0

      call H2O_MasseMolaire(H2O_m)

      call f_CpGaz(cp)
      call air_MasseMolaire(air_m)

      f = C(AIR_COMP)*cp*air_m*T + C(WATER_COMP)*ss*H2O_m

      dPf = 0.d0
      dTf = C(AIR_COMP)*cp*air_m + C(WATER_COMP)*H2O_m*dTss
      dCf(AIR_COMP) = cp*air_m*T
      dCf(WATER_COMP) = ss*H2O_m
      dSf = 0.d0

    else if(iph == LIQUID_PHASE)then

      a = -14.4319d+3
      b = +4.70915d+3
      cc = -4.87534
      d = 1.45008d-2
      T0 = 273.d0

      ss = a + b*(T-T0) + cc*(T-T0)**2.d0 + d*(T-T0)**3.d0
      dTss = b + 2.d0*cc*(T-T0) + 3.d0*d*(T-T0)**2.d0

      call H2O_MasseMolaire(H2O_m)

      f = ss * H2O_m

      dPf = 0.d0
      dTf = dTss * H2O_m
      dCf = 0.d0
      dSf = 0.d0
    endif

  end subroutine f_Enthalpie

  ! Specific Enthalpy (used in FreeFlow)
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_SpecificEnthalpy(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")

    ! input
    integer(c_int), value, intent(in) :: iph
    real(c_double), value, intent(in) :: P, T
    real(c_double), intent(in) :: C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp), &
                                  dCf(NbComp, NbComp), dSf(NbComp, NbPhase)

    real(c_double) :: H2O_m, air_m
    real(c_double) :: a,b,cc,d,Ts,T0,ss,dTss,cp

    if(iph == GAS_PHASE)then
      a = 1990.89d+3
      b = 190.16d+3
      cc = -1.91264d+3
      d = 0.2997d+3

      Ts = T/100.d0

      ss = a + b*Ts + cc*Ts**2.d0 + d*Ts**3.d0
      dTss = (b + 2.d0*cc*Ts + 3.d0*d*Ts**2.d0)/100.d0


      call H2O_MasseMolaire(H2O_m)

      call f_CpGaz(cp)
      call air_MasseMolaire(air_m)

      f(AIR_COMP) = cp*air_m*T
      f(WATER_COMP) = ss*H2O_m

      dPf = 0.d0
      dTf(AIR_COMP) = cp*air_m
      dTf(WATER_COMP) = H2O_m*dTss
      dCf = 0.d0
      dSf = 0.d0

    else if(iph == LIQUID_PHASE)then

      a = -14.4319d+3
      b = +4.70915d+3
      cc = -4.87534
      d = 1.45008d-2
      T0 = 273.d0

      ss = a + b*(T-T0) + cc*(T-T0)**2.d0 + d*(T-T0)**3.d0
      dTss = b + 2.d0*cc*(T-T0) + 3.d0*d*(T-T0)**2.d0

      call H2O_MasseMolaire(H2O_m)

      f(:) = ss * H2O_m

      dPf = 0.d0
      dTf(:) = dTss * H2O_m
      dCf = 0.d0
      dSf = 0.d0
    endif

  end subroutine f_SpecificEnthalpy

  subroutine f_CpGaz(c)

    real(c_double), intent(out) :: c

    c = 1000.d0
  end subroutine

  ! Compute Psat(T)
  subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

    real(c_double), value, intent(in) :: T
    real(c_double), intent(out) :: Psat, dT_PSat

! valid between -50C and 200C
    Psat = 100.d0 * dexp(46.784d0 - 6435.d0/T - 3.868d0 * dlog(T))
    dT_PSat = (6435.d0/T**2.d0 - 3.868d0/T) * Psat

  end subroutine FluidThermodynamics_Psat


  ! Compute Tsat(P)
  subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

    real(c_double), value, intent(in) :: P
    real(c_double), intent(out) :: Tsat, dP_Tsat

    Tsat = 0.d0
    dP_Tsat = 0.d0

    call CommonMPI_abort('entered in Tsat, not implemented')

  end subroutine FluidThermodynamics_Tsat

#endif

end module Thermodynamics
