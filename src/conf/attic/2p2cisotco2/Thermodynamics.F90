!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp thermal, MCP=(1,1; 1, 1)

! 1: Gas
! 2: Water

module Thermodynamics
 
   use DefModel, only: NbPhase, NbComp, IndThermique


   implicit none

     public :: &
      f_Fugacity,           &  ! Fucacity
      f_DensiteMolaire,     &  ! \xi^alpha(P,T,C,S)
      f_DensiteMassique,    &  ! \rho^alpha(P,T,C,S)
      f_Viscosite,          &  ! \mu^alpha(P,T,C,S)
      f_PermRel,            &  ! k_{r_alpha}(S)
      f_PressionCapillaire, &  ! P_{c,alpha}(S)
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat, &
      CO2_henry,            &
      liquid_pressure

!#ifdef _THERMIQUE_

  public :: &
      f_EnergieInterne,     &
      f_Enthalpie
!#endif

contains

  FUNCTION liquid_pressure(z_ref, p_ref, rho, g, z)
    DOUBLE PRECISION, INTENT(IN) :: z_ref
    DOUBLE PRECISION, INTENT(IN) :: p_ref
    DOUBLE PRECISION, INTENT(IN) :: rho
    DOUBLE PRECISION, INTENT(IN) :: g
    DOUBLE PRECISION, INTENT(IN) :: z
    DOUBLE PRECISION :: liquid_pressure

    liquid_pressure = p_ref - rho * g * (z - z_ref)
  END FUNCTION liquid_pressure

  ! *** Physics *** !

  ! Fugacity coefficient
  ! f * c_i
  ! P = Pg pression de reference
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_Fugacity(rt,iph,icp,P,T,C,S,f,DPf,DTf,DCf,DSf)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    integer, intent(in) :: iph, icp
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, DPf, DTf, DCf(NbComp),DSf(NbPhase)

    double precision :: PSat, dTSat, Pc, DSPc(NbPhase)
    DOUBLE PRECISION :: RZetal

    RZetal = 8.314d0 * 1000.d0 / 0.018d0

    IF(iph == GAS_PHASE)THEN
       f = P

       dPf = 1.d0
       dTf = 0.d0
       dCf = 0.d0
       dSf = 0.d0
    ELSE IF(iph == LIQUID_PHASE)THEN
      IF(icp==1)THEN
        CALL CO2_henry(T,f, dTf)
        dPf = 0.d0
        dCf = 0.d0
        dSf = 0.d0
      ELSE
        CALL FluidThermodynamics_Psat(T, Psat, dTSat)

        f = Psat

        dPf = 0.d0
        dTf = dTSat
        dCf = 0.d0
        dSf = 0.d0
      ENDIF
    END IF
   end subroutine f_Fugacity


   SUBROUTINE CO2_henry(T,H, H_dt)

     DOUBLE PRECISION, INTENT(IN) :: T
     DOUBLE PRECISION, INTENT(OUT) :: H, H_dt

     DOUBLE PRECISION :: T1
     DOUBLE PRECISION :: T2
     DOUBLE PRECISION :: H1
     DOUBLE PRECISION :: H2

     T1 = 298.d0
     T2 = 353.d0

     H1 = 1.6e8
     H2 = 5.6e8

     H = H1 + (H2-H1)*(T-T1)/(T2-T1)
     H_dt = (H2-H1)/(T2-T1)

   END SUBROUTINE


  ! Densite molaire
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_DensiteMolaire(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_density")

    ! input
    INTEGER, INTENT(IN) :: iph
    DOUBLE PRECISION, INTENT(IN) :: P, T, C(NbComp), S(NbPhase)

    ! output
    DOUBLE PRECISION, INTENT(OUT) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    DOUBLE PRECISION :: Rgp
    DOUBLE PRECISION :: Pg0, T0

    Rgp = 8.314d0

    IF(iph==GAS_PHASE)THEN
      f = P/(Rgp*T)

      dPf = 1/(Rgp*T)
      dTf = -P/Rgp/T**2
      dCf = 0.d0
      dSf = 0.d0
    ELSE
      f = 1000.d0/0.018d0

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0
    ENDIF
  end subroutine f_DensiteMolaire


  SUBROUTINE CO2_MasseMolaire(m)

    DOUBLE PRECISION, INTENT(OUT) :: m

    m = 44d-3
  END SUBROUTINE CO2_MasseMolaire


  SUBROUTINE H2O_MasseMolaire(m)

    DOUBLE PRECISION, INTENT(OUT) :: m

    m = 18d-3
  END SUBROUTINE H2O_MasseMolaire


  ! Densite Massique
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_DensiteMassique(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: zeta, CO2_m, H2O_m, m

    CALL f_DensiteMolaire(iph,P,T,C,S,zeta,dPf,dTf,dCf,dSf)

    CALL CO2_MasseMolaire(CO2_m)
    CALL H2O_MasseMolaire(H2O_m)

    m = C(1)*CO2_m + C(2)*H2O_m

    f = zeta * m

    dPf = dPf * m
    dTf = dTf * m
    dCf(1) = dCf(1) * m + zeta * CO2_m
    dCf(2) = dCf(2) * m + zeta * H2O_m
    dSf = dSf * m

  end subroutine f_DensiteMassique


  ! Viscosities
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_Viscosite(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: a, ss, ds, Cs, T1

    IF(iph==GAS_PHASE)THEN
      f = 15.d-6    ! https://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html
    ELSE
      f = 1.0d-3
    ENDIF

    dPf = 0.d0
    dTf = 0.d0
    dCf = 0.d0
    dSf = 0.d0

  end subroutine f_Viscosite


  ! Permeabilites = S**2
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_PermRel(rt,iph,S,f,DSf)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    double precision :: Slrk, Sgrk, rNvgk, rMvgk
    double precision :: Sbl, ss, Sblc, rkrlc, ds
    double precision :: c

    dSf = 0.d0

    if (iph == GAS_PHASE) then
       f = S(iph)**2
       dSf(iph) = 2*S(iph)
    else if (iph == LIQUID_PHASE) then
       f  = S(iph)**2
       dSf(iph) = 2* S(iph)
    else
       print*, 'Error in f_PermRel, unknown phase'
    endif
    
  END SUBROUTINE f_PermRel


  ! Pressions Capillaires des Phases et leurs derivees
  subroutine f_PressionCapillaire(rt,iph,S,f,DSf)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    f = 0.d0
    dSf = 0.d0
  END SUBROUTINE f_PressionCapillaire




!#ifdef _THERMIQUE_

  ! EnergieInterne
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  subroutine f_EnergieInterne(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    CALL f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

  end subroutine f_EnergieInterne


  ! Enthalpie
  ! iph is an identificator for each phase:
  ! GAS_PHASE = 1; LIQUID_PHASE = 2
  ! If Enthalpide depends on the compositon C, change DefFlash.F90
  subroutine f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: a,b,cc,d,Ts,T0,ss,cp
    double precision :: H2O_m, air_m, m

    print*, 'Enthalpie should no tbe called'
    stop 431
  end subroutine f_Enthalpie


  SUBROUTINE f_CpGaz(c)

    DOUBLE PRECISION, INTENT(OUT) :: c

    c = 1000.d0
  END SUBROUTINE

!#endif


  ! Compute Psat(T)
  subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

    double precision, intent(in) :: T
    double precision, intent(out) :: Psat, dT_PSat

    Psat = 1.013E+5*dexp(13.7d0 -5120.d0/T)
    dT_PSat = 5120.d0*Psat/T**2

  end subroutine FluidThermodynamics_Psat


  ! Compute Tsat(P)
  subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

    double precision, intent(in) :: P
    double precision, intent(out) :: Tsat, dP_Tsat

    Tsat = - 5120.d0 / (DLOG(P/1.013E+5) - 13.7d0)
    dP_Tsat = 5120.d0 * ( 1.013E+5  / P ) / (DLOG(P/1.013E+5) - 13.7d0)**2

  end subroutine FluidThermodynamics_Tsat

end module Thermodynamics
