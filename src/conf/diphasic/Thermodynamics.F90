!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Thermodynamics

  use DefModel

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
      Thermodynamics_Psat, &
      Thermodynamics_Tsat
#endif

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
  ! PHASE_GAS = 1; PHASE_WATER = 2
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

    IF(iph == PHASE_GAS)THEN
       f = P

       dPf = 1.d0
       dTf = 0.d0
       dCf = 0.d0
       dSf = 0.d0
    ELSE IF(iph == PHASE_WATER)THEN
      IF(icp==1)THEN
        CALL air_henry(T,f)
        dPf = 0.d0
        CALL air_henry_dT(dTf)
        dCf = 0.d0
        dSf = 0.d0
      ELSE
        CALL f_PressionCapillaire(rt,iph,S,Pc,DSPc)
        CALL Thermodynamics_Psat(T, Psat, dTSat)

        f = Psat * DEXP(Pc/(T*RZetal))

        dPf = 0.d0
        dTf = (dTSat - Psat*Pc/RZetal/(T**2))*DEXP(Pc/(T*RZetal))
        dCf = 0.d0
        dSf = DSPc * f/(T*RZetal)
      ENDIF
    END IF
   end subroutine f_Fugacity


   SUBROUTINE air_henry(T,H)

     DOUBLE PRECISION, INTENT(IN) :: T
     DOUBLE PRECISION, INTENT(OUT) :: H

     DOUBLE PRECISION :: T1
     DOUBLE PRECISION :: T2
     DOUBLE PRECISION :: H1
     DOUBLE PRECISION :: H2

     T1 = 293.d0
     T2 = 353.d0

     H1 = 6.d+9
     H2 = 10.d+9

     H = H1 + (H2-H1)*(T-T1)/(T2-T1)

   END SUBROUTINE


   SUBROUTINE air_henry_dT(H_dt)

     DOUBLE PRECISION, INTENT(OUT) :: H_dt

     DOUBLE PRECISION :: T1
     DOUBLE PRECISION :: T2
     DOUBLE PRECISION :: H1
     DOUBLE PRECISION :: H2

     T1 = 293.d0
     T2 = 353.d0

     H1 = 6.d+9
     H2 = 10.d+9

     H_dt = (H2-H1)/(T2-T1)

   END SUBROUTINE


  ! Densite molaire
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_DensiteMolaire(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    INTEGER, INTENT(IN) :: iph
    DOUBLE PRECISION, INTENT(IN) :: P, T, C(NbComp), S(NbPhase)

    ! output
    DOUBLE PRECISION, INTENT(OUT) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    DOUBLE PRECISION :: Rgp
    DOUBLE PRECISION :: Pg0, T0

    Rgp = 8.314d0

    IF(iph==PHASE_GAS)THEN
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


  SUBROUTINE air_MasseMolaire(m)

    DOUBLE PRECISION, INTENT(OUT) :: m

    m = 29d-3
  END SUBROUTINE air_MasseMolaire


  SUBROUTINE H2O_MasseMolaire(m)

    DOUBLE PRECISION, INTENT(OUT) :: m

    m = 18d-3
  END SUBROUTINE H2O_MasseMolaire


  ! Densite Massique
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_DensiteMassique(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: zeta, air_m, H2O_m, m

    CALL f_DensiteMolaire(iph,P,T,C,S,zeta,dPf,dTf,dCf,dSf)

    CALL air_MasseMolaire(air_m)
    CALL H2O_MasseMolaire(H2O_m)

    m = C(1)*air_m + C(2)*H2O_m

    f = zeta * m

    dPf = dPf * m
    dTf = dTf * m
    dCf(1) = dCf(1) * m + zeta * air_m
    dCf(2) = dCf(2) * m + zeta * H2O_m
    dSf = dSf * m

  end subroutine f_DensiteMassique


  ! Viscosities
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_Viscosite(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: a, ss, ds, Cs, T1

    IF(iph==PHASE_GAS)THEN
      f = 18.51d-6
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
  ! PHASE_GAS = 1; PHASE_WATER = 2
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

    IF(rt(1) == 1)THEN
      Slrk = 0.4d0
      Sgrk = 0.d0
      rNvgk = 1.49d0
      rMvgk = 1.d0-1.d0/rNvgk
    ELSEIF( rt(1) == 2 )THEN
      Slrk = 0.01d0
      Sgrk = 0.d0
      rNvgk = 1.54d0
      rMvgk = 1.d0-1.d0/rNvgk
    ELSE
      PRINT*, 'error in f_PermRel, unknown rocktype'
      STOP
    ENDIF

    IF(iph==PHASE_GAS)THEN
      Sbl = (1.d0-S(iph)-Slrk)/(1.d0-Slrk-Sgrk)
      ss = 1.d0 - (Sbl)**(1.d0/rMvgk)
      if (S(iph).le.Sgrk) then
        f = 0.d0
        dSf(iph) = 0.d0
      else if (S(iph).ge.1.d0-Slrk) then
        f = 1.d0
        dSf(iph) = 0.d0
      else
        f = dsqrt(1.d0-Sbl)*( 1.d0-(Sbl)**(1.d0/rMvgk) )**(rMvgk*2.d0)

        dSf(iph) = ( 1.d0/2.d0/dsqrt(1.d0-Sbl)*ss**(rMvgk*2.d0) &
          + 2.d0*dsqrt(1.d0-Sbl)*ss**(2.d0*rMvgk-1.d0) &
          *Sbl**(1.d0/rMvgk-1) )/(1.d0-Slrk-Sgrk)
      endif
    ELSE
      Sbl = (S(iph)-Slrk)/(1.d0-Slrk-Sgrk)
      ss = 1.d0 - ( 1.d0 - (Sbl)**(1.d0/rMvgk) )**rMvgk
      c = 1.0E-2
      Sblc = (1.d0-c-Slrk)/(1.d0-Slrk-Sgrk)
      rkrlc = dsqrt(Sblc)*( 1.d0 - (1.d0-(Sblc)**(1/rMvgk) )**rMvgk)**2
      if (S(iph).le.Slrk) then
         f = 0.d0
         dSf(iph) = 0.d0
      else if ( (S(iph).ge.Slrk).and.(S(iph).le.1-c) ) then
         f = dsqrt(Sbl)*(1.d0 - (1.d0-(Sbl)**(1.d0/rMvgk))**rMvgk)**2
         ds = Sbl**(1.d0/rMvgk-1.d0) &
           *( 1.d0-(Sbl)**(1.d0/rMvgk) )**(rMvgk-1.d0)
         dSf(iph) = ( ss**2/(2.d0*dsqrt(Sbl)) + 2.d0*ss*ds*dsqrt(Sbl) ) &
               /(1.d0-Slrk-Sgrk)
      else if ( (S(iph).gt.1.0-c).and.(S(iph).lt.1.0-Sgrk) ) then
        f = (1.d0 - rkrlc)/c*S(iph) - (1.d0 - c - rkrlc)/c
        dSf(iph) = (1.d0 - rkrlc )/c
      else
        f = 1.d0
        dSf(iph) = 0.d0
      endif
    ENDIF
  END SUBROUTINE f_PermRel


  ! Pressions Capillaires des Phases et leurs derivees
  subroutine f_PressionCapillaire(rt,iph,S,f,DSf)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    double precision :: Slrk,Prvgk,rNvgk,rMvgk
    double precision :: Sbl, dSbl, Sbl1, Pc1, Sgrk

    IF(rt(1) == 1)THEN
      Slrk = 0.4d0
      Sgrk = 0.d0
      rNvgk = 1.49d0
      Prvgk = 15.d+6
      rMvgk = 1.d0-1.d0/rNvgk
    ELSEIF( rt(1) == 2 )THEN
      Slrk = 0.01d0
      Sgrk = 0.d0
      rNvgk = 1.54d0
      Prvgk = 2.d+6
      rMvgk = 1.d0-1.d0/rNvgk
    ELSE
      PRINT*, 'error in f_PressionCapillaire, unknown rocktype'
      STOP
    ENDIF

    dSf = 0.d0
    IF(iph==PHASE_GAS)THEN
      f = 0.d0
    ELSE
      Sbl = (S(iph)-Slrk)/(1.d0-Slrk-Sgrk)
      dSbl = 1.d0/(1.d0-Slrk-Sgrk)

      Sbl1 = 0.999d0    ! 0.97
      Pc1 = Prvgk*( Sbl1**(-1.d0/rMvgk) - 1.d0 )**(1.d0/rNvgk)

      if (Sbl < 1.0E-7) then
        write(*,*)' Sbl < 1.0E-7 ',Sbl,S(iph)
        stop
      ELSE IF (Sbl < Sbl1) THEN
        f = Prvgk*( Sbl**(-1.d0/rMvgk) - 1.d0 )**(1.d0/rNvgk)
        dSf(iph) = - Prvgk*dSbl/(rNvgk*rMvgk)*Sbl**(-1.d0/rMvgk-1.d0) &
          *( Sbl**(-1.d0/rMvgk) -1.d0 )**(1.d0/rNvgk-1.d0)
      ELSE IF  ( Sbl <= 1.0 ) THEN
        f = Pc1/(Sbl1-1.d0)*(Sbl-Sbl1) + Pc1
        dSf(iph) = Pc1/(Sbl1-1.d0)*dSbl
      ELSE
        f = 0.d0
        dSf(iph) = 0.d0
      ENDIF

      f = -f
      dSf(iph) = -dSf(iph)
    ENDIF
  END SUBROUTINE f_PressionCapillaire


  SUBROUTINE f_Sl(rt,Pc,Sl)

    INTEGER, INTENT(IN) :: rt(IndThermique+1)
    DOUBLE PRECISION, INTENT(IN) :: Pc

    DOUBLE PRECISION, INTENT(OUT) :: Sl

    DOUBLE PRECISION :: Slrk, Sgrk, Prvgk, rNvgk, rMvgk

    IF(rt(1) == 1)THEN
      Slrk = 0.4d0
      Sgrk = 0.d0
      rNvgk = 1.49d0
      Prvgk = 15.d+6
      rMvgk = 1.d0-1.d0/rNvgk
    ELSEIF( rt(1) == 2 )THEN
      Slrk = 0.01d0
      Sgrk = 0.d0
      rNvgk = 1.54d0
      Prvgk = 2.d+6
      rMvgk = 1.d0-1.d0/rNvgk
    ELSE
      PRINT*, 'error'
      STOP
    ENDIF


    IF(Pc < 0.d0) THEN
      Sl = 1.d0
    ELSE
      Sl = Slrk + (1.d0-Slrk)*( 1.d0 + (Pc/Prvgk)**rNvgk )**(-rMvgk)
    ENDIF
  END SUBROUTINE


#ifdef _THERMIQUE_

  ! EnergieInterne
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
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
  ! PHASE_GAS = 1; PHASE_WATER = 2
  ! If Enthalpide depends on the compositon C, change DefFlash.F90
  subroutine f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: a,b,cc,d,Ts,T0,ss,cp
    double precision :: H2O_m, air_m, m

    IF(iph==PHASE_GAS)THEN
      a = 1990.89d+3
      b = 190.16d+3
      cc = -1.91264d+3
      d = 0.2997d+3

      Ts = T/100.d0

      ss = a + b*Ts + cc*Ts**2 + d*Ts**3

      CALL H2O_MasseMolaire(H2O_m)

      CALL f_CpGaz(cp)
      CALL air_MasseMolaire(air_m)

      f = C(1)*cp*air_m*T + C(2)*ss*H2O_m

      dPf = 0.d0
      dTf = C(1)*cp*air_m + C(2)*H2O_m*(b + 2.d0*cc*Ts + 3.d0*d*Ts**2)/100.d0
      dCf(1) = cp*air_m*T
      dCf(2) = ss*H2O_m
      dSf = 0.d0
    ELSE
      a = -14.4319d+3
      b = +4.70915d+3
      cc = -4.87534
      d = 1.45008d-2
      T0 = 273.d0

      ss = a + b*(T-T0) + cc*(T-T0)**2 + d*(T-T0)**3

      CALL H2O_MasseMolaire(m)

      f = ss * m

      ss = b + 2*cc*(T-T0) + 3*d*(T-T0)**2

      dPf = 0.d0
      dTf = ss * m
      dCf = 0.d0
      dSf = 0.d0
    ENDIF
  end subroutine f_Enthalpie


  SUBROUTINE f_CpGaz(c)

    DOUBLE PRECISION, INTENT(OUT) :: c

    c = 1000.d0
  END SUBROUTINE

#endif

#ifdef _THERMIQUE_
  ! Compute Psat(T)
  subroutine Thermodynamics_Psat(T, Psat, dT_PSat)

    double precision, intent(in) :: T
    double precision, intent(out) :: Psat, dT_PSat

    Psat = 1.013E+5*dexp(13.7d0 -5120.d0/T)
    dT_PSat = 5120.d0*Psat/T**2

  end subroutine Thermodynamics_Psat


  ! Compute Tsat(P)
  subroutine Thermodynamics_Tsat(P, Tsat, dP_Tsat)

    double precision, intent(in) :: P
    double precision, intent(out) :: Tsat, dP_Tsat

    Tsat = - 5120.d0 / (DLOG(P/1.013E+5) - 13.7d0)
    dP_Tsat = 5120.d0 * ( 1.013E+5  / P ) / (DLOG(P/1.013E+5) - 13.7d0)**2

  end subroutine Thermodynamics_Tsat

#endif

end module Thermodynamics
