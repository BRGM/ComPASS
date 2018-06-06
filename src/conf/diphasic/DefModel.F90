!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp thermal, MCP=(1,1,1,1)
!
! 1: Gas
! 2: Liquid
!
! 1: Air
! 2: H2O

module DefModel

  use CommonType

  implicit none

  ! ! ****** Model ****** ! !

  integer, parameter :: &
       NbComp = ComPASS_NUMBER_OF_COMPONENTS, &
       NbPhase = ComPASS_NUMBER_OF_PHASES

  integer, parameter :: &
      NbContexte = 2**NbPhase - 1

  ! MCP
  integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = RESHAPE( &
      (/ 1, 1, 1, 1 /), (/NbComp, NbPhase/))

  ! Phase PHASE_WATER is Liquid; PHASE_GAS is gas
  integer, parameter :: PHASE_GAS = 1
  integer, parameter :: PHASE_WATER = 2

  !FIXME: this is used for wells which are monophasic
  integer, parameter :: LIQUID_PHASE = PHASE_WATER

  ! Thermique
#ifdef _THERMIQUE_
  integer, parameter :: IndThermique = 1
#else
  integer, parameter :: IndThermique = 0
#endif


  ! ! ****** Constants derived from model (do not edit) ****** ! !

  ! Nombre Max d'eq d'equilibre
  !	       d'eq de fermeture thermodynamique
  !	       d'inc P (T) C
  !	       d'inc P (T) C primaires
  integer, parameter :: &
       NbEqEquilibreMax  = NbComp*(NbPhase-1),           & !< Max number of balance equations
       NbEqFermetureMax  = NbPhase + NbEqEquilibreMax,   & !< Max number of closure laws
       NbIncPTCMax       = 1 + IndThermique + sum(MCP),  &
       NbIncPTCSecondMax = NbEqFermetureMax,             &
       NbIncPTCSMax      = NbIncPTCMax + NbPhase,        &
       NbIncPTCSPrimMax  = NbComp + IndThermique,        &
       NbCompThermique   = NbComp + IndThermique


  ! ! ****** How to choose primary variables ****** ! !

  ! Served in module LoisthermoHydro.F90

  ! pschoice=1: manually
  !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
  !
  ! WARNING
  ! Il faut mettre les Sprim en dernier sinon il y a un pb qui reste a comprendre
  ! P est forcement primaire et en numero 1
  ! Si T est primaire elle doit etre en numero 2

  ! pschoice=2: Glouton method
  !     the matrix psprim and pssecd are defined formally for compile

  ! pschoise=3: Gauss method
  !     the matrix psprim and pssecd are defined formally for compile

  integer, parameter :: pschoice = 1



  integer, parameter, dimension( NbIncPTCSPrimMax, NbContexte) :: &
    psprim = RESHAPE( (/ &
      1, 2, 3, & ! ic=1
      1, 2, 4, & ! ic=2
      1, 2, 7  & ! ic=3
      /), (/ NbIncPTCSPrimMax, NbContexte /))

  integer, parameter, dimension( NbIncPTCSecondMax, NbContexte) :: &
       pssecd = RESHAPE( (/ &
       4, 0, 0, 0, & ! ic=1
       3, 0, 0, 0, & ! ic=2
       3, 4, 5, 6  & ! ic=3
       /), (/ NbIncPTCSecondMax, NbContexte/))

  ! ! ****** Alignment method ****** ! !

  ! Served in module Jacobian.F90

  ! aligmethod=1, manually
  !     it is necessary to give a three-dimension matrix: aligmat
  !     aligmat(:,:,ic) is the alignment matrix for context ic
  !     the index order of aligmat(:,:,ic) is (col,row), it allows us
  !     to define aligmat(:,:,ic) without considering that the matrix
  !     in Fortran is column-major
  !
  ! aligmethod=2, inverse diagnal
  !     it is necessary to define aligmat formally for compile

  integer, parameter :: aligmethod = 1

  double precision, parameter, &
       dimension( NbCompThermique, NbCompThermique, NbContexte) :: &
       aligmat = RESHAPE( (/ &
       1.d0, 1.d0, 0.d0, & ! ic=1
       0.d0, 0.d0, 1.d0, &
       1.d0, 0.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=2
       0.d0, 0.d0, 1.d0, &
       0.d0, 1.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=3
       0.d0, 0.d0, 1.d0, &
       0.d0, 1.d0, 0.d0  &
       /), (/ NbCompThermique, NbCompThermique, NbContexte /) )

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
      f_Enthalpie
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
        CALL DefModel_Psat(T, Psat, dTSat)

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

  !> \brief User set permeability
  !!
  !! \param[in] NbCellG,NbFracG Global number of cell and fracture face
  !! \param[in] IdCellG it is possible to set different permeability by rocktype
  !! \param[in,out] PermCellG Permeability tensor for each cell
  !! \param[in,out] PermFracG Permeability constant for each fracture face
  subroutine DefModel_SetPerm( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PermCellG, PermFracG)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG
    ! ouptuts:
    double precision, dimension(:,:,:), allocatable, intent(inout) :: &
      PermCellG
    double precision, dimension(:), allocatable, intent(inout) :: &
      PermFracG

    integer :: i

    ! allocate and set values
    allocate(PermCellG(3,3,NbCellG))
    do i=1, NbCellG
      PermCellG(:,:,i) = 0.d0

     IF(CellRocktypeG(1,i) == 1) THEN
        PermCellG(1,1,i) = 5.d-20
        PermCellG(2,2,i) = 5.d-20
        PermCellG(3,3,i) = 5.d-20
     ELSEIF(CellRocktypeG(1,i) == 2)THEN
       PermCellG(1,1,i) = 1.d-18
       PermCellG(2,2,i) = 1.d-18
       PermCellG(3,3,i) = 1.d-18
     ELSE
       PRINT*, 'error DefModel_SetPerm, unknow rocktype'
       PRINT*, i, CellRocktypeG(1,i)
       STOP
     ENDIF
    end do

    allocate(PermFracG(NbFracG))
    PermFracG(:) = 1.d-11
  end subroutine DefModel_SetPerm


  subroutine DefModel_SetPorosite( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      PorositeCell, PorositeFrac)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG
    ! ouptuts:
    double precision, dimension(:), allocatable, intent(inout) :: &
      PorositeCell
    double precision, dimension(:), allocatable, intent(inout) :: &
      PorositeFrac

    integer :: i

    allocate(PorositeCell(NbCellG))
    do i=1, NbCellG
      IF(CellRocktypeG(1,i) == 1) THEN
        PorositeCell(i) = 0.15d0
      ELSEIF(CellRocktypeG(1,i) == 2)THEN
        PorositeCell(i) = 0.3d0
      ELSE
        PRINT*, 'error in DefModel_SetPorosite, unknow CellRocktype'
        PRINT*, i, CellRocktypeG(1,i)
        STOP
      ENDIF
    end do

    allocate(PorositeFrac(NbFracG))
    do i=1, NbFracG
      IF(FracRocktypeG(1,i) == 1) THEN
        PorositeFrac(i) = 0.15d0
      ELSEIF(FracRocktypeG(1,i) == 2)THEN
        PorositeFrac(i) = 0.3d0
      ELSE
        PRINT*, 'error in DefModel_SetPorosite, unknow FracRocktype'
        PRINT*, i, FracRocktypeG(1,i)
        STOP
      ENDIF
    ENDDO
  end subroutine DefModel_SetPorosite


#ifdef _THERMIQUE_

  subroutine DefModel_SetCondThermique( &
      NbCellG, CellRocktypeG, NbFracG, FracRocktypeG, &
      CondThermalCellG, CondThermalFracG)

    integer, intent(in) :: NbCellG
    integer, dimension(:,:), intent(in) :: CellRocktypeG
    integer, intent(in) :: NbFracG
    integer, dimension(:,:), intent(in) :: FracRocktypeG

    ! ouptuts:
    double precision, dimension(:,:,:), allocatable, intent(inout) :: &
      CondThermalCellG
    double precision, dimension(:), allocatable, intent(inout) :: &
      CondThermalFracG

    integer :: i

    allocate(CondThermalCellG(3,3,NbCellG))
    do i=1, NbCellG
      CondThermalCellG(:,:,i) = 0.d0
      IF(CellRocktypeG(2,i) == 1) THEN
        CondThermalCellG(1,1,i) = 2.d0
        CondThermalCellG(2,2,i) = 2.d0
        CondThermalCellG(3,3,i) = 2.d0
      ELSEIF(CellRocktypeG(2,i) == 2)THEN
        CondThermalCellG(1,1,i) = 2.d0
        CondThermalCellG(2,2,i) = 2.d0
        CondThermalCellG(3,3,i) = 2.d0
      ELSE
        PRINT*, 'error in DefModel_SetCondThermique, unknow FracRocktype'
        PRINT*, i, FracRocktypeG(2,i)
        STOP
      ENDIF
    end do

    allocate(CondThermalFracG(NbFracG))
    CondThermalFracG(:) = 2.d0

  end subroutine DefModel_SetCondThermique


  !SUBROUTINE DefModel_SetThermalSource( &
  !    NbCell, &
  !    CellThermalSourceType, &
  !    NbFrac, &
  !    FracThermalSourceType, &
  !    CellThermalSource, &
  !    FracThermalSource)
  !
  !  INTEGER, INTENT(IN) :: NbCell
  !  INTEGER, DIMENSION(:), INTENT(IN) :: CellThermalSourceType
  !  INTEGER, INTENT(IN) :: NbFrac
  !  INTEGER, DIMENSION(:), INTENT(IN) :: FracThermalSourceType
  !  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CellThermalSource
  !  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: FracThermalSource
  !
  !  ALLOCATE(CellThermalSource(NbCell))
  !  CellThermalSource = MERGE(30.d0, 0.d0, CellThermalSourceType == 1)
  !
  !  ALLOCATE(FracThermalSource(NbFrac))
  !  FracThermalSource = 0.d0
  !END SUBROUTINE DefModel_SetThermalSource


  ! Compute Psat(T)
  subroutine DefModel_Psat(T, Psat, dT_PSat)

    double precision, intent(in) :: T
    double precision, intent(out) :: Psat, dT_PSat

    Psat = 1.013E+5*dexp(13.7d0 -5120.d0/T)
    dT_PSat = 5120.d0*Psat/T**2

  end subroutine DefModel_Psat


  ! Compute Tsat(P)
  subroutine DefModel_Tsat(P, Tsat, dP_Tsat)

    double precision, intent(in) :: P
    double precision, intent(out) :: Tsat, dP_Tsat

    Tsat = - 5120.d0 / (DLOG(P/1.013E+5) - 13.7d0)
    dP_Tsat = 5120.d0 * ( 1.013E+5  / P ) / (DLOG(P/1.013E+5) - 13.7d0)**2

  end subroutine DefModel_Tsat

#endif

end module DefModel
