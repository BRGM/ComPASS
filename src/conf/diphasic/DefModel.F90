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
       NbComp = 2, & !< Number of Component
       NbPhase = 2   !< Number of Phase

  integer, parameter :: &
      NbContexte = 2**NbPhase - 1

  ! MCP 
  integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = RESHAPE( &
      (/ 1, 1, 1, 1 /), (/NbComp, NbPhase/))

  ! Phase PHASE_WATER is Liquid; PHASE_GAS is gas
  integer, parameter :: PHASE_GAS = 1
  integer, parameter :: PHASE_WATER = 2

  ! Gravite
  double precision, parameter :: Gravite = 0.d0 !< Gravity constant
  
  ! CpRoche
  double precision, parameter :: CpRoche = 2000.d0*1000.d0 !< en volumique

  ! thickness of frac
  double precision, parameter :: Thickness = 1.d0 !< Thickness of the fractures
  
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
  ! pschoice=2: Glouton method
  !     the matrix psprim and pssecd are defined formally for compile
  ! pschoise=3: Gauss method
  !     the matrix psprim and pssecd are defined formally for compile

  integer, parameter :: pschoice = 1

  integer, parameter, dimension( NbIncPTCSPrimMax, NbContexte) :: &
    psprim = RESHAPE( (/ &
      1, 2, 4, & ! ic=1
      1, 2, 3, & ! ic=2
      1, 7, 5  & ! ic=3
      /), (/ NbIncPTCSPrimMax, NbContexte /))

  integer, parameter, dimension( NbIncPTCSecondMax, NbContexte) :: &
       pssecd = RESHAPE( (/ &
       3, 0, 0, 0, & ! ic=1
       4, 0, 0, 0, & ! ic=2
       2, 3, 4, 6  & ! ic=3
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

  integer, parameter :: aligmethod = 2

  double precision, parameter, &
       dimension( NbCompThermique, NbCompThermique, NbContexte) :: &
       aligmat = RESHAPE( (/ &
       1.d0, 1.d0, 0.d0, & ! ic=1
       0.d0, 0.d0, 1.d0, &
       0.d0, 1.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=2
       0.d0, 0.d0, 1.d0, &
       1.d0, 0.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=3
       0.d0, 0.d0, 1.d0, &
       1.d0, 0.d0, 0.d0  &
       /), (/ NbCompThermique, NbCompThermique, NbContexte /) )


  ! ! ******* Mesh type ****** ! !
  ! ! * cartesian-quad
  ! ! * hexahedron-quad
  ! ! * tetrahedron-triangle
  ! ! * wedge
  character(len=40), parameter :: &
      MESH_TYPE = "cartesian-quad"


  ! ! ****** Times ****** ! !

  ! One day
  double precision, parameter :: OneSecond = 1.d0
  double precision, parameter :: OneMinute = 60.d0 * OneSecond
  double precision, parameter :: OneHour = 60.d0 * OneMinute
  double precision, parameter :: OneDay = 24.d0 * OneHour
  double precision, parameter :: OneMonth = 30.5d0 * OneDay
  double precision, parameter :: OneYear = 3.6525d2 * OneDay

  ! time step init and max 
  ! FIXME: parameter is removed to assign variable from python
  double precision :: TimeFinal = 200 * OneYear

  double precision, parameter :: TimeStepInit = OneSecond
  double precision, parameter :: TimeStepMax1 = OneYear

  ! output_frequency for visu
  double precision, parameter :: output_frequency = OneSecond


  ! ! ****** Newton iters max and stop condition ****** ! !   
  integer, parameter :: NewtonNiterMax = 40
  double precision, parameter :: NewtonTol = 1.d-10


  ! ! ****** ksp linear solver iters max and stop condition ****** ! !
  integer, parameter :: KspNiterMax = 150        ! max nb of iterations
  double precision, parameter :: KspTol = 1.d-12  ! tolerance


  ! ! ****** Obj values used to compute Newton increment ****** ! !

  double precision, parameter :: NewtonIncreObj_P = 2.d5
  double precision, parameter :: NewtonIncreObj_T = 20.d0
  double precision, parameter :: NewtonIncreObj_C = 0.4d0
  double precision, parameter :: NewtonIncreObj_S = 0.3d0


  ! ! ****** Obj values used to compute next time step ****** ! !

  double precision, parameter :: TimeStepObj_P = 5.d5
  double precision, parameter :: TimeStepObj_T = 40.d0
  double precision, parameter :: TimeStepObj_C = 0.9d0
  double precision, parameter :: TimeStepObj_S = 0.9d0


  ! ! ****** Parameters of VAG schme (volume distribution) ****** ! !

  double precision, parameter :: &
      omegaDarcyCell = 0.075,    & ! darcy cell/frac
      omegaDarcyFrac = 0.15

  double precision, parameter :: & 
      omegaFourierCell = 0.075,  & ! fourier cell/frac 
      omegaFourierFrac = 0.15


  ! ! ****** Others ****** ! !

  ! eps
  double precision, parameter :: eps = 1.0d-10

  public :: &
      f_Fugacity,           &  ! Fucacity
      f_DensiteMolaire,     &  ! \xi^alpha(P,T,C,S)
      f_DensiteMassique,    &  ! \rho^alpha(P,T,C,S)
      f_Viscosite,          &  ! \mu^alpha(P,T,C,S)
      f_PermRel,            &  ! k_{r_alpha}(S)
      f_PressionCapillaire     ! P_{c,alpha}(S)

#ifdef _THERMIQUE_

  public :: &
      f_EnergieInterne,     &
      f_Enthalpie
#endif

contains

  ! *** Physics *** !

  ! Fugacity coefficient
  ! f * c_i
  ! P = Pg pression de reference
  ! iph is an identificator for each phase: 
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_Fugacity(rocktype, iph,icp,P,T,C,S,f,DPf,DTf,DCf,DSf)

    ! input
    integer, intent(in) :: rocktype, iph, icp
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
        CALL f_PressionCapillaire(rocktype,iph,S,Pc,DSPc)
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
    
    Rgp = 8.314d0

    IF(iph==PHASE_GAS)THEN
      f = P/(Rgp*T)

      dPf = 1/(Rgp*T)
      dTf = P/Rgp * (-1/T**2)
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

    m = 29.d-3
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
  subroutine f_PermRel(iph,S,f,DSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    double precision :: Slrk, Sgrk, rNvgk, rMvgk
    double precision :: Sbl, ss, Sblc, rkrlc, ds
    double precision :: c

    dSf = 0.d0

    Slrk = 0.4d0 
    Sgrk = 0.d0 
    rNvgk = 1.49d0
    rMvgk = 1.d0-1.d0/rNvgk

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
      dSf(PHASE_WATER) = 0.d0
      dSf(PHASE_WATER) = -dSf(PHASE_GAS)
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
      dSf(PHASE_GAS) = 0.d0
      dSf(PHASE_GAS) = -DSf(PHASE_WATER)
    ENDIF
  END SUBROUTINE f_PermRel


  ! Pressions Capillaires des Phases et leurs derivees
  subroutine f_PressionCapillaire(rocktype,iph,S,f,DSf)

    ! input
    integer, intent(in) :: rocktype
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    double precision :: Slrk,Prvgk,rNvgk,rMvgk
    double precision :: Sbl, dSbl, Sbl1, Pc1, Sgrk

    IF(iph==PHASE_GAS)THEN
      f = 0.d0
      dSf = 0.d0
    ELSE
      Slrk = 0.4d0
      Prvgk = 15.d+6
      rNvgk = 1.49d0
      rMvgk = 1.d0-1.d0/rNvgk

      Sgrk = 0.d0

      Sbl = (S(PHASE_WATER)-Slrk)/(1.d0-Slrk-Sgrk)
      dSbl = 1.d0/(1.d0-Slrk-Sgrk)

      Sbl1 = 0.999d0    ! 0.97
      Pc1 = Prvgk*( Sbl1**(-1.d0/rMvgk) - 1.d0 )**(1.d0/rNvgk)

      if (Sbl < 1.0E-7) then
        write(*,*)' Sbl < 1.0E-7 ',Sbl,S(PHASE_WATER)  
        stop
      ELSE IF (Sbl < Sbl1) THEN
        f = Prvgk*( Sbl**(-1.d0/rMvgk) - 1.d0 )**(1.d0/rNvgk) 
        dSf(PHASE_WATER) = - Prvgk*dSbl/(rNvgk*rMvgk)*Sbl**(-1.d0/rMvgk-1.d0) &
          *( Sbl**(-1.d0/rMvgk) -1.d0 )**(1.d0/rNvgk-1.d0)
      ELSE IF  ( Sbl <= 1.0 ) THEN
        f = Pc1/(Sbl1-1.d0)*(Sbl-Sbl1) + Pc1
        dSf(PHASE_WATER) = Pc1/(Sbl1-1.d0)*dSbl
      ELSE
        f = 0.d0
        dSf(PHASE_WATER) = 0.d0
      ENDIF    
      dSf(PHASE_GAS) = 0.d0
      dSf(PHASE_GAS) = -dSf(PHASE_WATER) 

      f = -f
      dSf = -dSf
    ENDIF
  END SUBROUTINE f_PressionCapillaire


  SUBROUTINE f_Sl(Pc,Sl)

    DOUBLE PRECISION, INTENT(IN) :: Pc
    
    DOUBLE PRECISION, INTENT(OUT) :: Sl

    DOUBLE PRECISION :: Slrk, Prvgk, rNvgk, rMvgk

    Slrk = 0.4d0
    Prvgk = 15.d+6
    rNvgk = 1.49d0
    rMvgk = 1.d0-1.d0/rNvgk

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

    double precision :: a,b,cc,d,T0,ss
    double precision :: m

    IF(iph==PHASE_GAS)THEN
      CALL f_CpGaz(cc)
      CALL air_MasseMolaire(m)

      f = cc * m * T

      dPf = 0.d0
      dTf =  cc * m
      dTf = 0.d0
      dCf = 0.d0
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
      f = 0.d0

      ss = b + 2*cc*(T-T0) + 3*d*(T-T0)**2

      dPf = 0.d0
      dTf = ss * m
      dTf = 0.d0
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
      NbCellG, IdCellG, NbFaceG, &
      PermCellG, PermFracG)

    integer, intent(in) :: NbCellG, NbFaceG
    integer, dimension(:), intent(in) :: IdCellG
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
      PermCellG(1,1,i) = 5.d-20
      PermCellG(2,2,i) = 5.d-20
      PermCellG(3,3,i) = 5.d-20
    end do

    allocate(PermFracG(NbFaceG))
    PermFracG(:) = 1.d-11
    
  end subroutine DefModel_SetPerm


#ifdef _THERMIQUE_

  subroutine DefModel_SetCondThermique( &
      NbCellLocal, IdCellLocal, NbFracLocal, &
      CondThermalCellLocal, CondThermalFracLocal)

    integer, intent(in) :: NbCellLocal, NbFracLocal
    integer, dimension(:), intent(in) :: IdCellLocal

    ! ouptuts:
    double precision, dimension(:,:,:), allocatable, intent(inout) :: &
        CondThermalCellLocal
    double precision, dimension(:), allocatable, intent(inout) :: &
        CondThermalFracLocal

    integer :: i

    allocate(CondThermalCellLocal(3,3,NbCellLocal))
    do i=1, NbCellLocal
      CondThermalCellLocal(:,:,i) = 0.d0
      CondThermalCellLocal(1,1,i) = 2.d0
      CondThermalCellLocal(2,2,i) = 2.d0
      CondThermalCellLocal(3,3,i) = 2.d0
    end do

    allocate(CondThermalFracLocal(NbFracLocal))
    CondThermalFracLocal(:) = 2.d0

  end subroutine DefModel_SetCondThermique

  ! Compute Psat(T)
  subroutine DefModel_Psat(T, Psat, dT_PSat)

    double precision, intent(in) :: T
    double precision, intent(out) :: Psat, dT_PSat

    Psat = 1.013E+5*exp(13.7d0 -5120.d0/T)
    dT_PSat = 1.013d+5*5120.d0*exp(13.7d0 -5120.d0/T)/T**2 

  end subroutine DefModel_Psat


  ! Compute Tsat(P)
  subroutine DefModel_Tsat(P, Tsat, dP_Tsat)

    double precision, intent(in) :: P
    double precision, intent(out) :: Tsat, dP_Tsat

    Tsat = -5120.d0 / (13.7d0 - LOG(P/1.013E+5)) 
    dP_Tsat = 1/P/1.013E+5 * -5120.d0/(13.7d0 - LOG(P/1.013E+5))**2 

  end subroutine DefModel_Tsat

#endif

end module DefModel
