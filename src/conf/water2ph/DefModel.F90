!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp thermal, MCP=(1,1)

! 1: Gas
! 2: Water

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
      MCP = transpose(reshape( &
      (/ 1, 1 /), (/NbPhase, NbComp/)))

  ! Phase PHASE_WATER is Liquid; PHASE_GAS is gas
  integer, parameter :: PHASE_GAS = 1
  integer, parameter :: PHASE_WATER = 2

  ! Gravite
  double precision, parameter :: Gravite = 10.d0 !< Gravity constant

  ! CpRoche
  double precision, parameter :: CpRoche = 800.d0 * 2000.d0 !< ???

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
       psprim = reshape( (/ &
       1, 2, & ! ic=1
       1, 2, & ! ic=2
       1, 5  & ! ic=3
       /), (/ NbIncPTCSPrimMax, NbContexte/))

  integer, parameter, dimension( NbIncPTCSecondMax, NbContexte) :: &
       pssecd = reshape( (/ &
       3, 4, 0, & ! ic=1
       3, 4, 0, & ! ic=2
       2, 3, 4  & ! ic=3
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
       aligmat = reshape( (/ &
       1.d0, 0.d0, & ! ic=1
       0.d0, 1.d0, &
       !
       1.d0, 0.d0, & ! ic=2
       0.d0, 1.d0, &
       !
       1.d0, 0.d0, & ! ic=3
       0.d0, 1.d0  &
       /), (/ NbCompThermique, NbCompThermique, NbContexte /))


  ! ! ******* Mesh type ****** ! !
  ! ! * cartesian-quad
  ! ! * hexahedron-quad
  ! ! * tetrahedron-triangle
  ! ! * wedge
  character(len=40) :: &
      MESH_TYPE = "cartesian-quad"


  ! ! ****** Times ****** ! !

  ! One day
  double precision, parameter :: OneSecond = 1.d0
  double precision, parameter :: OneDay = 24.d0 * 3600.d0
  double precision, parameter :: OneMonth = 3.d1 * OneDay
  double precision, parameter :: OneYear = 3.6525d2 * OneDay

  ! time step init and max
  ! FIXME: parameter is removed to assign variable from python
  double precision :: TimeFinal = 30 * OneYear

  double precision :: TimeStepInit = OneDay
  double precision :: TimeStepMax = OneYear

  ! output_frequency for visu
  double precision, parameter :: output_frequency = OneYear


  ! ! ****** Newton iters max and stop condition ****** ! !
  integer, parameter :: NewtonNiterMax = 40
  double precision, parameter :: NewtonTol = 1.d-5


  ! ! ****** ksp linear solver iters max and stop condition ****** ! !
  integer, parameter :: KspNiterMax = 150        ! max nb of iterations
  double precision, parameter :: KspTol = 1.d-6  ! tolerance


  ! ! ****** Obj values used to compute Newton increment ****** ! !

  double precision, parameter ::  &
      NewtonIncreObj_P = 5.d5,   &
      NewtonIncreObj_T = 20.d0,  &
      NewtonIncreObj_C = 1.d0,   &
      NewtonIncreObj_S = 0.2d0


  ! ! ****** Obj values used to compute next time step ****** ! !

  double precision, parameter :: &
      TimeStepObj_P = 5.d5,     &
      TimeStepObj_T = 20.d0,    &
      TimeStepObj_C = 1.d0,     &
      TimeStepObj_S = 0.6d0


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

  ! Fugacity
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

    if(iph==1) then

       f = P
       dPf = 1.d0
       dTf = 0.d0

    else if(iph==2) then

       call DefModel_Psat(T, Psat, dTSat)
       f = Psat
       dPf = 0.d0
       dTf = dTSat
    end if

    dCf(:) = 0.d0
    dSf(:) = 0.d0

   end subroutine f_Fugacity


  ! Densite molaire
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_DensiteMolaire(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    double precision :: Cs, rho0, a, b, a1, a2, b1, b2, c1, c2, cw, dcwp, dcwt
    double precision :: u, R, T0, d, Z, ds, ss, rs, dZ
    double precision :: Psat, dT_Psat

    if(iph==1) then

       u = 0.018016d0
       R = 8.3145d0
       T0 = 273.d0

       a = 1.102d0
       b = -4.461d-4

       Z = 1.d0
       dZ = 0.d0

       f = P*u/(R*T*Z)
       dPf = u/(R*T*Z)
       dTf = - P*u/(R*T*Z)**2 * (R*T*dZ+R*Z)

    else if(iph==2) then

       rs = 0.d0

       call DefModel_Psat(T, Psat, dT_Psat)

       Cs = 35.d0
       rho0 = 780.83795d0
       a = 1.6269192d0
       b = -3.0635410d-3

       a1 = 2.4638d-9
       a2 = 1.1343d-17
       b1 = -1.2171d-11
       b2 = 4.8695d-20
       c1 = 1.8452d-14
       c2 = -5.9978d-23

       ss = (rho0 + a*T + b*T**2)*(1.d0 + 6.51d-4*Cs)
       ds = (a + b*T*2.d0)*(1.d0 + 6.51d-4*Cs)

       cw = (1.d0+5.d-2*rs) &
            *( a1 + a2*P + T*(b1+b2*P) + T**2*(c1+c2*P) )

       dcwp = (1.d0+5.d-2*rs) *( a2 + T*b2 + T**2*c2 )
       dcwt = (1.d0+5.d-2*rs) * ( (b1+b2*P) + T*2.d0*(c1+c2*P) )

       f = ss*( 1.d0 + cw*(P-Psat) )
       dPf = ss*dcwp*(P-Psat) + ss*cw
       dTf = ds*( 1.d0 + cw*(P-Psat) ) + ss*dcwt*(P-Psat) - ss*cw*dT_Psat
    else
       print*, "densite molaire error: f_DensiteMolaire"
    end if

    dCf(:) = 0.d0
    dSf(:) = 0.d0

  end subroutine f_DensiteMolaire


  ! Densite Massique
  ! iph is an identificator for each phase:
  ! PHASE_GAS = 1; PHASE_WATER = 2
  subroutine f_DensiteMassique(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

    ! input
    integer, intent(in) :: iph
    double precision, intent(in) :: P, T, C(NbComp), S(NbPhase)

    ! output
    double precision, intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

    call f_DensiteMolaire(iph,P,T,C,S,f,dPf,dTf,dCf,dSf)

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

    ! outputs
    if(iph==1) then

       f = (0.361d0*T - 10.2d0)*1.d-7
       dPf = 0.d0
       dTf = 0.361*1.d-7

    else if(iph==2) then

       T1 = 300.d0
       Cs = 0.04d0
       a = 1.d0 + Cs*1.34d0 + 6.12d0*Cs**2
       ss = 0.021482*( T-273.d0 - 8.435d0 &
            + sqrt(8078.4d0+(T-273.d0-8.435d0)**2) )
       ss = ss - 1.2
       ds = 0.021482d0 * ( 1.d0 + (T-273.d0-8.435d0) &
            / sqrt(8078.4d0 + (T-273.d0-8.435d0)**2) )

       f= 1.d-3*a/ss
       dPf = 0.d0
       dTf = -a*1.d-3*ds/ss**2
    else
       print*, "viscosite error: f_Viscosite"
    end if

    dCf(:) = 0.d0
    dSf(:) = 0.d0

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

    if(iph==1) then
       f = S(1)**2
       dSf(1) = 2.d0*S(1)
       dSf(2) = 0.d0
    else if(iph==2) then
       f = S(2)**2
       dSf(1) = 0.d0
       dSf(2) = 2.d0*S(2)
    else
       print*, "Perm Rel error"
    end if

  end subroutine f_PermRel


  ! Pressions Capillaires des Phases et leurs derivees
  subroutine f_PressionCapillaire(rt,iph,S,f,DSf)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    integer, intent(in) :: iph
    double precision, intent(in) :: S(NbPhase)

    ! output
    double precision, intent(out) :: f, DSf(NbPhase)

    f = 0.d0
    dSf(:) = 0.d0

  end subroutine f_PressionCapillaire


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

    call f_Enthalpie(iph,P,T,C,S,f,dPf,dTf,dCf, dSf)

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

    double precision :: a, b, cc, d, T0

    if(iph==PHASE_GAS) then
      a = 1990.89d+3
      b = 190.16d+3

      f = a + T*b/100.d0
      dPf = 0.d0
      dTf = b/100.d0

    else if(iph==PHASE_WATER) then
      a = -14.4319d+3
      b = 4.70915d+3
      cc = -4.87534d0
      d = 1.45008d-2
      T0 = 273.d0

      f = a + b*(T-T0) + cc*(T-T0)**2 + d*(T-T0)**3
      dPf = 0.d0
      dTf = b + 2.d0*cc*(T-T0) + 3.d0*d*(T-T0)**2
    end if

    dCf(:) = 0.d0
    dSf(:) = 0.d0

  end subroutine f_Enthalpie

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
      PermCellG(1,1,i) = 1.d-14
      PermCellG(2,2,i) = 1.d-14
      PermCellG(3,3,i) = 1.d-14
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


  allocate(PorositeCell(NbCellG))
  allocate(PorositeFrac(NbFracG))

  PorositeCell(:) = 1.d-1
  PorositeFrac(:) = 4.d-1
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
      CondThermalCellG(1,1,i) = 2.d0
      CondThermalCellG(2,2,i) = 2.d0
      CondThermalCellG(3,3,i) = 2.d0
    end do

    allocate(CondThermalFracG(NbFracG))
    CondThermalFracG(:) = 2.d0

  end subroutine DefModel_SetCondThermique


  SUBROUTINE DefModel_SetThermalSource( &
      NbCell, &
      CellThermalSourceType, &
      NbFrac, &
      FracThermalSourceType, &
      CellThermalSource, &
      FracThermalSource)

    INTEGER, INTENT(IN) :: NbCell
    INTEGER, DIMENSION(:), INTENT(IN) :: CellThermalSourceType
    INTEGER, INTENT(IN) :: NbFrac
    INTEGER, DIMENSION(:), INTENT(IN) :: FracThermalSourceType
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: CellThermalSource
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: FracThermalSource

    ALLOCATE(CellThermalSource(NbCell))
    CellThermalSource = 0.d0

    ALLOCATE(FracThermalSource(NbFrac))
    FracThermalSource = 0.d0
  END SUBROUTINE DefModel_SetThermalSource


  ! Compute Psat(T)
  subroutine DefModel_Psat(T, Psat, dT_PSat)

    double precision, intent(in) :: T
    double precision, intent(out) :: Psat, dT_PSat

    Psat = (T-273.d0)**4.d0 / 1.0d3
    dT_PSat = 4.d0 * (T-273.d0)**3.d0 / 1.0d3

  end subroutine DefModel_Psat


  ! Compute Tsat(P)
  subroutine DefModel_Tsat(P, Tsat, dP_Tsat)

    double precision, intent(in) :: P
    double precision, intent(out) :: Tsat, dP_Tsat

    Tsat = 100.d0 * (P/1.d5)**0.25d0 + 273.d0
    dP_Tsat = 0.25d0 * 1.d-3 * (P/1.d5)**(-0.75d0)

  end subroutine DefModel_Tsat

#endif

end module DefModel
