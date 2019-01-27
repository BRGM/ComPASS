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

  integer, parameter :: GAS_CONTEXT = 1
  integer, parameter :: LIQUID_CONTEXT = 2
  integer, parameter :: DIPHASIC_CONTEXT = 3

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

  ! Used in module LoisthermoHydro.F90

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

#ifdef _THERMIQUE_
   integer, parameter, private :: P=1, T=2
#else
   integer, parameter, private :: P=1
#endif


  ! Sum Salpha =1 was already eliminated
  integer, parameter, dimension( NbIncPTCSPrimMax, NbContexte) :: &
    psprim = RESHAPE( (/ &
#ifdef _THERMIQUE_
      P, T, 3, & ! ic=1 GAS_CONTEXT       Cga=3, Cgw=4
      P, T, 4, & ! ic=2 LIQUID_CONTEXT    Cla=3, Clw=4
      P, T, 7  & ! ic=3 DIPHASIC_CONTEXT  Cga=3, Cgw=4, Cla=5, Clw=6, Sg=7, Sl=8
#else
      P, 2, & ! ic=1 GAS_CONTEXT        Cga=2, Cgw=3
      P, 3, & ! ic=2 LIQUID_CONTEXT     Cla=2, Clw=3
      P, 6  & ! ic=3 DIPHASIC_CONTEXT   Cga=2, Cgw=3, Cla=4, Clw=5, Sg=6, Sl=7
#endif
      /), (/ NbIncPTCSPrimMax, NbContexte /))

  integer, parameter, dimension( NbIncPTCSecondMax, NbContexte) :: &
       pssecd = RESHAPE( (/ &
#ifdef _THERMIQUE_
       4, 0, 0, 0, & ! ic=1 GAS_CONTEXT       Cga=3, Cgw=4
       3, 0, 0, 0, & ! ic=2 LIQUID_CONTEXT    Cla=3, Clw=4
       3, 4, 5, 6  & ! ic=3 DIPHASIC_CONTEXT  Cga=3, Cgw=4, Cla=5, Clw=6, Sg=7, Sl=8
#else
       3, 0, 0, 0, & ! ic=1 GAS_CONTEXT       Cga=2, Cgw=3
       2, 0, 0, 0, & ! ic=2 LIQUID_CONTEXT    Cla=2, Clw=3
       2, 3, 4, 5  & ! ic=3 DIPHASIC_CONTEXT  Cga=2, Cgw=3, Cla=4, Clw=5, Sg=6, Sl=7
#endif
       /), (/ NbIncPTCSecondMax, NbContexte/))

  ! ! ****** Alignment method ****** ! !

  ! Used in module Jacobian.F90
  ! The idea is to have postive diagonal using linear combinations
  ! (alternative is to used inverse of block = LC of)
  ! good for LU O (pas bonne pour amg)
  ! not used if not preconditionner (but avoid pivoting)
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

    ! ic=1 ! we sum conservation equations to have non degenerate conservation equation -> pressure block in CPR-AMG
    ! combine row corresponding to equations (must be invertible)
    ! diagonal element of the Jacobian must be non null
  double precision, parameter, &
       dimension( NbCompThermique, NbCompThermique, NbContexte) :: &
       aligmat = RESHAPE( (/ &
       1.d0, 1.d0, 0.d0, & ! ic=1 GAS_CONTEXT
       0.d0, 0.d0, 1.d0, &
       1.d0, 0.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=2 LIQUID_CONTEXT
       0.d0, 0.d0, 1.d0, &
       0.d0, 1.d0, 0.d0, &
       1.d0, 1.d0, 0.d0, & ! ic=3 DIPHASIC_CONTEXT
       0.d0, 0.d0, 1.d0, &
       0.d0, 1.d0, 0.d0  &
       /), (/ NbCompThermique, NbCompThermique, NbContexte /) )

end module DefModel
