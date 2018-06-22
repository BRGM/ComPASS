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
                      (/1, 1/), (/NbPhase, NbComp/)))


   !FIXME: this is used for wells which are monophasic
   integer, parameter :: LIQUID_PHASE = 1

   ! Thermique
#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

   ! ! ****** Constants derived from model (do not edit) ****** ! !

   ! Nombre Max d'eq d'equilibre
   !               d'eq de fermeture thermodynamique
   !               d'inc P (T) C
   !               d'inc P (T) C primaires
   integer, parameter :: &
      NbEqEquilibreMax = NbComp*(NbPhase - 1), & !< Max number of balance equations
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax, & !< Max number of closure laws
      NbIncPTCMax = 1 + IndThermique + sum(MCP), &
      NbIncPTCSecondMax = NbEqFermetureMax, &
      NbIncPTCSMax = NbIncPTCMax + NbPhase, &
      NbIncPTCSPrimMax = NbComp + IndThermique, &
      NbCompThermique = NbComp + IndThermique

   ! ! ****** How to choose primary variables ****** ! !

   ! Used in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
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
   integer, parameter :: P=1, T=2, Cw=3, Ct=4, Sl=5
#else
   integer, parameter :: P=1, T=-1, Cw=2, Ct=3, Sl=4
#endif
   private :: P, T, Cw, Ct, Sl
   
   integer, parameter, dimension(NbIncPTCSPrimMax, NbContexte) :: &
      psprim = reshape((/ &
#ifdef _THERMIQUE_
                       P, T, Ct & ! only one context
#else
                       P, Ct & ! only one context
#endif
                       /), (/NbIncPTCSPrimMax, NbContexte/))

  ! Sum Salpha =1 was already eliminated
   ! Sl is deduced from Sg: Sl=1-Sg
   integer, parameter, dimension(NbIncPTCSecondMax, NbContexte) :: &
      pssecd = reshape((/ &
                       Cw & ! only one context
                       /), (/NbIncPTCSecondMax, NbContexte/))

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
   ! aligmethod=2, inverse diagonal
   !     it is necessary to define aligmat formally for compile

   integer, parameter :: aligmethod = 1

   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ &
#ifdef _THERMIQUE_
            1.d0, 1.d0, 0.d0,  & ! only one context
            0.d0, 0.d0, 1.d0,  &
            0.d0, 1.d0, 0.d0   &
#else
            1.d0, 0.d0,  & ! only one context
            0.d0, 1.d0   &
#endif
/), (/NbCompThermique, NbCompThermique, NbContexte/))

end module DefModel
