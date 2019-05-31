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

   use iso_c_binding, only: c_int, c_bool
   use CommonType, only: CSR, type_IdNode, ModelConfiguration
   use CommonMPI, only: CommonMPI_abort

   implicit none

   ! ! ****** Model ****** ! !

   integer, parameter :: NbComp = ComPASS_NUMBER_OF_COMPONENTS
   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   integer, parameter :: GAS_PHASE = ComPASS_GAS_PHASE
   integer, parameter :: LIQUID_PHASE = ComPASS_LIQUID_PHASE

   integer, parameter :: NbContexte = ComPASS_NUMBER_OF_CONTEXTS
   integer, parameter :: GAS_CONTEXT = ComPASS_GAS_CONTEXT
   integer, parameter :: LIQUID_CONTEXT = ComPASS_LIQUID_CONTEXT
   integer, parameter :: DIPHASIC_CONTEXT = ComPASS_DIPHASIC_CONTEXT

  ! Number of phases that are present in each context
  integer, parameter, dimension(NbContexte) :: &
      NbPhasePresente_ctx = (/ &
        1, & ! GAS_CONTEXT
        1, & ! LIQUID_CONTEXT
        2  & ! DIPHASIC_CONTEXT
      /)
  ! Numero of the phase(s) which is/are present in each context
  ! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
  integer, parameter, dimension(NbPhase, NbContexte) :: &
      NumPhasePresente_ctx = &
        reshape((/ &
            GAS_PHASE, 0,            & ! GAS_CONTEXT
            LIQUID_PHASE, 0,         & ! LIQUID_CONTEXT
            GAS_PHASE, LIQUID_PHASE  & ! DIPHASIC_CONTEXT
        /), (/NbPhase, NbContexte/))


   ! MCP
   integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = transpose(reshape( &
                      (/1, 1/), (/NbPhase, NbComp/)))

   ! Thermique
#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

#include "../common/DefModel_constants.F90"
  integer, parameter :: &
       NbEqFermetureMax  = NbPhase + NbEqEquilibreMax,   & !< Max number of closure laws
       NbIncTotalMax     = NbIncPTCMax + NbPhase           !< Max number of unknowns P (T) C S

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

! Global unknowns depending on the context (in the thermal case)
! Careful: the index of unknowns must coincide with the lines of there derivatives in IncPrimSecd.F90
! ic=1 GAS_CONTEXT:       P=1, T=2, Cg=3
! ic=2 LIQUID_CONTEXT:    P=1, T=2, Cl=3
! ic=3 DIPHASIC_CONTEXT:  P=1, T=2, Cg=3, Cl=4, Sprincipal=5        (Sg+Sl=1 is not a closure law, it is forced in the implementation)
#ifdef _THERMIQUE_
   integer, parameter, private :: P=1, T=2
#else
   integer, parameter, private :: P=1
#endif
   
   integer, parameter, dimension(NbIncTotalPrimMax, NbContexte) :: &
      psprim = reshape((/ &
#ifdef _THERMIQUE_
                       P, T, & ! ic=1 GAS_CONTEXT
                       P, T, & ! ic=2 LIQUID_CONTEXT
                       P, 5  & ! ic=3 DIPHASIC_CONTEXT     P=1, T=2, Cg=3, Cl=4, Sprincipal=5
#else
                       P, & ! ic=1 GAS_CONTEXT
                       P, & ! ic=2 LIQUID_CONTEXT
                       4  & ! ic=3 DIPHASIC_CONTEXT        P=1, Cg=2, Cl=3, Sprincipal=4
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

  ! Sum Salpha =1 was already eliminated
   ! Sl is deduced from Sg: Sl=1-Sg
   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = reshape((/ &
#ifdef _THERMIQUE_
                       3, 0, 0, & ! ic=1 GAS_CONTEXT           P=1, T=2, Cg=3
                       3, 0, 0, & ! ic=2 LIQUID_CONTEXT        P=1, T=2, Cl=3
                       T, 3, 4  & ! ic=3 DIPHASIC_CONTEXT      P=1, T=2, Cg=3, Cl=4, Sprincipal=5
#else
                       2, 0, 0, & ! ic=1 GAS_CONTEXT        P=1, Cg=2
                       2, 0, 0, & ! ic=2 LIQUID_CONTEXT     P=1, Cl=2
                       0, 0, 0  & ! ic=3 DIPHASIC_CONTEXT is MEANINGLESS here
#endif
                       /), (/NbEqFermetureMax, NbContexte/))

   ! ! ****** Alignment method ****** ! !

  ! Used in module Jacobian.F90
   ! aligmethod=1, manually
   ! The idea is to have postive diagonal using linear combinations
   ! (alternative is to used inverse of block = LC of)
   ! good for LU O (not good for amg)
   ! not used if not preconditionner (but avoid pivoting)
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
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ &   ! NOT USED because aligmethod = 2
#ifdef _THERMIQUE_
                        1.d0, 0.d0, & ! GAS_CONTEXT=1
                        0.d0, 1.d0, &
                        1.d0, 0.d0, & ! LIQUID_CONTEXT=2
                        0.d0, 1.d0, &
                        1.d0, 0.d0, & ! DIPHASIC_CONTEXT=3
                        0.d0, 1.d0 &
#else
                        1.d0, & ! GAS_CONTEXT=1
                        1.d0, & ! LIQUID_CONTEXT=2
                        0.d0  & ! DIPHASIC_CONTEXT=3 is MEANINGLESS here
#endif
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

#include "../common/DefModel_common.F90"

end module DefModel
