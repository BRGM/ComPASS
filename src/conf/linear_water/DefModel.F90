!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 1 phase 1 comp thermal, MCP

! 1: Water

module DefModel

   use iso_c_binding, only: c_int, c_bool
   use CommonType, only: CSR, type_IdNode, ModelConfiguration
   use CommonMPI, only: CommonMPI_abort

   implicit none

   integer, parameter :: NbComp = ComPASS_NUMBER_OF_COMPONENTS
   integer, parameter :: WATER_COMP = ComPASS_WATER_COMPONENT
   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   integer, parameter :: NbContexte = 1

   ! Number of phases that are present in each context
   integer, parameter, dimension(NbContexte) :: NbPhasePresente_ctx = (/1/)
   ! Phase(s) that is/are present in each context
   ! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
   integer, parameter, dimension(NbPhase, NbContexte) :: NumPhasePresente_ctx = &
                                                         transpose(reshape((/1/), (/NbContexte, NbPhase/)))

   !FIXME: Asssume that the latest phase is the liquid phase
   integer, parameter :: LIQUID_PHASE = NbPhase
   integer, parameter :: GAS_PHASE = -1 ! FIXME: dummy value for IncPrimSecdFreeFlow

   ! MCP
   integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = transpose(reshape( &
                      (/1/), (/NbPhase, NbComp/)))

   ! Thermique
#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

#include "../common/DefModel_constants.F90"
   ! FIXME: NbIncTotalMax is duplicated l.46 in src/wrappers/NewtonIncrements.h
   integer, parameter :: &
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax, & !< Max number of closure laws
      NbIncTotalMax = NbIncPTCMax + NbPhase           !< Max number of unknowns P (T) C S

   ! ! ****** How to choose primary variables ****** ! !

   ! Used in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
   ! pschoice=2: Glouton method
   !     the matrix psprim and pssecd are defined formally for compile
   ! pschoise=3: Gauss method
   !     the matrix psprim and pssecd are defined formally for compile

   integer, parameter :: pschoice = 1

   ! it is generaly not possible to define flags for the index because the numbering depends on the context (here only one context).
#ifdef _THERMIQUE_
   integer, parameter :: P = 1, T = 2, C = 3, S = 4
#else
   integer, parameter :: P = 1, T = 0, C = 2, S = 3
#endif
   private :: P, T, C, S

   integer, parameter, dimension(NbIncTotalPrimMax, NbContexte) :: &
      psprim = reshape((/ &
#ifdef _THERMIQUE_
                       P, T & ! only one context ic=1
#else
                       P & ! only one context ic=1
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

   ! Sum Salpha =1 was already eliminated (only one phase in this physic)
   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = reshape((/C/), (/NbEqFermetureMax, NbContexte/))

   ! ! ****** Alignment method ****** ! !
   ! Used in module Jacobian.F90
   !
   ! aligmethod=1, manually using aligmat
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
   ! aligmethod=2, inverse diagonal
   !     it is necessary to define aligmat formally for compile

   integer, parameter :: aligmethod = 2

   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ & ! NOT USED because aligmethod = 2
#ifdef _THERMIQUE_
                        1.d0, 0.d0, & ! only one context ic=1    P: sum(component conservation)
                        0.d0, 1.d0 & !                          T: energy conservation
#else
                        1.d0 &        ! only one context ic=1    P: sum(component conservation)
#endif
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

#include "../common/DefModel_common.F90"

end module DefModel
