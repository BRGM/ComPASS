!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 1 phase 2 components thermal

module DefModel

   use iso_c_binding, only: c_int, c_bool
   use CommonType, only: ModelConfiguration
   use CommonMPI, only: CommonMPI_abort

   implicit none

   ! ! ****** Model ****** ! !

   integer, parameter :: NbComp = ComPASS_NUMBER_OF_COMPONENTS
   integer, parameter :: SALT_COMP = ComPASS_SALT_COMPONENT ! (defined in cmake.conf)
   integer, parameter :: WATER_COMP = ComPASS_WATER_COMPONENT ! (defined in cmake.conf) CHECKME: does the water has to be at the end?
   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   integer, parameter :: NbContexte = 1

   ! Number of phases that are present in each context
   integer, parameter, dimension(NbContexte) :: &
      NbPhasePresente_ctx = (/1/)
   ! Numero of the phase(s) which is/are present in each context
   ! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
   integer, parameter, dimension(NbPhase, NbContexte) :: NumPhasePresente_ctx = &
                                                         transpose(reshape((/1/), (/NbContexte, NbPhase/)))

   !FIXME: Asssume that the latest phase is the liquid phase
   integer, parameter :: LIQUID_PHASE = NbPhase
   integer, parameter :: GAS_PHASE = -1 ! FIXME: dummy value for IncPrimSecdFreeFlow

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
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax, & !< Max number of closure laws
      NbIncTotalMax = NbIncPTCMax + NbPhase           !< Max number of unknowns P (T) C S

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
#ifdef _THERMIQUE_
   integer, parameter, private :: P = 1, T = 2, Cs = 3, Cw = 4
#else
   integer, parameter, private :: P = 1, T = 0, Cs = 2, Cw = 3
#endif

   integer, parameter, dimension(NbIncTotalPrimMax, NbContexte) :: &
      psprim = reshape((/ &
#ifdef _THERMIQUE_
                       P, T, Cs & ! single context
#else
                       P, Cs & ! single context
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = reshape((/Cw/), (/NbEqFermetureMax, NbContexte/))

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

   integer, parameter :: aligmethod = 1

   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ &   ! NOT USED because aligmethod = 2
#ifdef _THERMIQUE_
                        1.d0, 1.d0, 0.d0, &
                        0.d0, 0.d0, 1.d0, &
                        0.d0, 1.d0, 0.d0 &
#else
                        1.d0, 1.d0, &
                        0.d0, 1.d0 &
#endif
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

#include "../common/DefModel_common.F90"

end module DefModel
