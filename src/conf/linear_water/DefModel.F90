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
   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   integer, parameter :: NbContexte = 1

  ! Number of phases that are present in each context
  integer, parameter, dimension(NbContexte) :: NbPhasePresente_ctx = (/ 1 /)
  ! Phase(s) that is/are present in each context
  ! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
  integer, parameter, dimension(NbPhase, NbContexte) :: NumPhasePresente_ctx = &
    transpose(reshape((/ 1 /), (/NbContexte, NbPhase/)))

   !FIXME: Asssume that the latest phase is the liquid phase
   integer, parameter :: LIQUID_PHASE = NbPhase

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

   include '../common/DefModel_constants.F90'

   ! ! ****** How to choose primary variables ****** ! !

   ! Used in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
   ! pschoice=2: Glouton method
   !     the matrix psprim and pssecd are defined formally for compile
   ! pschoise=3: Gauss method
   !     the matrix psprim and pssecd are defined formally for compile

   integer, parameter :: pschoice = 1

#ifdef _THERMIQUE_
   integer, parameter :: P=1, T=2, C=3, Sg=4, Sl=5
#else
   integer, parameter :: P=1, T=0, C=2, Sg=3, Sl=4
#endif
   private :: P, T, C, Sg, Sl
   
   integer, parameter, dimension(NbIncTotalPrimMax, NbContexte) :: &
      psprim = reshape((/ &
#ifdef _THERMIQUE_
                       P, T & ! only one context ic=1 
#else
                       P    & ! only one context ic=1
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

   ! Sum Salpha =1 was already eliminated
   ! Sl is deduced from Sg: Sl=1-Sg
   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = reshape((/ &
#ifdef _THERMIQUE_
                       C, Sg, 0 & ! only one context ic=1
#else
                       C, Sg, 0 & ! only one context ic=1
#endif
                       /), (/NbEqFermetureMax, NbContexte/))

   ! ! ****** Alignment method ****** ! !
   ! Used in module Jacobian.F90
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
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = reshape((/ &
#ifdef _THERMIQUE_
                        1.d0, 0.d0, & ! only one context ic=1
                        0.d0, 1.d0  &
#else
                        1.d0 &        ! only one context ic=1
#endif
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

    include '../common/DefModel_common.F90'

end module DefModel
