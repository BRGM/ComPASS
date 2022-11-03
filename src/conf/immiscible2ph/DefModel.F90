!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp thermal,
! only water in liquid and air in gas phase
! MCP=(1,0,
!      0,1)
!

module DefModel

   use iso_c_binding, only: c_int, c_bool
   use CommonType, only: ModelConfiguration
   use CommonMPI, only: CommonMPI_abort

   implicit none

   ! --------------------------------------------------------------
   !      Definition of the components
   integer, parameter :: NbComp = ComPASS_NUMBER_OF_COMPONENTS
   integer, parameter :: AIR_COMP = ComPASS_AIR_COMPONENT ! =1 (defined in cmake.conf)
   integer, parameter :: WATER_COMP = ComPASS_WATER_COMPONENT ! =2 (defined in cmake.conf) CHECKME: does the water has to be at the end?

   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   !FIXME: Asssume that the latest phase is the liquid phase (wells)
   integer, parameter :: LIQUID_PHASE = ComPASS_LIQUID_PHASE
   integer, parameter :: GAS_PHASE = ComPASS_GAS_PHASE

   ! --------------------------------------------------------------
   !      Definition of the contexts (sets of present phases)
   integer, parameter :: NbContexte = ComPASS_NUMBER_OF_CONTEXTS ! NbContexte = 2**NbPhase - 1

   integer, parameter :: GAS_CONTEXT = ComPASS_GAS_CONTEXT
   integer, parameter :: LIQUID_CONTEXT = ComPASS_LIQUID_CONTEXT
   integer, parameter :: DIPHASIC_CONTEXT = ComPASS_DIPHASIC_CONTEXT

   ! Number of phases that are present in each context
! careful, the lines must coincide with index defined in cmake.conf
   integer, parameter, dimension(NbContexte) :: NbPhasePresente_ctx = (/ &
                                                1, & ! GAS_CONTEXT
                                                1, & ! LIQUID_CONTEXT
                                                2 & ! DIPHASIC_CONTEXT
                                                /)
   ! Index of the phase(s) which is/are present in each context
   ! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
   integer, parameter, dimension(NbPhase, NbContexte) :: NumPhasePresente_ctx = &
                                                         reshape((/ &
                                                                 GAS_PHASE, 0, & ! GAS_CONTEXT
                                                                 LIQUID_PHASE, 0, & ! LIQUID_CONTEXT
                                                                 GAS_PHASE, LIQUID_PHASE & ! DIPHASIC_CONTEXT
                                                                 /), (/NbPhase, NbContexte/))

   ! MCP: the components are potentially present in which phase(s) ?
   integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = RESHAPE((/ &
                    1, 0, & ! AIR_COMP = 1 is present only in gas phase
                    0, 1 & ! WATER_COMP = 2 is present only in liquid phase
                    /), (/NbComp, NbPhase/))

#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

   ! ! ****** Constants derived from model (do not edit) ****** ! !

#include "../common/DefModel_constants.F90"
   ! FIXME: NbIncTotalMax is duplicated l.46 in src/wrappers/NewtonIncrements.h
   integer, parameter :: &
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax, & !< Max number of closure laws
      NbIncTotalMax = NbIncPTCMax + NbPhase           !< Max number of unknowns P (T) C S

   ! ! ****** How to choose primary variables ****** ! !

   ! Used in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
   !
   ! WARNING
   ! Il faut mettre les Sprim en dernier sinon pb dans les liens entre IncPTC et IncTotal
   ! P est forcement primaire et en numero 1
   ! Si T est primaire elle doit etre en numero 2
   ! C est numéroté ensuite (en fonction des phases présentes et des comp dans chaque phase)
   ! S principale

   ! pschoise=3: Gauss method
   !     the matrix psprim and pssecd are defined formally for compile

   integer, parameter :: pschoice = 1
! Global unknowns depending on the context (in the thermal case)
! IF THE COMPONENT CANNOT BE PRESENT IN THE PHASE (due to MCP), IT HAS NO NUMBER.
   ! Careful: the index of unknowns must coincide with the lines of there derivatives in IncPrimSecd.F90
! ic=1 GAS_CONTEXT:       P=1, T=2, Cga=3
! ic=2 LIQUID_CONTEXT:    P=1, T=2, Clw=3
! ic=3 DIPHASIC_CONTEXT:  P=1, T=2, Cga=3, Clw=4, Sprincipal=5        (Sg+Sl=1 is not a closure law, it is forced in the implementation)
!
#ifdef _THERMIQUE_
   integer, parameter, private :: P = 1, T = 2
#else
   integer, parameter, private :: P = 1
#endif

   ! Sum Salpha =1 was already eliminated
   ! (Sg+Sl=1 is not a closure law, it is forced in the implementation)
   integer, parameter, dimension(NbIncTotalPrimMax, NbContexte) :: &
      psprim = RESHAPE((/ &
#ifdef _THERMIQUE_
                       P, T, 0, & ! GAS_CONTEXT=1   "0" in psprim means Ctilde ???
                       P, T, 0, & ! LIQUID_CONTEXT=2   "0" in psprim means Ctilde ???
                       P, T, 5 & ! DIPHASIC_CONTEXT=3  Cga=3, Clw=4, Sprincipal=5
#else
                       P, 2, & ! GAS_CONTEXT=1        Cga=2
                       P, 2, & ! LIQUID_CONTEXT=2     Clw=2
                       P, 4 & ! DIPHASIC_CONTEXT=3   Cga=2, Clw=3, Sprincipal=4
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = RESHAPE((/ &
#ifdef _THERMIQUE_
                       3, 0, 0, 0, & ! GAS_CONTEXT=1       Cga=3
                       3, 0, 0, 0, & ! LIQUID_CONTEXT=2    Clw=3
                       3, 4, 0, 0 & ! DIPHASIC_CONTEXT=3  Cga=3, Clw=4, Sprincipal=5
#else
                       0, 0, 0, 0, & ! GAS_CONTEXT=1
                       0, 0, 0, 0, & ! LIQUID_CONTEXT=2
                       2, 3, 0, 0 & ! DIPHASIC_CONTEXT=3  Cga=2, Clw=3, Sprincipal=4
#endif
                       /), (/NbEqFermetureMax, NbContexte/))

   ! cf. ComPASS.simmulation.AlignmentMethod (in ComPASS.simmulation.__init__.py)
   ! Whatever the alignment method, it is necessary to define aligmat formally to compile
   ! GAS_CONTEXT we sum conservation equations to have non degenerate conservation equation -> pressure block in CPR-AMG
   ! combine row corresponding to equations (must be invertible)
   ! diagonal element of the Jacobian must be non null
   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = RESHAPE((/ &
                        1.d0, 1.d0, 0.d0, & ! GAS_CONTEXT=1            P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        0.d0, 1.d0, 0.d0, & !                          ???: water conservation
                        1.d0, 1.d0, 0.d0, & ! LIQUID_CONTEXT=2         P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        1.d0, 0.d0, 0.d0, & !                          ???: air conservation
                        1.d0, 1.d0, 0.d0, & ! DIPHASIC_CONTEXT=3       P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        1.d0, 0.d0, 0.d0 &  !                          Sg: air conservation
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

#include "../common/DefModel_common.F90"

end module DefModel
