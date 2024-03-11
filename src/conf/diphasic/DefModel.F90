!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 2 comp thermal, MCP=(1,1,1,1), freeflow interface
!
! Gas and Liquid
!
! 1: Air
! 2: H2O

module DefModel

   use iso_c_binding, only: c_int, c_bool
   use CommonType, only: ModelConfiguration
   use CommonMPI, only: CommonMPI_abort

   implicit none

   ! ! ****** Model ****** ! !

   integer, parameter :: NbComp = ComPASS_NUMBER_OF_COMPONENTS
   integer, parameter :: AIR_COMP = ComPASS_AIR_COMPONENT ! =1 (defined in cmake.conf)
   integer, parameter :: WATER_COMP = ComPASS_WATER_COMPONENT ! =2 (defined in cmake.conf) CHECKME: does the water has to be at the end?

   integer, parameter :: NbPhase = ComPASS_NUMBER_OF_PHASES
   !FIXME: Asssume that the latest phase is the liquid phase (wells)
   integer, parameter :: LIQUID_PHASE = ComPASS_LIQUID_PHASE
   integer, parameter :: GAS_PHASE = ComPASS_GAS_PHASE

   ! --------------------------------------------------------------
   !      Definition of the contexts (sets of present phases)
   integer, parameter :: NbContexte = ComPASS_NUMBER_OF_CONTEXTS ! 2**NbPhase - 1 contexts for the reservoir nodes + 3 for the freeflow contexts

   integer, parameter :: GAS_CONTEXT = ComPASS_GAS_CONTEXT
   integer, parameter :: LIQUID_CONTEXT = ComPASS_LIQUID_CONTEXT
   integer, parameter :: DIPHASIC_CONTEXT = ComPASS_DIPHASIC_CONTEXT
   ! Context with NO Liquid outflow in the FreeFlow BC
   integer, parameter :: GAS_FF_NO_LIQ_OUTFLOW_CONTEXT = ComPASS_GAS_FF_NO_LIQ_OUTFLOW_CONTEXT  ! GAS_CONTEXT with NO_LIQUID_OUTFLOW
   integer, parameter :: DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT = ComPASS_DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT ! DIPHASIC_CONTEXT with NO_LIQUID_OUTFLOW
   ! Context WITH Liquid outflow in the FreeFlow BC
   integer, parameter :: DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT = ComPASS_DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT ! DIPHASIC_CONTEXT with LIQUID_OUTFLOW_CONTEXT

! Number of phases that are present in each context
! careful, the lignes must coincide with index defined in cmake.conf
   integer, parameter, dimension(NbContexte) :: NbPhasePresente_ctx = (/ &
                                                1, & ! GAS_CONTEXT
                                                1, & ! LIQUID_CONTEXT
                                                2, & ! DIPHASIC_CONTEXT
                                                1, & ! GAS_FF_NO_LIQ_OUTFLOW_CONTEXT
                                                2, & ! DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT
                                                2 & ! DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
                                                /)
! Index of the phase(s) which is/are present in each context
! FIXME: NB: we could deduce NbPhasePresente_ctx from this array
   integer, parameter, dimension(NbPhase, NbContexte) :: NumPhasePresente_ctx = &
                                                         reshape((/ &
                                                                 GAS_PHASE, 0, & ! GAS_CONTEXT
                                                                 LIQUID_PHASE, 0, & ! LIQUID_CONTEXT
                                                                 GAS_PHASE, LIQUID_PHASE, & ! DIPHASIC_CONTEXT
                                                                 GAS_PHASE, 0, & ! GAS_FF_NO_LIQ_OUTFLOW_CONTEXT
                                                                 GAS_PHASE, LIQUID_PHASE, & ! DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT
                                                                 GAS_PHASE, LIQUID_PHASE & ! DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT
                                                                 /), (/NbPhase, NbContexte/))

   ! MCP: the components are potentially present in which phase(s) ?
   integer, parameter, dimension(NbComp, NbPhase) :: &
      MCP = RESHAPE((/ &
                    1, 1, & ! ! gas phase : both components
                    1, 1 & ! liquid phase : both components
                    /), (/NbComp, NbPhase/))

   ! Thermique
#ifdef _THERMIQUE_
   integer, parameter :: IndThermique = 1
#else
   integer, parameter :: IndThermique = 0
#endif

   ! ! ****** Constants derived from model (do not edit) ****** ! !

#include "../common/DefModel_constants.F90"
   ! FIXME: NbIncTotalMax is duplicated l.46 in src/wrappers/NewtonIncrements.h
   integer, parameter :: &
      NbEqFermetureMax = NbPhase + NbEqEquilibreMax + NbPhase - 1, & !< Max number of closure laws (+NbPhase-1 for the FreeFlow BC)
      NbIncTotalMax = NbIncPTCMax + 2*NbPhase                   !< Max number of unknowns P (T) C S FIXME: add q(NbPhase)

   ! ! ****** How to choose primary variables ****** ! !

   ! Used in module LoisthermoHydro.F90

   ! pschoice=1: manually
   !     it is necessary to give PTCS Prim and PTC Secd for each context: psprim
   !
   ! WARNING
   ! Il faut mettre les Sprim en dernier sinon pb dans les liens entre IncPTC et IncTotal
   ! P est forcement primaire et en numero 1   (not in this physic)
   ! Si T est primaire elle doit etre en numero 2
   ! C est numéroté ensuite (en fonction des phases présentes et des comp dans chaque phase)
   ! S principale
   ! freeflow flowrate of every phases
   ! Ctilde
   ! pschoise=3: Gauss method
   !     the matrix psprim and pssecd are defined formally for compile

   integer, parameter :: pschoice = 1

! Global unknowns depending on the context (in the thermal case)
! Careful: the index of unknowns must coincide with the lines of there derivatives in IncPrimSecd.F90
! ic=1 GAS_CONTEXT:       P=1, T=2, Cga=3, Cgw=4
! ic=2 LIQUID_CONTEXT:    P=1, T=2, Cla=3, Clw=4
! ic=3 DIPHASIC_CONTEXT:  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7        (Sg+Sl=1 is not a closure law, it is forced in the implementation)
!
!            --------------------------------------------------
! Freeflow dof, global unknowns depending on the context (in the thermal case)
!!!!!!!!!!!!!!!!
! ic=4 GAS_FF_NO_LIQ_OUTFLOW_CONTEXT: P=1, T=2, Cga=3, Cgw=4, freeflow_gas_flowrate=5
! ic=5 DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT: P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7, freeflow_gas_flowrate=8, Ctilde (there is no in this physical example)
! ic=6 DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT:   (the context is diphasic even if IndOutflow = 1, then Pc = 0 and Sg = 0)
!         P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, freeflow_gas_flowrate=8, freeflow_liquid_flowrate=9, Ctilde (there is no in this physical example)
!                                               7 is reserved for Sprincipal (in IncPrimSecd.F90) but is not used as S^g=0 is in hard
!
! Remark: Sg+Sl=1 is not a closure law, as well as Pphase=Pref+Pc(phase) and q^l_atm * S^g = 0, it is forced in the implementation

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
                       ! reservoir dof
                       P, T, 3, & ! GAS_CONTEXT=1       P=1, T=2, Cga=3, Cgw=4
                       P, T, 3, & ! LIQUID_CONTEXT=2    P=1, T=2, Cla=3, Clw=4
                       P, T, 7, & ! DIPHASIC_CONTEXT=3  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7
                       ! FF dof : NO_LIQUID_OUTFLOW
                       T, 3, 5, & ! GAS_FF_NO_LIQ_OUTFLOW_CONTEXT=4       P=1, T=2, Cga=3, Cgw=4, freeflow_gas_flowrate=5
                       3, 7, 8, & ! DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT=5  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7, freeflow_gas_flowrate=8, freeflow_liquid_flowrate=9 (=0.0 so eliminated)
                       ! FF dof : LIQUID_OUTFLOW_CONTEXT
                       T, 8, 9 & ! DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT=6  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, freeflow_gas_flowrate=8, freeflow_liquid_flowrate=9
#else
#error Freeflow BC is not implemented for isothermal case.
                       !    P, 2, & ! GAS_CONTEXT=1        Cga=2, Cgw=3
                       !    P, 3, & ! LIQUID_CONTEXT=2     Cla=2, Clw=3
                       !    P, 6 & ! DIPHASIC_CONTEXT=3   Cga=2, Cgw=3, Cla=4, Clw=5, Sprincipal=6
#endif
                       /), (/NbIncTotalPrimMax, NbContexte/))

   integer, parameter, dimension(NbEqFermetureMax, NbContexte) :: &
      pssecd = RESHAPE((/ &
#ifdef _THERMIQUE_
                       ! reservoir dof
                       4, 0, 0, 0, 0, & ! GAS_CONTEXT=1       P=1, T=2, Cga=3, Cgw=4
                       4, 0, 0, 0, 0, & ! LIQUID_CONTEXT=2    P=1, T=2, Cla=3, Clw=4
                       3, 4, 5, 6, 0, & ! DIPHASIC_CONTEXT=3  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7
                       ! FF dof : NO_LIQUID_OUTFLOW
                       P, 4, 0, 0, 0, & ! GAS_FF_NO_LIQ_OUTFLOW_CONTEXT=4       P=1, T=2, Cga=3, Cgw=4, freeflow_gas_flowrate=5
                       P, T, 4, 5, 6, & ! DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT=5  P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, Sprincipal=7, freeflow_gas_flowrate=8, freeflow_liquid_flowrate=9 (=0.0 so eliminated)
                       ! FF dof : LIQUID_OUTFLOW_CONTEXT
                       P, 3, 4, 5, 6 & ! DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT=6    P=1, T=2, Cga=3, Cgw=4, Cla=5, Clw=6, freeflow_gas_flowrate=8, freeflow_liquid_flowrate=9
#else
#error Freeflow BC is not implemented for isothermal case.
                       !    3, 0, 0, 0, 0, & ! GAS_CONTEXT=1       Cga=2, Cgw=3
                       !    2, 0, 0, 0, 0, & ! LIQUID_CONTEXT=2    Cla=2, Clw=3
                       !    2, 3, 4, 5, 0 & ! DIPHASIC_CONTEXT=3  Cga=2, Cgw=3, Cla=4, Clw=5, Sprincipal=6
#endif
                       /), (/NbEqFermetureMax, NbContexte/))

   ! cf. ComPASS.simmulation.AlignmentMethod (in ComPASS.simmulation.__init__.py)
   ! Whatever the alignment method, it is necessary to define aligmat formally to compile
   ! ic=1 ! we sum conservation equations to have non degenerate conservation equation -> pressure block in CPR-AMG
   ! combine row corresponding to equations (must be invertible)
   ! diagonal element of the Jacobian must be non null
   double precision, parameter, &
      dimension(NbCompThermique, NbCompThermique, NbContexte) :: &
      aligmat = RESHAPE((/ &
                        !  reservoir dof
                        1.d0, 1.d0, 0.d0, & ! GAS_CONTEXT=1            P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        1.d0, 0.d0, 0.d0, & !                          Cga: air conservation
                        1.d0, 1.d0, 0.d0, & ! LIQUID_CONTEXT=2         P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        0.d0, 1.d0, 0.d0, & !                          Clw: water conservation
                        1.d0, 1.d0, 0.d0, & ! DIPHASIC_CONTEXT=3       P: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & !                          T: energy conservation
                        0.d0, 1.d0, 0.d0, & !                          Sg: water conservation
                        ! FF dof : NO_LIQUID_OUTFLOW
                        0.d0, 0.d0, 1.d0, & ! GAS_FF_NO_LIQ_OUTFLOW_CONTEXT=4       T: energy conservation
                        1.d0, 0.d0, 0.d0, & !                                       Cga: air conservation
                        1.d0, 1.d0, 0.d0, & !                                       gas flowrate: sum(component conservation)
                        0.d0, 0.d0, 1.d0, & ! DIPHASIC_FF_NO_LIQ_OUTFLOW_CONTEXT=5    T: energy conservation
                        1.d0, 0.d0, 0.d0, & !                                         Cga: air conservation
                        1.d0, 1.d0, 0.d0, & !                                         gas flowrate: sum(component conservation)
                        ! FF dof : LIQUID_OUTFLOW_CONTEXT
                        0.d0, 0.d0, 1.d0, & ! DIPHASIC_FF_LIQ_OUTFLOW_CONTEXT=6   T: energy conservation
                        1.d0, 1.d0, 0.d0, & !                                     gas flowrate: sum(component conservation)
                        0.d0, 1.d0, 0.d0 & !                                     liq flowrate: water conservation
                        /), (/NbCompThermique, NbCompThermique, NbContexte/))

#include "../common/DefModel_common.F90"

end module DefModel
