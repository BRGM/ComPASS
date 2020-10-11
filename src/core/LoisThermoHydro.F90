!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LoisThermoHydro

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use CommonMPI, only: commRank, CommonMPI_abort
   use Thermodynamics, only: &
#ifdef _THERMIQUE_
      f_EnergieInterne, f_Enthalpie, f_SpecificEnthalpy, &
#endif
      f_Viscosite, f_DensiteMolaire, f_PressionCapillaire, f_DensiteMassique
   use RelativePermeabilities, only: f_PermRel
   use DefModel, only: &
#ifdef _WIP_FREEFLOW_STRUCTURES_
      WATER_COMP, &
#endif
      NbPhase, NbComp, IndThermique, LIQUID_PHASE, MCP, &
      NbEqEquilibreMax, NbIncPTCMax, NbIncTotalPrimMax, &
      NbIncTotalMax, NbEqFermetureMax, &
      NbPhasePresente_ctx, NumPhasePresente_ctx
   use NumbyContext, only: &
      NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, NumIncComp2NumIncPTC_ctx, &
      NbEqEquilibre_ctx, NbEqFermeture_ctx, NumCompCtilde_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
      NbCompCtilde_ctx, NbIncPTC_ctx, NbIncTotalPrim_ctx
   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, &
      IncNode, IncCell, IncFrac, &
      NbCellLocal_Ncpus, NbNodeLocal_Ncpus, NbFracLocal_Ncpus
   use IncCVWells, only: &
      PerfoWellInj, DataWellInjLocal, NodeByWellInjLocal, NbWellInjLocal_Ncpus
   use IncPrimSecd, only: &
      dXssurdXpCell, dXssurdXpNode, dXssurdXpFrac, &
      SmdXsCell, SmdXsNode, SmdXsFrac, SmFNode, SmFCell, SmFFrac, &
      NumIncTotalPrimCell, NumIncTotalPrimNode, NumIncTotalPrimFrac, &
      NumIncTotalSecondCell, NumIncTotalSecondNode, NumIncTotalSecondFrac
   use MeshSchema, only: &
#ifdef _WIP_FREEFLOW_STRUCTURES_
      IdFFNodeLocal, &
#endif
      NodeDatabyWellInjLocal, NbWellProdLocal_Ncpus, &
      CellRocktypeLocal, FracRocktypeLocal, NodeRocktypeLocal, &
      PhaseDOFFamilyArray, MeshSchema_allocate_PhaseDOFFamilyArray, MeshSchema_free_PhaseDOFFamilyArray, &
      CompPhaseDOFFamilyArray, MeshSchema_allocate_CompPhaseDOFFamilyArray, MeshSchema_free_CompPhaseDOFFamilyArray
#ifdef _WIP_FREEFLOW_STRUCTURES_
   use Physics, only: atm_comp, Hm, HT, atm_temperature, atm_flux_radiation, &
                      soil_emissivity, Stephan_Boltzmann_cst, atm_pressure
#endif

   implicit none

   ! densite massique
   ! Rq important: it contains values for all phases, not only present phases
   double precision, allocatable, dimension(:, :), protected :: &
      DensiteMassiqueCell, &
      DensiteMassiqueFrac, &
      DensiteMassiqueNode
   double precision, allocatable, dimension(:, :, :), protected :: &
      divDensiteMassiqueCell, &
      divDensiteMassiqueFrac, &
      divDensiteMassiqueNode
   double precision, allocatable, dimension(:, :), protected :: &
      SmDensiteMassiqueCell, &
      SmDensiteMassiqueFrac, &
      SmDensiteMassiqueNode

   ! pression
   double precision, allocatable, dimension(:, :), protected :: &
      divPressionCell, &
      divPressionFrac, &
      divPressionNode
   double precision, allocatable, dimension(:), protected :: &
      SmPressionCell, &
      SmPressionFrac, &
      SmPressionNode

   ! pression capillaire
   double precision, allocatable, dimension(:, :), protected :: &
      PressionCapCell, &
      PressionCapFrac, &
      PressionCapNode
   double precision, allocatable, dimension(:, :, :), protected :: &
      divPressionCapCell, &
      divPressionCapFrac, &
      divPressionCapNode

   ! Saturation
   double precision, allocatable, dimension(:, :, :), protected :: &
      divSaturationCell, &
      divSaturationFrac, &
      divSaturationNode

#ifdef _WIP_FREEFLOW_STRUCTURES_
   ! Freeflow phase molar flowrates
   double precision, allocatable, dimension(:, :, :), protected :: &
      divFreeFlowMolarFlowrateNode, &
      FreeFlowMolarFlowrateCompNode, &
      SmFreeFlowMolarFlowrateCompNode, &
      FreeFlowHmCompNode, &
      SmFreeFlowHmCompNode
   double precision, allocatable, dimension(:, :), protected :: &
      SmFreeFlowMolarFlowrateNode
   double precision, allocatable, dimension(:, :, :, :), protected :: &
      divFreeFlowMolarFlowrateCompNode, &
      divFreeFlowHmCompNode
   ! Thermal vectors
   double precision, allocatable, dimension(:), protected :: &
      FreeFlowHTTemperatureNetRadiationNode, &
      SmFreeFlowHTTemperatureNetRadiationNode
   double precision, allocatable, dimension(:, :), protected :: &
      FreeFlowMolarFlowrateEnthalpieNode, &
      SmFreeFlowMolarFlowrateEnthalpieNode, &
      divFreeFlowHTTemperatureNetRadiationNode, &
      AtmEnthalpieNode
   double precision, allocatable, dimension(:, :, :), protected :: &
      divFreeFlowMolarFlowrateEnthalpieNode
#endif

   ! DensiteMolaire*Kr/Viscosite*Comp
   double precision, allocatable, dimension(:, :, :), protected :: &
      DensiteMolaireKrViscoCompCell, &
      DensiteMolaireKrViscoCompFrac, &
      DensiteMolaireKrViscoCompNode
   double precision, allocatable, dimension(:, :, :, :), protected :: &
      divDensiteMolaireKrViscoCompCell, &
      divDensiteMolaireKrViscoCompFrac, &
      divDensiteMolaireKrViscoCompNode
   double precision, allocatable, dimension(:, :, :), protected :: &
      SmDensiteMolaireKrViscoCompCell, &
      SmDensiteMolaireKrViscoCompFrac, &
      SmDensiteMolaireKrViscoCompNode

   ! DensiteMolaire*Kr/Viscosite*Comp for wells (injection and production)
   double precision, allocatable, dimension(:, :), protected :: &
      DensiteMolaireKrViscoCompWellInj
   double precision, allocatable, dimension(:, :), protected :: &
      divDensiteMolaireKrViscoCompWellInj

   ! DensiteMolaire * Sat * Comp
   type(CompPhaseDOFFamilyArray), target, protected :: DensiteMolaireSatComp
   double precision, allocatable, dimension(:, :, :, :), protected :: &
      divDensiteMolaireSatCompCell, &
      divDensiteMolaireSatCompFrac, &
      divDensiteMolaireSatCompNode
   type(CompPhaseDOFFamilyArray), target, protected :: SmDensiteMolaireSatComp

   ! temperature
   double precision, allocatable, dimension(:, :), protected :: &
      divTemperatureCell, &
      divTemperatureFrac, &
      divTemperatureNode
   double precision, allocatable, dimension(:), protected :: &
      SmTemperatureCell, &
      SmTemperatureFrac, &
      SmTemperatureNode

   ! DensiteMolaire * PermRel / Viscosite * Enthalpie
   double precision, allocatable, dimension(:, :), protected :: &
      DensiteMolaireKrViscoEnthalpieCell, &
      DensiteMolaireKrViscoEnthalpieFrac, &
      DensiteMolaireKrViscoEnthalpieNode
   double precision, allocatable, dimension(:, :, :), protected :: &
      divDensiteMolaireKrViscoEnthalpieCell, &
      divDensiteMolaireKrViscoEnthalpieFrac, &
      divDensiteMolaireKrViscoEnthalpieNode
   double precision, allocatable, dimension(:, :), protected :: &
      SmDensiteMolaireKrViscoEnthalpieCell, &
      SmDensiteMolaireKrViscoEnthalpieFrac, &
      SmDensiteMolaireKrViscoEnthalpieNode

   ! DensiteMolaire * PermRel / Viscosite * Enthalpie for injection wells
   double precision, allocatable, dimension(:), protected :: &
      DensiteMolaireKrViscoEnthalpieWellInj
   double precision, allocatable, dimension(:), protected :: &
      divDensiteMolaireKrViscoEnthalpieWellInj

   ! densitemolaire * energieinterne * Saturation
   type(PhaseDOFFamilyArray), target, protected :: DensiteMolaireEnergieInterneSat
   double precision, allocatable, dimension(:, :, :), protected :: &
      divDensiteMolaireEnergieInterneSatCell, &
      divDensiteMolaireEnergieInterneSatFrac, &
      divDensiteMolaireEnergieInterneSatNode
   double precision, allocatable, dimension(:, :), protected :: &
      SmDensiteMolaireEnergieInterneSatCell, &
      SmDensiteMolaireEnergieInterneSatFrac, &
      SmDensiteMolaireEnergieInterneSatNode

   ! work arrays - with wa_ prefix
   real(c_double), allocatable, target, dimension(:, :) :: wa_UnsurViscosite
   real(c_double), allocatable, target, dimension(:, :, :) :: wa_divUnsurViscosite
   real(c_double), allocatable, target, dimension(:, :) :: wa_SmUnsurViscosite   ! tmp values to simpfy notations of numerotation

   ! tmp values to simpfy notations of numerotation
   ! ex. NbPhasePresente = NbPhasePresente_ctx(inc%ic)
   type ContextInfo
      integer :: &
         NbPhasePresente, NbCompCtilde, &
         NbEqFermeture, NbEqEquilibre, &
         NbIncPTC, NbIncPTCPrim, &
         NbIncTotalPrim, &
         !
         NumPhasePresente(NbPhase), &
         NumCompCtilde(NbComp), &
         NumCompEqEquilibre(NbEqEquilibreMax), &
         NumIncPTC2NumIncComp_comp(NbIncPTCMax), &
         NumIncPTC2NumIncComp_phase(NbIncPTCMax), &
         !
         Num2PhasesEqEquilibre(2, NbEqEquilibreMax), &
         NumIncComp2NumIncPTC(NbComp, NbPhase)
   end type ContextInfo

   public :: &
      LoisThermoHydro_allocate, &
      LoisThermoHydro_free, &
      LoisThermoHydro_compute, &
      LoisThermoHydro_divP_wellinj, &
      LoisThermoHydro_divPrim_nodes

   private :: &
      LoisThermoHydro_divPrim_cv, & ! main function for prim divs for each control volume (cv)
#ifdef _WIP_FREEFLOW_STRUCTURES_
      LoisThermoHydro_divPrim_FreeFlow_cv, & ! prim divs for Molar flowrates in Freeflow dof
      LoisThermoHydro_FreeFlowMolarFlowrateComp_cv, & ! FreeFlowMolarFlowrate * Comp
      LoisThermoHydro_FreeFlowHmComp_cv, & ! Hm * (Ci - Ci_atm)
#ifdef _THERMIQUE_
      LoisThermoHydro_FreeFlowMolarFlowrateEnthalpie_cv, & ! FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)  ; FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)
      LoisThermoHydro_FreeFlowHTTemperatureNetRadiation_cv, & ! HT * (T - atm_temperature) - net Radiation (which is a factor of T**4)
      LoisThermoHydro_AtmEnthalpie_cv, & ! SpecificEnthalpy(water, gas) of the far field atmosphere ; Enthalpie(liquid) of the far field atmosphere
#endif
#endif
      LoisThermoHydro_init_cv, & ! init infos according to ic (context) for each control volume (cv)
      !
      LoisThermoHydro_fill_gradient_dfdX, & ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
      LoisThermoHydro_dfdX_ps, & ! fill dfdX_prim/dfdX_secd with the derivatives w.r.t. the primary/secondary unknowns
      LoisThermoHydro_DensiteMolaire_cv, & ! prim divs: densitemolaire
      LoisThermoHydro_Viscosite_cv, & !            1/viscosite
      LoisThermoHydro_PermRel_cv, & !            Permrel
      LoisThermoHydro_Inc_cv, & !            called with Pression and Temperature
      LoisThermoHydro_PressionCapillaire_cv, & !            Pressioncapillaire
      LoisThermoHydro_Saturation_cv, & !            Saturation
      !
      LoisThermoHydro_DensiteMassique_cv, & !          densitemassique
      LoisThermoHydro_DensiteMolaireKrViscoComp_cv, & !  densitemolaire * Permrel / viscosite * Comp
      LoisThermoHydro_DensiteMolaireSatComp_cv, & !  densitemolaire * Saturation * Comp
      LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv, & !  densitemolaire * Permrel / viscosite * Enthalpie
      LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv  !  densitemolaire * energieinterne * Saturation

#ifdef _THERMIQUE_

   private :: &
      LoisThermoHydro_EnergieInterne_cv, & !  Enthalpie
      LoisThermoHydro_Enthalpie_cv, &
      LoisThermoHydro_SpecificEnthalpy_cv
#endif

contains

   ! main subroutine of this module
   ! compute all prim div for future use
   ! loop of cell/frac/node
   ! subroutine Loisthermohydro_divPrim_cv compute all prim div for one control volume
   subroutine LoisThermoHydro_compute() &
      bind(C, name="LoisThermoHydro_compute")

      ! cell
      call LoisThermoHydro_divPrim_cv(NbCellLocal_Ncpus(commRank + 1), IncCell, &
                                      CellRocktypeLocal, &
                                      !
                                      dXssurdXpCell, &
                                      SmdXsCell, &
                                      SmFCell, &
                                      !
                                      NumIncTotalPrimCell, &
                                      NumIncTotalSecondCell, &
                                      !
                                      DensiteMassiqueCell, &
                                      divDensiteMassiqueCell, &
                                      SmDensiteMassiqueCell, &
                                      !
                                      divPressionCell, &
                                      SmPressionCell, &
                                      !
                                      divTemperatureCell, &
                                      SmTemperatureCell, &
                                      !
                                      divSaturationCell, &
                                      !
                                      PressionCapCell, &
                                      divPressionCapCell, &
                                      !
                                      DensiteMolaireSatComp%cells, &
                                      divDensiteMolaireSatCompCell, &
                                      SmDensiteMolaireSatComp%cells, &
                                      !
                                      DensiteMolaireKrViscoCompCell, &
                                      divDensiteMolaireKrViscoCompCell, &
                                      SmDensiteMolaireKrViscoCompCell, &
                                      !
                                      DensiteMolaireEnergieInterneSat%cells, &
                                      divDensiteMolaireEnergieInterneSatCell, &
                                      SmDensiteMolaireEnergieInterneSatCell, &
                                      !
                                      DensiteMolaireKrViscoEnthalpieCell, &
                                      divDensiteMolaireKrViscoEnthalpieCell, &
                                      SmDensiteMolaireKrViscoEnthalpieCell)

      ! frac
      call LoisThermoHydro_divPrim_cv(NbFracLocal_Ncpus(commRank + 1), IncFrac, &
                                      FracRocktypeLocal, &
                                      !
                                      dXssurdXpFrac, &
                                      SmdXsFrac, &
                                      SmFFrac, &
                                      !
                                      NumIncTotalPrimFrac, &
                                      NumIncTotalSecondFrac, &
                                      !
                                      DensiteMassiqueFrac, &
                                      divDensiteMassiqueFrac, &
                                      SmDensiteMassiqueFrac, &
                                      !
                                      divPressionFrac, &
                                      SmPressionFrac, &
                                      !
                                      divTemperatureFrac, &
                                      SmTemperatureFrac, &
                                      !
                                      divSaturationFrac, &
                                      !
                                      PressionCapFrac, &
                                      divPressionCapFrac, &
                                      !
                                      DensiteMolaireSatComp%fractures, &
                                      divDensiteMolaireSatCompFrac, &
                                      SmDensiteMolaireSatComp%fractures, &
                                      !
                                      DensiteMolaireKrViscoCompFrac, &
                                      divDensiteMolaireKrViscoCompFrac, &
                                      SmDensiteMolaireKrViscoCompFrac, &
                                      !
                                      DensiteMolaireEnergieInterneSat%fractures, &
                                      divDensiteMolaireEnergieInterneSatFrac, &
                                      SmDensiteMolaireEnergieInterneSatFrac, &
                                      !
                                      DensiteMolaireKrViscoEnthalpieFrac, &
                                      divDensiteMolaireKrViscoEnthalpieFrac, &
                                      SmDensiteMolaireKrViscoEnthalpieFrac)

      ! node
      call LoisThermoHydro_divPrim_cv(NbNodeLocal_Ncpus(commRank + 1), IncNode, &
                                      NodeRocktypeLocal, &
                                      !
                                      dXssurdXpNode, &
                                      SmdXsNode, &
                                      SmFNode, &
                                      !
                                      NumIncTotalPrimNode, &
                                      NumIncTotalSecondNode, &
                                      !
                                      DensiteMassiqueNode, &
                                      divDensiteMassiqueNode, &
                                      SmDensiteMassiqueNode, &
                                      !
                                      divPressionNode, &
                                      SmPressionNode, &
                                      !
                                      divTemperatureNode, &
                                      SmTemperatureNode, &
                                      !
                                      divSaturationNode, &
                                      !
                                      PressionCapNode, &
                                      divPressionCapNode, &
                                      !
                                      DensiteMolaireSatComp%nodes, &
                                      divDensiteMolaireSatCompNode, &
                                      SmDensiteMolaireSatComp%nodes, &
                                      !
                                      DensiteMolaireKrViscoCompNode, &
                                      divDensiteMolaireKrViscoCompNode, &
                                      SmDensiteMolaireKrViscoCompNode, &
                                      !
                                      DensiteMolaireEnergieInterneSat%nodes, &
                                      divDensiteMolaireEnergieInterneSatNode, &
                                      SmDensiteMolaireEnergieInterneSatNode, &
                                      !
                                      DensiteMolaireKrViscoEnthalpieNode, &
                                      divDensiteMolaireKrViscoEnthalpieNode, &
                                      SmDensiteMolaireKrViscoEnthalpieNode)

#ifdef _WIP_FREEFLOW_STRUCTURES_
      ! FreeFlow nodes
      call LoisThermoHydro_divPrim_FreeFlow_cv(NbNodeLocal_Ncpus(commRank + 1), IncNode, &
                                               NodeRocktypeLocal, &
                                               !
                                               dXssurdXpNode, &
                                               SmdXsNode, &
                                               SmFNode, &
                                               !
                                               NumIncTotalPrimNode, &
                                               NumIncTotalSecondNode, &
                                               !
                                               divTemperatureNode, &
                                               SmTemperatureNode, &
                                               !
                                               divFreeFlowMolarFlowrateNode, &
                                               SmFreeFlowMolarFlowrateNode, &
                                               !
                                               FreeFlowMolarFlowrateCompNode, &
                                               divFreeFlowMolarFlowrateCompNode, &
                                               SmFreeFlowMolarFlowrateCompNode, &
                                               !
                                               FreeFlowHmCompNode, &
                                               divFreeFlowHmCompNode, &
                                               SmFreeFlowHmCompNode, &
                                               !
                                               FreeFlowMolarFlowrateEnthalpieNode, &
                                               divFreeFlowMolarFlowrateEnthalpieNode, &
                                               SmFreeFlowMolarFlowrateEnthalpieNode, &
                                               !
                                               FreeFlowHTTemperatureNetRadiationNode, &
                                               divFreeFlowHTTemperatureNetRadiationNode, &
                                               SmFreeFlowHTTemperatureNetRadiationNode, &
                                               !
                                               AtmEnthalpieNode)
#endif

      ! well injection
      !   compute q_{w,s,i} (and directive of pression)
      !   for each local well and its nodes
      call LoisThermoHydro_divP_wellinj(NbWellInjLocal_Ncpus(commRank + 1))

      ! SmDensiteMolaireKrViscoEnthalpieNode = 0.d0
      ! SmDensiteMolaireEnergieInterneSatNode = 0.d0
      ! SmDensiteMolaireKrViscoCompNode = 0.d0
      ! SmDensiteMolaireSatComp%nodes = 0.d0

      ! SmDensiteMolaireKrViscoEnthalpieFrac = 0.d0
      ! SmDensiteMolaireEnergieInterneSatFrac = 0.d0
      ! SmDensiteMolaireKrViscoCompFrac = 0.d0
      ! SmDensiteMolaireSatComp%fractures = 0.d0

      ! SmDensiteMolaireKrViscoEnthalpieCell = 0.d0
      ! SmDensiteMolaireEnergieInterneSatCell = 0.d0
      ! SmDensiteMolaireKrViscoCompCell = 0.d0
      ! SmDensiteMolaireSatComp%cells = 0.d0

      ! SmDensiteMassiqueCell = 0.d0
      ! SmDensiteMassiqueNode = 0.d0
      ! SmDensiteMassiqueFrac = 0.d0

   end subroutine LoisThermoHydro_compute

   ! all operations for one cv
   subroutine LoisThermoHydro_divPrim_cv( &
      NbIncLocal, inc, &
      rt, &
      dXssurdXp, SmdXs, SmF, &
      NumIncTotalPrimCV, NumIncTotalSecondCV, &
      DensiteMassique, divDensiteMassique, SmDensiteMassique, &
      divPression, SmPression, &
      divTemperature, SmTemperature, &
      divSaturation, &
      PressionCap, divPressionCap, &
      DensiteMolaireSatComp, divDensiteMolaireSatComp, SmDensiteMolaireSatComp, &
      DensiteMolaireKrViscoComp, divDensiteMolaireKrViscoComp, SmDensiteMolaireKrViscoComp, &
      DensiteMolaireEnergieInterneSat, divDensiteMolaireEnergieInterneSat, SmDensiteMolaireEnergieInterneSat, &
      DensiteMolaireKrViscoEnthalpie, divDensiteMolaireKrViscoEnthalpie, SmDensiteMolaireKrViscoEnthalpie)

      ! input
      integer, intent(in) :: NbIncLocal

      type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

      integer, intent(in) :: rt(IndThermique + 1, NbIncLocal)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax, NbIncLocal), &
         NumIncTotalSecondCV(NbEqFermetureMax, NbIncLocal)

      double precision, intent(in) :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs(NbEqFermetureMax, NbIncLocal), &
         SmF(NbEqFermetureMax, NbIncLocal)

      ! output

      double precision, intent(out) :: &
         !
         DensiteMassique(NbPhase, NbIncLocal), &
         divDensiteMassique(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMassique(NbPhase, NbIncLocal), &
         !
         divPression(NbIncTotalPrimMax, NbIncLocal), &
         SmPression(NbIncLocal), &
         !
         divTemperature(NbIncTotalPrimMax, NbIncLocal), &
         SmTemperature(NbIncLocal), &
         !
         divSaturation(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         !
         PressionCap(NbPhase, NbIncLocal), &
         divPressionCap(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         !
         DensiteMolaireSatComp(NbComp, NbPhase, NbIncLocal), &
         divDensiteMolaireSatComp(NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmDensiteMolaireSatComp(NbComp, NbPhase, NbIncLocal), &
         !
         DensiteMolaireKrViscoComp(NbComp, NbPhase, NbIncLocal), &
         divDensiteMolaireKrViscoComp(NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmDensiteMolaireKrViscoComp(NbComp, NbPhase, NbIncLocal), &
         !
         DensiteMolaireEnergieInterneSat(NbPhase, NbIncLocal), &
         divDensiteMolaireEnergieInterneSat(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMolaireEnergieInterneSat(NbPhase, NbIncLocal), &
         !
         DensiteMolaireKrViscoEnthalpie(NbPhase, NbIncLocal), &
         divDensiteMolaireKrViscoEnthalpie(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMolaireKrViscoEnthalpie(NbPhase, NbIncLocal)

      ! tmp
      double precision :: &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase), &
         PermRel(NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase), &
         !
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      double precision :: &
         EnergieInterne(NbPhase), &
         divEnergieInterne(NbIncTotalPrimMax, NbPhase), &
         SmEnergieInterne(NbPhase), &
         !
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncTotalPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)

      integer :: k, i, icp, iph
      type(ContextInfo) :: ctxinfo

      do k = 1, NbIncLocal
         call LoisThermoHydro_viscosite_cv( &
            inc(k), dXssurdXp(:, :, k), SmdXs(:, k), NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
            wa_UnsurViscosite(:, k), wa_divUnsurViscosite(:, :, k), wa_SmUnSurViscosite(:, k))
      end do

      do k = 1, NbIncLocal

         ! init tmp values for each cv
         call LoisThermoHydro_init_cv(inc(k), ctxinfo)

         ! densite massique
         call LoisThermoHydro_densitemassique_cv(inc(k), ctxinfo, dXssurdXp(:, :, k), SmdXs(:, k), &
                                                 NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                                 DensiteMassique(:, k), divDensiteMassique(:, :, k), SmDensiteMassique(:, k))

         ! deniste molaire
         call LoisThermoHydro_densitemolaire_cv(inc(k), ctxinfo, dXssurdXp(:, :, k), SmdXs(:, k), &
                                                NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                                DensiteMolaire, divDensiteMolaire, SmDensiteMolaire)

         ! PermRel
         call LoisThermoHydro_PermRel_cv(inc(k), ctxinfo, rt(1, k), dXssurdXp(:, :, k), SmdXs(:, k), &
                                         NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                         PermRel, divPermRel)

         ! Reference Pressure (unknown index is 1)
         call LoisThermoHydro_Inc_cv(1, inc(k), ctxinfo, &
                                     NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                     dXssurdXp(:, :, k), SmdXs(:, k), &
                                     divPression(:, k), SmPression(k))

         ! Comp
         do i = 2 + IndThermique, ctxinfo%NbIncPTC ! loop over index of Components
            iph = ctxinfo%NumIncPTC2NumIncComp_phase(i) ! phase corresponding to unknown i
            icp = ctxinfo%NumIncPTC2NumIncComp_comp(i) ! component corresponding to unknown i
            call LoisThermoHydro_Inc_cv(i, inc(k), ctxinfo, &
                                        NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                        dXssurdXp(:, :, k), SmdXs(:, k), &
                                        divComp(:, icp, iph), SmComp(icp, iph))
         enddo

         ! Pression Capillaire
         call LoisThermoHydro_PressionCapillaire_cv(rt(:, k), inc(k), ctxinfo, &
                                                    NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                                    dXssurdXp(:, :, k), &
                                                    PressionCap(:, k), divPressionCap(:, :, k))

         ! Saturation div FIXME: not done with LoisThermoHydro_Inc_cv because last saturation is eliminated
         call LoisThermoHydro_Saturation_cv(inc(k), ctxinfo, &
                                            NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                            dXssurdXp(:, :, k), &
                                            divSaturation(:, :, k))

         ! term: DensiteMolaire * PermRel / Viscosite * Comp
         call LoisThermoHydro_DensiteMolaireKrViscoComp_cv( &
            inc(k), ctxinfo, &
            DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
            PermRel, divPermRel, &
            wa_UnsurViscosite(:, k), wa_divUnsurViscosite(:, :, k), wa_SmUnSurViscosite(:, k), &
            divComp, SmComp, &
            NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
            dXssurdXp(:, :, k), SmdXs(:, k), &
            DensiteMolaireKrViscoComp(:, :, k), &
            divDensiteMolaireKrViscoComp(:, :, :, k), &
            SmDensiteMolaireKrViscoComp(:, :, k))

         ! term: DensiteMolaire * Saturation * Comp
         call LoisThermoHydro_DensiteMolaireSatComp_cv( &
            inc(k), ctxinfo, divSaturation(:, :, k), &
            DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
            divComp, SmComp, &
            NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
            dXssurdXp(:, :, k), SmdXs(:, k), &
            DensiteMolaireSatComp(:, :, k), &
            divDensiteMolaireSatComp(:, :, :, k), &
            SmDensiteMolaireSatComp(:, :, k))

#ifdef _THERMIQUE_
         ! FIXME: Temperature (unknown index is 2)
         call LoisThermoHydro_Inc_cv(2, inc(k), ctxinfo, &
                                     NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                     dXssurdXp(:, :, k), SmdXs(:, k), &
                                     divTemperature(:, k), SmTemperature(k))

         ! energie interne
         call LoisThermoHydro_EnergieInterne_cv(inc(k), ctxinfo, dXssurdXp(:, :, k), SmdXs(:, k), &
                                                NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                                EnergieInterne, divEnergieInterne, SmEnergieInterne)

         ! Enthalpie
         call LoisThermoHydro_Enthalpie_cv(inc(k), ctxinfo, dXssurdXp(:, :, k), SmdXs(:, k), &
                                           NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                           Enthalpie, divEnthalpie, SmEnthalpie)

         ! term: DensiteMolaire * Energieinterne * Saturation
         call LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv( &
            inc(k), ctxinfo, divSaturation(:, :, k), &
            DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
            EnergieInterne, divEnergieInterne, SmEnergieInterne, &
            DensiteMolaireEnergieInterneSat(:, k), &
            divDensiteMolaireEnergieInterneSat(:, :, k), &
            SmDensiteMolaireEnergieInterneSat(:, k))

         ! term: DensiteMolaire * PermRel / Viscosite * Enthalpie
         call LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv( &
            ctxinfo, &
            DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
            PermRel, divPermRel, &
            wa_UnsurViscosite(:, k), wa_divUnsurViscosite(:, :, k), wa_SmUnSurViscosite(:, k), &
            Enthalpie, divEnthalpie, SmEnthalpie, &
            DensiteMolaireKrViscoEnthalpie(:, k), &
            divDensiteMolaireKrViscoEnthalpie(:, :, k), &
            SmDensiteMolaireKrViscoEnthalpie(:, k))
#endif

      end do

   end subroutine LoisThermoHydro_divPrim_cv

#ifdef _WIP_FREEFLOW_STRUCTURES_
   ! Compute derivative of phase molar flowrates in the Freeflow dof
   subroutine LoisThermoHydro_divPrim_FreeFlow_cv( &
      NbIncLocal, inc, &
      rt, &
      dXssurdXp, SmdXs, SmF, &
      NumIncTotalPrimCV, NumIncTotalSecondCV, &
      divTemperature, SmTemperature, &
      divFreeFlowMolarFlowrate, SmFreeFlowMolarFlowrate, &
      FreeFlowMolarFlowrateComp, divFreeFlowMolarFlowrateComp, SmFreeFlowMolarFlowrateComp, &
      FreeFlowHmComp, divFreeFlowHmComp, SmFreeFlowHmComp, &
      FreeFlowMolarFlowrateEnthalpie, divFreeFlowMolarFlowrateEnthalpie, SmFreeFlowMolarFlowrateEnthalpie, &
      FreeFlowHTTemperatureNetRadiation, divFreeFlowHTTemperatureNetRadiation, SmFreeFlowHTTemperatureNetRadiation, &
      AtmEnthalpie)

      ! input
      integer, intent(in) :: NbIncLocal

      type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

      integer, intent(in) :: rt(IndThermique + 1, NbIncLocal)

      double precision, intent(in) :: &
         divTemperature(NbIncTotalPrimMax, NbIncLocal), &
         SmTemperature(NbIncLocal)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax, NbIncLocal), &
         NumIncTotalSecondCV(NbEqFermetureMax, NbIncLocal)

      double precision, intent(in) :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs(NbEqFermetureMax, NbIncLocal), &
         SmF(NbEqFermetureMax, NbIncLocal)

      ! output
      double precision, intent(out) :: &
         divFreeFlowMolarFlowrate(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmFreeFlowMolarFlowrate(NbPhase, NbIncLocal), &
         FreeFlowMolarFlowrateComp(NbComp, NbPhase, NbIncLocal), &
         divFreeFlowMolarFlowrateComp(NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmFreeFlowMolarFlowrateComp(NbComp, NbPhase, NbIncLocal), &
         FreeFlowHmComp(NbComp, NbPhase, NbIncLocal), &
         divFreeFlowHmComp(NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmFreeFlowHmComp(NbComp, NbPhase, NbIncLocal), &
         FreeFlowMolarFlowrateEnthalpie(NbPhase, NbIncLocal), &
         divFreeFlowMolarFlowrateEnthalpie(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmFreeFlowMolarFlowrateEnthalpie(NbPhase, NbIncLocal), &
         FreeFlowHTTemperatureNetRadiation(NbIncLocal), &
         divFreeFlowHTTemperatureNetRadiation(NbIncTotalPrimMax, NbIncLocal), &
         SmFreeFlowHTTemperatureNetRadiation(NbIncLocal), &
         AtmEnthalpie(NbPhase, NbIncLocal)

      ! tmp
      double precision :: &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase), &
         SpecificEnthalpy(NbComp, NbPhase), &
         divSpecificEnthalpy(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmSpecificEnthalpy(NbComp, NbPhase)
      integer :: k, i, iph, icp
      type(ContextInfo) :: ctxinfo

      ! loop over each local element, called only with nodes
      do k = 1, NbIncLocal
         ! subroutine only called with nodes
         if (.not. IdFFNodeLocal(k)) cycle ! loop over Freeflow dof only, avoid reservoir context

         ! init tmp values for each cv
         call LoisThermoHydro_init_cv(inc(k), ctxinfo)

         ! Comp
         do i = 2 + IndThermique, ctxinfo%NbIncPTC ! loop over index of Components
            iph = ctxinfo%NumIncPTC2NumIncComp_phase(i) ! phase of Component i
            icp = ctxinfo%NumIncPTC2NumIncComp_comp(i) ! numero of the component i
            call LoisThermoHydro_Inc_cv(i, inc(k), ctxinfo, &
                                        NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                        dXssurdXp(:, :, k), SmdXs(:, k), &
                                        divComp(:, icp, iph), SmComp(icp, iph))
         enddo

         do i = 1, ctxinfo%NbPhasePresente
            ! phase molar flowrate (unknown index is NbIncPTC+NbPhasePresente-1+i) (FIXME: -1 because one saturation is eliminated)
            call LoisThermoHydro_Inc_cv(ctxinfo%NbIncPTC + ctxinfo%NbPhasePresente - 1 + i, inc(k), ctxinfo, &
                                        NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                        dXssurdXp(:, :, k), SmdXs(:, k), &
                                        divFreeFlowMolarFlowrate(:, i, k), SmFreeFlowMolarFlowrate(i, k))
         enddo

         ! term: FreeFlowMolarFlowrate * Comp
         call LoisThermoHydro_FreeFlowMolarFlowrateComp_cv(inc(k), ctxinfo, &
                                                           divFreeFlowMolarFlowrate(:, :, k), SmFreeFlowMolarFlowrate(:, k), &
                                                           divComp, SmComp, &
                                                           FreeFlowMolarFlowrateComp(:, :, k), &
                                                           divFreeFlowMolarFlowrateComp(:, :, :, k), &
                                                           SmFreeFlowMolarFlowrateComp(:, :, k))

         ! term: Hm * (Comp - atm_comp)
         call LoisThermoHydro_FreeFlowHmComp_cv(inc(k), ctxinfo, &
                                                divComp, SmComp, &
                                                FreeFlowHmComp(:, :, k), &
                                                divFreeFlowHmComp(:, :, :, k), &
                                                SmFreeFlowHmComp(:, :, k))

#ifdef _THERMIQUE_
         ! SpecificEnthalpy
         call LoisThermoHydro_SpecificEnthalpy_cv(inc(k), ctxinfo, dXssurdXp(:, :, k), SmdXs(:, k), &
                                                  NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                                  SpecificEnthalpy, divSpecificEnthalpy, SmSpecificEnthalpy)

         ! term: gas-> FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)
         !       liquid-> FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)
         call LoisThermoHydro_FreeFlowMolarFlowrateEnthalpie_cv(inc(k), ctxinfo, &
                                                                divFreeFlowMolarFlowrate(:, :, k), SmFreeFlowMolarFlowrate(:, k), &
                                                                SpecificEnthalpy, divSpecificEnthalpy, SmSpecificEnthalpy, &
                                                                divComp, SmComp, &
                                                                FreeFlowMolarFlowrateEnthalpie(:, k), &
                                                                divFreeFlowMolarFlowrateEnthalpie(:, :, k), &
                                                                SmFreeFlowMolarFlowrateEnthalpie(:, k))

         ! term: HT * (T - atm_temperature) - net Radiation (which is a factor of T**4)
         call LoisThermoHydro_FreeFlowHTTemperatureNetRadiation_cv(inc(k), ctxinfo, &
                                                                   divTemperature(:, k), SmTemperature(k), &
                                                                   FreeFlowHTTemperatureNetRadiation(k), &
                                                                   divFreeFlowHTTemperatureNetRadiation(:, k), &
                                                                   SmFreeFlowHTTemperatureNetRadiation(k))

         ! term: gas-> SpecificEnthalpy(water, gas) of the far field atmosphere
         !       liquid-> Enthalpie(liquid) of the far field atmosphere
         call LoisThermoHydro_AtmEnthalpie_cv(ctxinfo, AtmEnthalpie(:, k))
#endif
      end do ! k

   end subroutine LoisThermoHydro_divPrim_FreeFlow_cv

   ! term: FreeFlowMolarFlowrate * Comp
   subroutine LoisThermoHydro_FreeFlowMolarFlowrateComp_cv( &
      inc, ctxinfo, &
      divFreeFlowMolarFlowrate, SmFreeFlowMolarFlowrate, &
      divComp, SmComp, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Comp and FreeFlow_flowrate
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         divFreeFlowMolarFlowrate(NbIncTotalPrimMax, NbPhase), &
         SmFreeFlowMolarFlowrate(NbPhase), &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      ! output
      double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncTotalPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

      ! tmp
      integer :: i, iph, icp, k

      val(:, :) = 0.d0
      dval(:, :, :) = 0.d0
      Smval(:, :) = 0.d0

      ! 1. val: FreeFlowMolarFlowrate * Comp
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp
            ! only {alpha | alpha \in Q_k \cap P_i} is useful
            ! To understand better, change the order of the loop do i=.. and the loop do icp=..
            if (MCP(icp, iph) == 1) then ! P_i
               val(icp, i) = inc%FreeFlow_flowrate(iph)*inc%Comp(icp, iph)
            end if
         end do
      end do

      ! 2. div(FreeFlowMolarFlowrate * Comp)
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then

               do k = 1, ctxinfo%NbIncTotalPrim
                  dval(k, icp, i) = divFreeFlowMolarFlowrate(k, i)*inc%Comp(icp, iph) &
                                    + inc%FreeFlow_flowrate(iph)*divComp(k, icp, iph)
               end do

            end if
         end do ! end of P_i

      end do ! end of 2

      ! 3. Sm(FreeFlowMolarFlowrate * Comp)
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then

               Smval(icp, i) = SmFreeFlowMolarFlowrate(i)*inc%Comp(icp, iph) &
                               + inc%FreeFlow_flowrate(iph)*SmComp(icp, iph)
            end if
         end do
      end do ! end of 3.

   end subroutine LoisThermoHydro_FreeFlowMolarFlowrateComp_cv

   ! term: Hm * (Comp - atm_comp)
   subroutine LoisThermoHydro_FreeFlowHmComp_cv( &
      inc, ctxinfo, &
      divComp, SmComp, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Comp
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      ! output
      double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncTotalPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

      ! tmp
      integer :: i, iph, icp, k

      val(:, :) = 0.d0
      dval(:, :, :) = 0.d0
      Smval(:, :) = 0.d0

      ! 1. val: Hm * (Comp - atm_comp) if gas ; Hm(alpha)=0. if alpha not gas
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp
            ! only {alpha | alpha \in Q_k \cap P_i} is useful
            ! To understand better, change the order of the loop do i=.. and the loop do icp=..
            if (MCP(icp, iph) == 1) then ! P_i
               val(icp, i) = Hm(iph)*(inc%Comp(icp, iph) - atm_comp(icp, iph))
            end if
         end do
      end do

      ! 2. div(Hm * (Comp - atm_comp))
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then

               do k = 1, ctxinfo%NbIncTotalPrim
                  dval(k, icp, i) = Hm(iph)*divComp(k, icp, iph)
               end do

            end if
         end do ! end of P_i

      end do ! end of 2

      ! 3. Sm(Hm * (Comp - atm_comp))
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then

               Smval(icp, i) = Hm(iph)*SmComp(icp, iph)
            end if
         end do
      end do ! end of 3.

   end subroutine LoisThermoHydro_FreeFlowHmComp_cv

#ifdef _THERMIQUE_
   ! term: gas-> FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)
   !       liquid-> FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)
   subroutine LoisThermoHydro_FreeFlowMolarFlowrateEnthalpie_cv( &
      inc, ctxinfo, &
      divFreeFlowMolarFlowrate, SmFreeFlowMolarFlowrate, &
      SpecificEnthalpy, divSpecificEnthalpy, SmSpecificEnthalpy, &
      divComp, SmComp, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Temperature and Comp
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         divFreeFlowMolarFlowrate(NbIncTotalPrimMax, NbPhase), &
         SmFreeFlowMolarFlowrate(NbPhase), &
         SpecificEnthalpy(NbComp, NbPhase), &
         divSpecificEnthalpy(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmSpecificEnthalpy(NbComp, NbPhase), &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      ! output
      double precision, intent(out) :: &
         val(NbPhase), &
         dval(NbIncTotalPrimMax, NbPhase), &
         Smval(NbPhase)

      ! tmp
      integer :: i, iph, icp, k

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         if (iph == LIQUID_PHASE) then ! FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)

            do icp = 1, NbComp
               ! 1. val: FreeFlowMolarFlowrate(liquid) * sum_icp( SpecificEnthalpy(icp,i)*Comp(icp,iph) )
               val(i) = val(i) &
                        + inc%FreeFlow_flowrate(iph)*SpecificEnthalpy(icp, i)*inc%Comp(icp, iph)

               ! 2. dval
               do k = 1, ctxinfo%NbIncTotalPrim
                  dval(k, i) = dval(k, i) &
                               + divFreeFlowMolarFlowrate(k, i)*SpecificEnthalpy(icp, i)*inc%Comp(icp, iph) &
                               + inc%FreeFlow_flowrate(iph)*divSpecificEnthalpy(k, icp, i)*inc%Comp(icp, iph) &
                               + inc%FreeFlow_flowrate(iph)*SpecificEnthalpy(icp, i)*divComp(k, icp, iph)
               enddo ! k

               ! 3. Smval
               Smval(i) = Smval(i) &
                          + SmFreeFlowMolarFlowrate(i)*SpecificEnthalpy(icp, i)*inc%Comp(icp, iph) &
                          + inc%FreeFlow_flowrate(iph)*SmSpecificEnthalpy(icp, i)*inc%Comp(icp, iph) &
                          + inc%FreeFlow_flowrate(iph)*SpecificEnthalpy(icp, i)*SmComp(icp, iph)
            enddo ! icp

         else ! gaz phase

            ! 1. val: FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)
            val(i) = inc%FreeFlow_flowrate(iph)*SpecificEnthalpy(WATER_COMP, i)

            ! 2. dval
            do k = 1, ctxinfo%NbIncTotalPrim
               dval(k, i) = divFreeFlowMolarFlowrate(k, i)*SpecificEnthalpy(WATER_COMP, i) &
                            + inc%FreeFlow_flowrate(iph)*divSpecificEnthalpy(k, WATER_COMP, i)
            enddo ! k
            ! 3. Smval
            Smval(i) = SmFreeFlowMolarFlowrate(i)*SpecificEnthalpy(WATER_COMP, i) &
                       + inc%FreeFlow_flowrate(iph)*SmSpecificEnthalpy(WATER_COMP, i)

         endif ! phase
      enddo

   end subroutine LoisThermoHydro_FreeFlowMolarFlowrateEnthalpie_cv

   ! term: HT * (T - atm_temperature) - net Radiation (which is a factor of T**4)
   subroutine LoisThermoHydro_FreeFlowHTTemperatureNetRadiation_cv( &
      inc, ctxinfo, &
      divTemperature, SmTemperature, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Temperature
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         divTemperature(NbIncTotalPrimMax), &
         SmTemperature

      ! output
      double precision, intent(out) :: &
         val, &
         dval(NbIncTotalPrimMax), &
         Smval

      integer :: k

      val = 0.d0
      dval(:) = 0.d0
      Smval = 0.d0

      ! 1. val: HT * (T - atm_temperature)
      !         - atm_flux_radiation + soil_emissivity*Stephan_Boltzmann_cst*T**4
      val = HT*(inc%Temperature - atm_temperature) &
            - atm_flux_radiation + soil_emissivity*Stephan_Boltzmann_cst*inc%Temperature**4.d0

      ! 2. dval
      do k = 1, ctxinfo%NbIncTotalPrim
         dval(k) = HT*divTemperature(k) &
                   + soil_emissivity*Stephan_Boltzmann_cst*4.d0*divTemperature(k)*inc%Temperature**3.d0
      enddo

      ! 3. Smval
      Smval = HT*SmTemperature &
              + soil_emissivity*Stephan_Boltzmann_cst*4.d0*SmTemperature*inc%Temperature**3.d0

   end subroutine LoisThermoHydro_FreeFlowHTTemperatureNetRadiation_cv

   ! term: gas-> SpecificEnthalpy(water, gas) of the far field atmosphere
   !       liquid-> Enthalpie(liquid) of the far field atmosphere
   subroutine LoisThermoHydro_AtmEnthalpie_cv(ctxinfo, val)

      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(out) :: val(NbPhase)

      ! tmp
      double precision :: f(NbComp), Sat(NbPhase), dPf(NbComp), dTf(NbComp), &
         dCf(NbComp, NbComp), dSf(NbComp, NbPhase)
      integer :: i, iph, icp

      Sat = 0.d0
      val = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         call f_SpecificEnthalpy(iph, atm_pressure, atm_temperature, &
                                 atm_comp(:, iph), Sat, & ! Sat not used
                                 f, dPf, dTf, dCf, dSf)

         if (iph == LIQUID_PHASE) then ! sum_icp( specific_enthalpie(icp)*atm_comp(icp,iph) )
            do icp = 1, NbComp
               val(i) = val(i) + f(icp)*atm_comp(icp, iph)
            enddo
         else  ! gaz phase : specific_enthalpie(water)
            val(i) = f(WATER_COMP)
         endif

      enddo ! NbPhasePresente

   end subroutine LoisThermoHydro_AtmEnthalpie_cv

! _THERMIQUE_
#endif
! _WIP_FREEFLOW_STRUCTURES_
#endif

   subroutine LoisThermoHydro_local_Schur(nb_unknows, nb_closures, nb_phases, dXssurdXp, dfdX_secd, dval, SmdXs, Smval)
      integer, intent(in) :: nb_unknows
      integer, intent(in) :: nb_closures
      integer, intent(in) :: nb_phases
      double precision, intent(in) :: dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)
      double precision, intent(in) :: dfdX_secd(NbEqFermetureMax, nb_phases)
      double precision, intent(inout) :: dval(NbIncTotalPrimMax, nb_phases)
      double precision, optional, intent(in) :: SmdXs(NbEqFermetureMax)
      double precision, optional, intent(inout) :: Smval(nb_phases)

      ! dv/dXp - dv/dXs*dXs/dXp
      ! dval = dfdX_prim - dXssurdXp*dfdX_secd
      ! all the mats is in (col, row) index order, only need to consider as transpose
      call dgemm('N', 'N', &
                 nb_unknows, nb_phases, nb_closures, &
                 -1.d0, dXssurdXp, NbIncTotalPrimMax, dfdX_secd, &
                 NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

#ifndef NDEBUG
      if ((present(SmdXs) .and. .not. present(Smval)) .or. (present(Smval) .and. .not. present(SmdXs))) &
         call CommonMPI_abort("Both second members arguments should be provided!")
#endif

      ! - dv/dXs*SmdXs
      if (present(SmdXs)) then
         call dgemv('T', nb_closures, nb_phases, &
                    -1.d0, dfdX_secd, NbEqFermetureMax, &
                    SmdXs, 1, 0.d0, Smval, 1)
      end if

   end subroutine LoisThermoHydro_local_Schur

   !> Update thermo Laws of nodes
   subroutine LoisThermoHydro_divPrim_nodes

      call LoisThermoHydro_divPrim_cv(NbNodeLocal_Ncpus(commRank + 1), IncNode, &
                                      NodeRocktypeLocal, &
                                      !
                                      dXssurdXpNode, &
                                      SmdXsNode, &
                                      SmFNode, &
                                      !
                                      NumIncTotalPrimNode, &
                                      NumIncTotalSecondNode, &
                                      !
                                      DensiteMassiqueNode, &
                                      divDensiteMassiqueNode, &
                                      SmDensiteMassiqueNode, &
                                      !
                                      divPressionNode, &
                                      SmPressionNode, &
                                      !
                                      divTemperatureNode, &
                                      SmTemperatureNode, &
                                      !
                                      divSaturationNode, &
                                      !
                                      PressionCapNode, &
                                      divPressionCapNode, &
                                      !
                                      DensiteMolaireSatComp%nodes, &
                                      divDensiteMolaireSatCompNode, &
                                      SmDensiteMolaireSatComp%nodes, &
                                      !
                                      DensiteMolaireKrViscoCompNode, &
                                      divDensiteMolaireKrViscoCompNode, &
                                      SmDensiteMolaireKrViscoCompNode, &
                                      !
                                      DensiteMolaireEnergieInterneSat%nodes, &
                                      divDensiteMolaireEnergieInterneSatNode, &
                                      SmDensiteMolaireEnergieInterneSatNode, &
                                      !
                                      DensiteMolaireKrViscoEnthalpieNode, &
                                      divDensiteMolaireKrViscoEnthalpieNode, &
                                      SmDensiteMolaireKrViscoEnthalpieNode)

   end subroutine LoisThermoHydro_divPrim_nodes

   subroutine LoisThermoHydro_init_cv(inc, ctxinfo)

      type(TYPE_IncCVReservoir), intent(in) :: inc
      type(ContextInfo), intent(out) :: ctxinfo

      ctxinfo%NbPhasePresente = NbPhasePresente_ctx(inc%ic)
      ctxinfo%NbCompCtilde = NbCompCtilde_ctx(inc%ic)

      ctxinfo%NbEqFermeture = NbEqFermeture_ctx(inc%ic)
      ctxinfo%NbEqEquilibre = NbEqEquilibre_ctx(inc%ic)

      ctxinfo%NbIncPTC = NbIncPTC_ctx(inc%ic)
      ctxinfo%NbIncPTCPrim = ctxinfo%NbIncPTC - ctxinfo%NbEqFermeture

      ! ps. if there is only one phase, phase is secd
      ctxinfo%NbIncTotalPrim = NbIncTotalPrim_ctx(inc%ic)

      ctxinfo%NumPhasePresente(:) = NumPhasePresente_ctx(:, inc%ic)
      ctxinfo%NumCompCtilde(:) = NumCompCtilde_ctx(:, inc%ic)
      ctxinfo%NumCompEqEquilibre(:) = NumCompEqEquilibre_ctx(:, inc%ic)
      ctxinfo%NumIncPTC2NumIncComp_comp(:) = NumIncPTC2NumIncComp_comp_ctx(:, inc%ic)
      ctxinfo%NumIncPTC2NumIncComp_phase(:) = NumIncPTC2NumIncComp_phase_ctx(:, inc%ic)

      ctxinfo%Num2PhasesEqEquilibre(:, :) = Num2PhasesEqEquilibre_ctx(:, :, inc%ic)
      ctxinfo%NumIncComp2NumIncPTC(:, :) = NumIncComp2NumIncPTC_ctx(:, :, inc%ic)

   end subroutine LoisThermoHydro_init_cv

   ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
   subroutine LoisThermoHydro_fill_gradient_dfdX(context, iph, dPf, dTf, dCf, dSf, dfdX)
      integer(c_int), intent(in) :: context
      integer, intent(in) :: iph ! num of phase
      double precision, intent(in) :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision, intent(out) :: dfdX(NbIncTotalMax)

      integer(c_int) :: j, jc, jph, NbIncPTC, nb_phases
      double precision :: dfS_elim
      integer(c_int), pointer :: NumIncComp2NumIncPTC(:, :)

      NumIncComp2NumIncPTC => NumIncComp2NumIncPTC_ctx(:, :, context)

      dfdX = 0.d0

      dfdX(1) = dPf  ! P
#ifdef _THERMIQUE_
      dfdX(2) = dTf
#endif

      do j = 1, NbComp
         if (MCP(j, iph) == 1) then
            jc = NumIncComp2NumIncPTC(j, iph)
            dfdX(jc) = dCf(j)
         end if
      enddo

      nb_phases = NbPhasePresente_ctx(context)

#ifndef NDEBUG
      if (nb_phases < 1) &
         call CommonMPI_abort("LoisThermoHydro_fill_gradient_dfdX: Inconsistent number of phases.")
#endif

      if (nb_phases > 1) then
         dfS_elim = dSf(NumPhasePresente_ctx(nb_phases, context)) ! last saturation is eliminated
         NbIncPTC = NbIncPTC_ctx(context)
         do j = 1, nb_phases - 1
            jph = NumPhasePresente_ctx(j, context)
            jc = j + NbIncPTC
            dfdX(jc) = dSf(jph) - dfS_elim ! last saturation is eliminated
         enddo
         ! else nothing done because S=1, dS=0, and last saturation is eliminated
      endif

   end subroutine LoisThermoHydro_fill_gradient_dfdX

   subroutine LoisThermoHydro_dfdX_ps(context, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, dfdX_prim, dfdX_secd)
      integer(c_int), intent(in) :: context
      integer, intent(in) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(in) :: NumIncTotalSecondCV(NbEqFermetureMax)
      double precision, intent(in) :: dfdX(NbIncTotalMax)  ! dfdX = (df/dP, df/dT, df/dC, df/dS)
      double precision, intent(out) :: dfdX_prim(NbIncTotalPrimMax)
      double precision, intent(out) ::  dfdX_secd(NbEqFermetureMax)

      integer :: j, np, ns

      ! prim and secd part of dfdX
      np = NbIncTotalPrim_ctx(context)
      do j = 1, np
         dfdX_prim(j) = dfdX(NumIncTotalPrimCV(j))
      end do

      ! secd unknowns
      ns = NbEqFermeture_ctx(context)
      do j = 1, ns
         dfdX_secd(j) = dfdX(NumIncTotalSecondCV(j))
      end do

   end subroutine LoisThermoHydro_dfdX_ps

   subroutine LoisThermoHydro_densitemassique_cv(inc, ctxinfo, dXssurdXp, SmdXs, &
                                                 NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &    ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)
      double precision, intent(out) :: Smval(NbPhase)

      ! tmp
      double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase) !=NbEqFermetureMax

      integer :: iph, i

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

#ifdef DEBUG_LOISTHEMOHYDRO
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         write (*, *) 'Phase', iph, 'MCP row=', MCP(:, iph)
      end do
#endif

      dfdX_secd = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         call f_DensiteMassique(iph, inc%Pression, inc%Temperature, &
                                inc%Comp(:, iph), inc%Saturation, &
                                f, dPf, dTf, dCf, dSf)

         val(iph) = f ! val

         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf, dTf, dCf, dSf, dfdX)

         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
                                      dval(:, iph), dfdX_secd(:, iph))
      end do

      call LoisThermoHydro_local_Schur( &
         ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
         dXssurdXp, dfdX_secd, dval, SmdXs, Smval)

   end subroutine LoisThermoHydro_densitemassique_cv

   subroutine LoisThermoHydro_viscosite_cv( &
      inc, dXssurdXp, SmdXs, NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      double precision, intent(in) :: dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)
      double precision, intent(in) :: SmdXs(NbEqFermetureMax)
      integer, intent(in) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(in) :: NumIncTotalSecondCV(NbEqFermetureMax)
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)
      double precision, intent(out) :: Smval(NbPhase)

      integer :: i, iph, context, nb_phases
      double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase)

      val = 0.d0
      dval = 0.d0
      Smval = 0.d0
      dfdX_secd = 0.d0
      context = inc%ic
      nb_phases = NbPhasePresente_ctx(context)
      do i = 1, nb_phases
         iph = NumPhasePresente_ctx(i, context)
         call f_Viscosite(iph, inc%Pression, inc%Temperature, &
                          inc%Comp(:, iph), inc%Saturation, &
                          f, dPf, dTf, dCf, dSf)
         val(i) = 1.d0/f
         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS) - iph because we use MCP
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf, dTf, dCf, dSf, dfdX)
         dfdX(:) = -dfdX(:)/f**2
         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps( &
            inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, dval(:, i), dfdX_secd(:, i))
      end do

      call LoisThermoHydro_local_Schur( &
         NbIncTotalPrim_ctx(context), NbEqFermeture_ctx(context), NbPhase, &
         dXssurdXp, dfdX_secd, dval, SmdXs, Smval)

   end subroutine LoisThermoHydro_viscosite_cv

   subroutine LoisThermoHydro_densitemolaire_cv(inc, ctxinfo, dXssurdXp, SmdXs, &
                                                NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)
      double precision, intent(out) :: Smval(NbPhase)

      ! tmp
      double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase)

      integer :: i, iph

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      dfdX_secd(:, :) = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         call f_DensiteMolaire(iph, inc%Pression, inc%Temperature, &
                               inc%Comp(:, iph), inc%Saturation, &
                               f, dPf, dTf, dCf, dSf)

         val(i) = f ! val

         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf, dTf, dCf, dSf, dfdX)

         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
                                      dval(:, i), dfdX_secd(:, i))
      end do

      call LoisThermoHydro_local_Schur( &
         ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
         dXssurdXp, dfdX_secd, dval, SmdXs, Smval)

   end subroutine LoisThermoHydro_densitemolaire_cv

   subroutine LoisThermoHydro_all_relative_permeabilities(rock_type, phases, S, f, dSf)
      integer, intent(in) :: rock_type
      integer, intent(in) :: phases(:)
      double precision, intent(in) :: S(NbPhase)
      double precision, intent(out) :: f(NbPhase)
      double precision, intent(out) :: dSf(NbPhase, NbPhase)

      integer i

      do i = 1, size(phases)
         call f_PermRel(rock_type, phases(i), S, f(i), dSf(:, i))
      end do

   end subroutine LoisThermoHydro_all_relative_permeabilities

   subroutine LoisThermoHydro_PermRel_cv(inc, ctxinfo, rock_type, dXssurdXp, SmdXs, &
                                         NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval)
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      integer, intent(in) :: rock_type
      double precision, intent(in) :: dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)
      double precision, intent(in) :: SmdXs(NbEqFermetureMax)
      integer, intent(in) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(in) :: NumIncTotalSecondCV(NbEqFermetureMax)
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)

      double precision :: f(NbPhase), dSf(NbPhase, NbPhase), dCf(NbComp)
      double precision :: dfdX(NbIncTotalMax, NbPhase)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase)
      integer :: i, iph

      val(:) = 0.d0
      dval(:, :) = 0.d0
      dfdX_secd(:, :) = 0.d0

      ! WARNING: no variations with respect to molar fractions
      dCf = 0.d0

      call LoisThermoHydro_all_relative_permeabilities( &
         rock_type, ctxinfo%NumPhasePresente(1:ctxinfo%NbPhasePresente), &
         inc%Saturation, val, dSf)

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, 0.d0, 0.d0, dCf, dSf(:, i), dfdX(:, i))
         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX(:, i), &
                                      dval(:, i), dfdX_secd(:, i))
      end do

      ! FIXME: why not calling on RHS?
      !        because we don't use directly kr variations in Jacobian but divKrVisco...
      call LoisThermoHydro_local_Schur( &
         ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
         dXssurdXp, dfdX_secd, dval)

   end subroutine LoisThermoHydro_PermRel_cv

   subroutine LoisThermoHydro_Inc_cv(index_inc, inc, ctxinfo, &
                                     NumIncTotalPrimCV, NumIncTotalSecondCV, &
                                     dXssurdXp, SmdXs, &
                                     dval, Smval)

      ! input
      integer, intent(in) :: index_inc
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: dval(NbIncTotalPrimMax)
      double precision, intent(out) :: Smval

      integer :: i

      dval(:) = 0.d0
      Smval = 0.d0

      IF (ANY(NumIncTotalPrimCV == index_inc)) THEN ! index_inc (P or T) is prim
         do i = 1, NbIncTotalPrimMax
            if (NumIncTotalPrimCV(i) == index_inc) then
               dval(i) = 1.d0
            endif
         enddo
      ELSE IF (ANY(NumIncTotalSecondCV == index_inc)) THEN ! index_inc (P or T) is secd
         do i = 1, ctxinfo%NbEqFermeture
            if (NumIncTotalSecondCV(i) == index_inc) then
               dval(:) = -dXssurdXp(:, i)
               Smval = -SmdXs(i)
            endif
         enddo
      ELSE ! index_inc (P or T) not found
         write (*, *) "Looking for unknown ", index_inc, &
            "(NumIncTotalPrimCV=", NumIncTotalPrimCV, &
            " NumIncTotalSecondCV=", NumIncTotalSecondCV, ")"
         if (index_inc == 1) call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, P not found ')
#ifdef _THERMIQUE_
         if (index_inc == 1 + IndThermique) call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, T not found ')
#endif
         if (index_inc > 1 + IndThermique .and. index_inc <= ctxinfo%NbIncPTC) &
            call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, C not found ')
      ENDIF

   end subroutine LoisThermoHydro_Inc_cv

   subroutine LoisThermoHydro_Saturation_cv(inc, ctxinfo, &
                                            NumIncTotalPrimCV, NumIncTotalSecondCV, &
                                            dXssurdXp, &
                                            dval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)

      ! output
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)

      ! tmp
      integer :: i, k

      dval(:, :) = 0.d0

      ! FIXME: Elimination of the last present phase (sum S =1 forced in the code)
      DO i = 1, ctxinfo%NbPhasePresente - 1

         IF (ANY(NumIncTotalPrimCV == i + ctxinfo%NbIncPTC)) THEN ! S is prim
            do k = 1, NbIncTotalPrimMax
               if (NumIncTotalPrimCV(k) == i + ctxinfo%NbIncPTC) then
                  dval(k, i) = 1.d0

                  dval(k, ctxinfo%NbPhasePresente) = -1.d0
               endif
            enddo
         ELSE if (inc%ic < 2**NbPhase) then ! FIXME: avoid freeflow nodes: S not found in reservoir dof
            call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, S not found ')
         ENDIF
      ENDDO

   end subroutine LoisThermoHydro_Saturation_cv

   ! In the sens P(iph) = Pref + f_PressionCapillaire(iph)
   subroutine LoisThermoHydro_PressionCapillaire_cv(rt, inc, ctxinfo, &
                                                    NumIncTotalPrimCV, NumIncTotalSecondCV, &
                                                    dXssurdXp, &
                                                    val, dval)

      ! input
      integer, intent(in) :: rt(IndThermique + 1)
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)

      ! tmp
      double precision :: f, dSf(NbPhase), dfS_secd
      integer :: iph, j, jph, k

      val(:) = 0.d0
      dval(:, :) = 0.d0

      do iph = 1, NbPhase

         call f_PressionCapillaire(rt, iph, inc%Saturation, f, dSf)

         val(iph) = f

         dfS_secd = dSf(ctxinfo%NumPhasePresente(ctxinfo%NbPhasePresente))

         do j = 1, ctxinfo%NbPhasePresente - 1
            ! Look for S, is it primary or secondary unknowns ?
            ! FIXME: Elimination of the last present phase (sum S =1 forced in the code)
            IF (ANY(NumIncTotalPrimCV == j + ctxinfo%NbIncPTC)) THEN ! S is prim
               jph = ctxinfo%NumPhasePresente(j)
               do k = 1, NbIncTotalPrimMax
                  if (NumIncTotalPrimCV(k) == j + ctxinfo%NbIncPTC) then
                     dval(k, iph) = dSf(jph) - dfS_secd
                  endif
               enddo
            ELSE if (inc%ic < 2**NbPhase) then ! FIXME: avoid freeflow nodes: S not found in reservoir dof
               call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, S not found ')
            ENDIF
         enddo

      end do ! iph

   end subroutine LoisThermoHydro_PressionCapillaire_cv

#ifdef _THERMIQUE_

   subroutine LoisThermoHydro_EnergieInterne_cv(inc, ctxinfo, dXssurdXp, SmdXs, &
                                                NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)
      double precision, intent(out) :: Smval(NbPhase)

      ! tmp
      double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase)

      integer :: i, iph

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      dfdX_secd(:, :) = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         call f_EnergieInterne(iph, inc%Pression, inc%Temperature, &
                               inc%Comp(:, iph), inc%Saturation, &
                               f, dPf, dTf, dCf, dSf)

         val(i) = f ! val

         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf, dTf, dCf, dSf, dfdX)

         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
                                      dval(:, i), dfdX_secd(:, i))
      end do

      call LoisThermoHydro_local_Schur( &
         ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
         dXssurdXp, dfdX_secd, dval, SmdXs, Smval)

   end subroutine LoisThermoHydro_EnergieInterne_cv

   subroutine LoisThermoHydro_Enthalpie_cv(inc, ctxinfo, dXssurdXp, SmdXs, &
                                           NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)
      double precision, intent(out) :: Smval(NbPhase)

      ! tmp
      double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbPhase)

      integer :: i, iph

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      dfdX_secd(:, :) = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         call f_Enthalpie(iph, inc%Pression, inc%Temperature, &
                          inc%Comp(:, iph), inc%Saturation, &
                          f, dPf, dTf, dCf, dSf)

         val(i) = f ! val

         ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
         call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf, dTf, dCf, dSf, dfdX)

         ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
         ! and dfdX_secd w.r.t. the secondary unknowns
         call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
                                      dval(:, i), dfdX_secd(:, i))
      end do

      call LoisThermoHydro_local_Schur( &
         ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
         dXssurdXp, dfdX_secd, dval, SmdXs, Smval)

   end subroutine LoisThermoHydro_Enthalpie_cv

   ! Specific enthalpy
   subroutine LoisThermoHydro_SpecificEnthalpy_cv(inc, ctxinfo, dXssurdXp, SmdXs, &
                                                  NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: val(NbComp, NbPhase)
      double precision, intent(out) :: dval(NbIncTotalPrimMax, NbComp, NbPhase)
      double precision, intent(out) :: Smval(NbComp, NbPhase)

      ! tmp
      double precision :: f(NbComp), dPf(NbComp), dTf(NbComp), &
         dCf(NbComp, NbComp), dSf(NbComp, NbPhase)
      double precision :: dfdX(NbIncTotalMax)
      double precision :: dfdX_secd(NbEqFermetureMax, NbComp, NbPhase)

      integer :: i, iph, icp

      val = 0.d0
      dval = 0.d0
      Smval = 0.d0

      dfdX_secd = 0.d0

      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         call f_SpecificEnthalpy(iph, inc%Pression, inc%Temperature, &
                                 inc%Comp(:, iph), inc%Saturation, &
                                 f, dPf, dTf, dCf, dSf)

         do icp = 1, NbComp
            val(icp, i) = f(icp) ! val
            ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
            call LoisThermoHydro_fill_gradient_dfdX(inc%ic, iph, dPf(icp), dTf(icp), dCf(icp, :), dSf(icp, :), dfdX)
            ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
            ! and dfdX_secd w.r.t. the secondary unknowns
            call LoisThermoHydro_dfdX_ps(inc%ic, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
                                         dval(:, icp, i), dfdX_secd(:, icp, i))
         end do
      end do

      do icp = 1, NbComp
         call LoisThermoHydro_local_Schur( &
            ctxinfo%NbIncTotalPrim, ctxinfo%NbEqFermeture, NbPhase, &
            dXssurdXp, dfdX_secd(:, icp, :), dval(:, icp, :), SmdXs, Smval(icp, :))
      end do

   end subroutine LoisThermoHydro_SpecificEnthalpy_cv

#endif

   ! LoisThermoHydro_DensiteMolaireSat_cv is never called (cf. issue #292)
   ! term: desitemolaire * Saturation
!    subroutine LoisThermoHydro_DensiteMolaireSat_cv( &
!       inc, ctxinfo, divSaturation, &
!       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
!       val, dval, Smval)

!       ! input
!       type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation
!       type(ContextInfo), intent(in) :: ctxinfo
!       double precision, intent(in) :: &
!          divSaturation(NbIncTotalPrimMax, NbPhase), &
!          DensiteMolaire(NbPhase), &
!          divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
!          SmDensiteMolaire(NbPhase)

!       ! output
!       double precision, intent(out) :: &
!          val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

!       ! tmp
!       integer :: i, iph, k

!       val(:) = 0.d0
!       dval(:, :) = 0.d0
!       Smval(:) = 0.d0

!       ! 1. val
!       do i = 1, ctxinfo%NbPhasePresente
!          iph = ctxinfo%NumPhasePresente(i)
!          val(i) = DensiteMolaire(i)*inc%Saturation(iph)
!       end do

!       ! 2. dval: d(DensiteMolaire * Saturation)
!       do i = 1, ctxinfo%NbPhasePresente
!          iph = ctxinfo%NumPhasePresente(i)
!          do k = 1, ctxinfo%NbIncTotalPrim
!             dval(k, i) = &
!                divDensiteMolaire(k, i)*inc%Saturation(iph) &
!                + divSaturation(k, i)*DensiteMolaire(i)
!          end do
!       end do

!       ! 3. Smval
!       do i = 1, ctxinfo%NbPhasePresente
!          iph = ctxinfo%NumPhasePresente(i)
!          Smval(i) = SmDensiteMolaire(i)*inc%Saturation(iph)
!       end do

!    end subroutine LoisThermoHydro_DensiteMolaireSat_cv

   ! term: desitemolaire * Saturation * Comp
   subroutine LoisThermoHydro_DensiteMolaireSatComp_cv( &
      inc, ctxinfo, divSaturation, &
      DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
      divComp, SmComp, &
      NumIncTotalPrimCV, NumIncTotalSecondCV, &
      dXssurdXp, SmdXs, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         divSaturation(NbIncTotalPrimMax, NbPhase), &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase), &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncTotalPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

      ! tmp
      integer :: i, iph, icp, k, j, jcp, jph, numj, s
      double precision :: dv, tmp_val

      double precision :: dvi(NbIncTotalPrimMax)

      val(:, :) = 0.d0
      dval(:, :, :) = 0.d0
      Smval(:, :) = 0.d0

      ! 1. val
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         do icp = 1, NbComp
            ! only {alpha | alpha \in Q_k \cap P_i} is useful
            ! To understant better, change the order of the loop do i=.. and the loop do icp=..
            if (MCP(icp, iph) == 1) then ! P_i
               val(icp, i) = DensiteMolaire(i)*inc%Saturation(iph) &
                             *inc%Comp(icp, iph)
            end if
         end do
      end do

      ! 2.1 div(DensiteMolaire * Saturation * C_i^alpha)
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         ! 2.1.1 compute dvi, tmp vector, used in 2.1.2
         do k = 1, ctxinfo%NbIncTotalPrim
            dvi(k) = &
               divDensiteMolaire(k, i)*inc%Saturation(iph) &
               + divSaturation(k, i)*DensiteMolaire(i)
         end do
         tmp_val = DensiteMolaire(i)*inc%Saturation(iph)

         ! 2.1.2
         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then
               do k = 1, ctxinfo%NbIncTotalPrim
                  dval(k, icp, i) = dvi(k)*inc%Comp(icp, iph) &
                                    + tmp_val*divComp(k, icp, iph)
               end do
            end if
         end do ! end of P_i

      end do ! end of 2.1

      ! 3. Smval
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         dv = SmDensiteMolaire(i)*inc%Saturation(iph)
         tmp_val = DensiteMolaire(i)*inc%Saturation(iph)

         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then

               Smval(icp, i) = dv*inc%Comp(icp, iph) &
                               + tmp_val*SmComp(icp, iph)
            end if
         end do
      end do ! end of 2.2

   end subroutine LoisThermoHydro_DensiteMolaireSatComp_cv

   ! FIXME: LoisThermoHydro_DensiteMolaireKrVisco_cv is never called (cf. issue #291)
   ! div and Sm of term: DensiteMolaire*Kr/Viscosite
!    subroutine LoisThermoHydro_DensiteMolaireKrVisco_cv( &
!       ctxinfo, &
!       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
!       PermRel, divPermRel, &
!       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
!       val, dval, Smval)

!       ! input
!       type(ContextInfo), intent(in) :: ctxinfo
!       double precision, intent(in) :: &
!          DensiteMolaire(NbPhase), &
!          PermRel(NbPhase), &
!          UnsurViscosite(NbPhase), &
!          !
!          divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
!          divPermRel(NbIncTotalPrimMax, NbPhase), &
!          divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
!          !
!          SmDensiteMolaire(NbPhase), &
!          SmUnSurViscosite(NbPhase)

!       ! output
!       double precision, intent(out) :: &
!          val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

!       ! tmp
!       integer :: i, k

!       val(:) = 0.d0
!       dval(:, :) = 0.d0
!       Smval(:) = 0.d0

!       ! 1. val
!       do i = 1, ctxinfo%NbPhasePresente
!          val(i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
!       end do

!       ! 2. dval: d(DensiteMolaire*Kr/Viscosite)
!       do i = 1, ctxinfo%NbPhasePresente
!          do k = 1, ctxinfo%NbIncTotalPrim
!             dval(k, i) = &
!                divDensiteMolaire(k, i)*PermRel(i)*UnsurViscosite(i) &
!                + divPermRel(k, i)*DensiteMolaire(i)*UnsurViscosite(i) &
!                + divUnsurViscosite(k, i)*DensiteMolaire(i)*PermRel(i)
!          end do
!       end do

!       ! 3. Smval
!       do i = 1, ctxinfo%NbPhasePresente
!          Smval(i) = &
!             SmDensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
!             + SmUnSurViscosite(i)*DensiteMolaire(i)*PermRel(i)
!       end do

!    end subroutine LoisThermoHydro_DensiteMolaireKrVisco_cv

   ! div and Sm of term: DensiteMolaire*Kr/Viscosite*Comp
   subroutine LoisThermoHydro_DensiteMolaireKrViscoComp_cv(inc, ctxinfo, &
                                                           DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
                                                           PermRel, divPermRel, &
                                                           UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
                                                           divComp, SmComp, &
                                                           NumIncTotalPrimCV, NumIncTotalSecondCV, &
                                                           dXssurdXp, SmdXs, &
                                                           val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase), &
         !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase), &
         divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         !
         SmDensiteMolaire(NbPhase), &
         SmUnSurViscosite(NbPhase), &
         SmComp(NbComp, NbPhase)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncTotalPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

      ! tmp
      integer :: i, iph, icp, k, j, jcp, jph, numj, s
      double precision :: dv, tmp_val

      double precision :: dvi(NbIncTotalPrimMax)

      val(:, :) = 0.d0
      dval(:, :, :) = 0.d0
      Smval(:, :) = 0.d0

      ! 1. val
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         do icp = 1, NbComp
            ! only {alpha | alpha \in Q_k \cap P_i} is useful
            ! To understand better, change the order of the loop do i=.. and the loop do icp=..
            if (MCP(icp, iph) == 1) then ! P_i
               val(icp, i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
                             *inc%Comp(icp, iph)
            end if
         end do
      end do

      ! 2.1 div(DensiteMolaire*Kr/Viscosite*C_i^alpha)
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         ! 2.1.1 compute dvi=div(DensiteMolaire*Kr/Viscosite), tmp vector
         do k = 1, ctxinfo%NbIncTotalPrim
            dvi(k) = &
               divDensiteMolaire(k, i)*PermRel(i)*UnsurViscosite(i) &
               + divPermRel(k, i)*DensiteMolaire(i)*UnsurViscosite(i) &
               + divUnsurViscosite(k, i)*DensiteMolaire(i)*PermRel(i)
         end do

         tmp_val = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)

         ! 2.1.2
         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then
               do k = 1, ctxinfo%NbIncTotalPrim
                  dval(k, icp, i) = dvi(k)*inc%Comp(icp, iph) &
                                    + tmp_val*divComp(k, icp, iph)
               end do
            end if
         end do ! end of P_i

      end do ! end of 2.1

      ! 2.2. Sm(DensiteMolaire*Kr/Viscosite*C_i^alpha)
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         ! 2.2.1 compute dv=Sm(DensiteMolaire*Kr/Viscosite), tmp vector
         dv = &
            SmDensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
            + SmUnSurViscosite(i)*DensiteMolaire(i)*PermRel(i)
         tmp_val = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
         ! 2.2.2
         do icp = 1, NbComp ! P_i
            if (MCP(icp, iph) == 1) then
               Smval(icp, i) = dv*inc%Comp(icp, iph) &
                               + tmp_val*SmComp(icp, iph)
            end if
         end do
      end do ! end of 2.2

   end subroutine LoisThermoHydro_DensiteMolaireKrViscoComp_cv

   subroutine LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv( &
      inc, ctxinfo, divSaturation, &
      DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
      EnergieInterne, divEnergieInterne, SmEnergieInterne, &
      val, dval, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation
      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         EnergieInterne(NbPhase), &
         !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divEnergieInterne(NbIncTotalPrimMax, NbPhase), &
         divSaturation(NbIncTotalPrimMax, NbPhase), &
         !
         SmDensiteMolaire(NbPhase), &
         SmEnergieInterne(NbPhase)

      ! output
      double precision, intent(out) :: &
         val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

      ! tmp
      integer :: i, iph, k

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      ! 1. val
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         val(i) = DensiteMolaire(i)*EnergieInterne(i)*inc%Saturation(iph)
      end do

      ! 2. dval
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)

         do k = 1, ctxinfo%NbIncTotalPrim

            dval(k, i) = &
               divDensiteMolaire(k, i)*EnergieInterne(i)*inc%Saturation(iph) &
               + divEnergieInterne(k, i)*DensiteMolaire(i)*inc%Saturation(iph) &
               + divSaturation(k, i)*DensiteMolaire(i)*EnergieInterne(i)
         end do
      end do

      ! 3. Smval
      do i = 1, ctxinfo%NbPhasePresente
         iph = ctxinfo%NumPhasePresente(i)
         Smval(i) = &
            SmDensiteMolaire(i)*EnergieInterne(i)*inc%Saturation(iph) &
            + SmEnergieInterne(i)*DensiteMolaire(i)*inc%Saturation(iph)
      end do

   end subroutine LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv

   subroutine LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv( &
      ctxinfo, &
      DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
      PermRel, divPermRel, &
      UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
      Enthalpie, divEnthalpie, SmEnthalpie, &
      val, dval, Smval)

      type(ContextInfo), intent(in) :: ctxinfo
      double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase), &
         !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase), &
         divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
         !
         SmDensiteMolaire(NbPhase), &
         SmUnSurViscosite(NbPhase), &
         !
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncTotalPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)

      ! output
      double precision, intent(out) :: &
         val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

      ! tmp
      integer :: i, k

      val(:) = 0.d0
      dval(:, :) = 0.d0
      Smval(:) = 0.d0

      ! 1. val
      do i = 1, ctxinfo%NbPhasePresente
         val(i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)*Enthalpie(i)
      end do

      ! 2. dval
      do i = 1, ctxinfo%NbPhasePresente

         do k = 1, ctxinfo%NbIncTotalPrim

            dval(k, i) = &
               divDensiteMolaire(k, i)*PermRel(i)*UnsurViscosite(i)*Enthalpie(i) &
               + divPermRel(k, i)*DensiteMolaire(i)*UnsurViscosite(i)*Enthalpie(i) &
               + divUnsurViscosite(k, i)*DensiteMolaire(i)*PermRel(i)*Enthalpie(i) &
               + divEnthalpie(k, i)*DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
         end do
      end do

      ! 3. Smval
      do i = 1, ctxinfo%NbPhasePresente

         Smval(i) = &
            SmDensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)*Enthalpie(i) &
            + SmUnSurViscosite(i)*DensiteMolaire(i)*PermRel(i)*Enthalpie(i) &
            + SmEnthalpie(i)*DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
      end do

   end subroutine LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv

   subroutine LoisThermoHydro_divP_wellinj(NbIncLocal)

      integer, intent(in) :: NbIncLocal

      double precision :: Pws, Tw, Sw(NbPhase), Cw(NbComp)
      double precision :: &
         DensiteMolaire, dP_DensiteMolaire, &
         Viscosite, dP_Viscosite, &
         Enthalpie, dP_Enthalpie

      double precision :: dSf(NbPhase), dTf, dCf(NbComp), PermRel
      integer :: rt(IndThermique + 1)
      integer :: s, i, k

      do k = 1, NbIncLocal
         Sw(:) = 0.d0
         Sw(LIQUID_PHASE) = 1.d0

         ! node of well k
         do s = NodeDatabyWellInjLocal%Pt(k) + 1, NodeDatabyWellInjLocal%Pt(k + 1)

            Pws = PerfoWellInj(s)%Pression ! P_{w,s}

            Tw = DataWellInjLocal(k)%InjectionTemperature     ! T_w
            Cw(:) = DataWellInjLocal(k)%CompTotal(:) ! C_w

            rt = NodeRocktypeLocal(:, s)

            ! Permrel
            call f_PermRel(rt(1), LIQUID_PHASE, Sw, PermRel, dSf)

            ! Molar density
            call f_DensiteMolaire(LIQUID_PHASE, Pws, Tw, Cw, Sw, &
                                  DensiteMolaire, dP_DensiteMolaire, dTf, dCf, dSf)

            ! Viscosite
            call f_Viscosite(LIQUID_PHASE, Pws, Tw, Cw, Sw, &
                             Viscosite, dP_Viscosite, dTf, dCf, dSf)

#ifdef _THERMIQUE_
            ! Enthalpie
            call f_Enthalpie(LIQUID_PHASE, Pws, Tw, Cw, Sw, &
                             Enthalpie, dP_Enthalpie, dTf, dCf, dSf)
#endif

            do i = 1, NbComp

               ! value
               DensiteMolaireKrViscoCompWellInj(i, s) = Cw(i)*PermRel*DensiteMolaire/Viscosite

               ! div of pression
               divDensiteMolaireKrViscoCompWellInj(i, s) = &
                  Cw(i)*PermRel*(dP_DensiteMolaire/Viscosite - DensiteMolaire*dP_Viscosite/(Viscosite**2))
            end do

#ifdef _THERMIQUE_
            DensiteMolaireKrViscoEnthalpieWellInj(s) = Enthalpie*PermRel*DensiteMolaire/Viscosite

            divDensiteMolaireKrViscoEnthalpieWellInj(s) = &
               +dP_DensiteMolaire/Viscosite*Enthalpie &
               + dP_Enthalpie*DensiteMolaire/Viscosite &
               - dP_Viscosite/(Viscosite**2)*DensiteMolaire*Enthalpie
#endif
         end do ! nodes of well k
      end do ! wells

   end subroutine LoisThermoHydro_divP_wellinj

   ! allocate
   subroutine LoisThermoHydro_allocate

      integer :: nbCell, nbFrac, nbNode, nbNodeInj, max_nb_control_volumes

      nbCell = NbCellLocal_Ncpus(commRank + 1)
      nbFrac = NbFracLocal_Ncpus(commRank + 1)
      nbNode = NbNodeLocal_Ncpus(commRank + 1)
      max_nb_control_volumes = max(nbNode, nbFrac, nbCell)
      nbNodeInj = NodeByWellInjLocal%Pt(NodebyWellInjLocal%Nb + 1)
      ! print*, 'LoisThermoHydro_allocate', nbCell, nbFrac, nbNode, nbNodeInj

      ! densite massique
      allocate (DensiteMassiqueCell(NbPhase, nbCell))
      allocate (DensiteMassiqueFrac(NbPhase, nbFrac))
      allocate (DensiteMassiqueNode(NbPhase, nbNode))

      allocate (divDensiteMassiqueCell(NbIncTotalPrimMax, NbPhase, nbCell))
      allocate (divDensiteMassiqueFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
      allocate (divDensiteMassiqueNode(NbIncTotalPrimMax, NbPhase, nbNode))

      allocate (SmDensiteMassiqueCell(NbPhase, nbCell))
      allocate (SmDensiteMassiqueFrac(NbPhase, nbFrac))
      allocate (SmDensiteMassiqueNode(NbPhase, nbNode))

      ! pression
      allocate (divPressionCell(NbIncTotalPrimMax, nbCell))
      allocate (divPressionFrac(NbIncTotalPrimMax, nbFrac))
      allocate (divPressionNode(NbIncTotalPrimMax, nbNode))

      allocate (SmPressionCell(nbCell))
      allocate (SmPressionFrac(nbFrac))
      allocate (SmPressionNode(nbNode))

      ! pression capillaire
      allocate (PressionCapCell(NbPhase, nbCell))
      allocate (PressionCapFrac(NbPhase, nbFrac))
      allocate (PressionCapNode(NbPhase, nbNode))

      allocate (divPressionCapCell(NbIncTotalPrimMax, NbPhase, nbCell))
      allocate (divPressionCapFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
      allocate (divPressionCapNode(NbIncTotalPrimMax, NbPhase, nbNode))

      ! Saturation
      allocate (divSaturationCell(NbIncTotalPrimMax, NbPhase, nbCell))
      allocate (divSaturationFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
      allocate (divSaturationNode(NbIncTotalPrimMax, NbPhase, nbNode))
#ifdef _WIP_FREEFLOW_STRUCTURES_
      ! Freeflow phase molar flowrates
      allocate (divFreeFlowMolarFlowrateNode(NbIncTotalPrimMax, NbPhase, nbNode))
      allocate (SmFreeFlowMolarFlowrateNode(NbPhase, nbNode))

      ! Freeflow phase molar flowrates * Comp
      allocate (FreeFlowMolarFlowrateCompNode(NbComp, NbPhase, nbNode))
      allocate (divFreeFlowMolarFlowrateCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))
      allocate (SmFreeFlowMolarFlowrateCompNode(NbComp, NbPhase, nbNode))

      ! FIXME: Hm * (Comp - atm_comp) if gas ; 0. otherwise
      allocate (FreeFlowHmCompNode(NbComp, NbPhase, nbNode))
      allocate (divFreeFlowHmCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))
      allocate (SmFreeFlowHmCompNode(NbComp, NbPhase, nbNode))

      ! if gas-> FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)
      !    liquid-> FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)
      allocate (FreeFlowMolarFlowrateEnthalpieNode(NbPhase, nbNode))
      allocate (divFreeFlowMolarFlowrateEnthalpieNode(NbIncTotalPrimMax, NbPhase, nbNode))
      allocate (SmFreeFlowMolarFlowrateEnthalpieNode(NbPhase, nbNode))

      ! HT * (T - atm_temperature) + net Radiation (which is a factor of T**4)
      allocate (FreeFlowHTTemperatureNetRadiationNode(nbNode))
      allocate (divFreeFlowHTTemperatureNetRadiationNode(NbIncTotalPrimMax, nbNode))
      allocate (SmFreeFlowHTTemperatureNetRadiationNode(nbNode))

      ! Atmospheric Enthalpy: OF THE WATER if gas, global Enthalpy if liquid
      allocate (AtmEnthalpieNode(NbPhase, nbNode))
#endif
      ! DensiteMolaire*Kr/Viscosite * Comp
      allocate (DensiteMolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
      allocate (DensiteMolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
      allocate (DensiteMolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

      allocate (divDensiteMolaireKrViscoCompCell(NbIncTotalPrimMax, NbComp, NbPhase, nbCell))
      allocate (divDensiteMolaireKrViscoCompFrac(NbIncTotalPrimMax, NbComp, NbPhase, nbFrac))
      allocate (divDensiteMolaireKrViscoCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))

      allocate (SmDensiteMolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
      allocate (SmDensiteMolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
      allocate (SmDensiteMolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

      ! DensiteMolaire * Saturation * Comp
      call MeshSchema_allocate_CompPhaseDOFFamilyArray(DensiteMolaireSatComp)
      call MeshSchema_allocate_CompPhaseDOFFamilyArray(SmDensiteMolaireSatComp)

      allocate (divDensiteMolaireSatCompCell(NbIncTotalPrimMax, NbComp, NbPhase, nbCell))
      allocate (divDensiteMolaireSatCompFrac(NbIncTotalPrimMax, NbComp, NbPhase, nbFrac))
      allocate (divDensiteMolaireSatCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))

      ! well inj
      allocate (DensiteMolaireKrViscoCompWellInj(NbComp, nbNodeInj))
      allocate (divDensiteMolaireKrViscoCompWellInj(NbComp, nbNodeInj))

      allocate (DensiteMolaireKrViscoEnthalpieWellInj(nbNodeInj))
      allocate (divDensiteMolaireKrViscoEnthalpieWellInj(nbNodeInj))

      ! the following arrays must be allocated even if there is no energy transfer
      allocate (SmTemperatureCell(nbCell))
      allocate (SmTemperatureFrac(nbFrac))
      allocate (SmTemperatureNode(nbNode))

#ifdef _THERMIQUE_

      ! temperature
      allocate (divTemperatureCell(NbIncTotalPrimMax, nbCell))
      allocate (divTemperatureFrac(NbIncTotalPrimMax, nbFrac))
      allocate (divTemperatureNode(NbIncTotalPrimMax, nbNode))

      ! DensiteMolaire * PermRel / Viscosite * Enthalpie
      allocate (DensiteMolaireKrViscoEnthalpieCell(NbPhase, nbCell))
      allocate (DensiteMolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
      allocate (DensiteMolaireKrViscoEnthalpieNode(NbPhase, nbNode))

      allocate (divDensiteMolaireKrViscoEnthalpieCell(NbIncTotalPrimMax, NbPhase, nbCell))
      allocate (divDensiteMolaireKrViscoEnthalpieFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
      allocate (divDensiteMolaireKrViscoEnthalpieNode(NbIncTotalPrimMax, NbPhase, nbNode))

      allocate (SmDensiteMolaireKrViscoEnthalpieCell(NbPhase, nbCell))
      allocate (SmDensiteMolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
      allocate (SmDensiteMolaireKrViscoEnthalpieNode(NbPhase, nbNode))

      ! densitemolaire * energieinterne * Saturation
      call MeshSchema_allocate_PhaseDOFFamilyArray(DensiteMolaireEnergieInterneSat)

      allocate (divDensiteMolaireEnergieInterneSatCell(NbIncTotalPrimMax, NbPhase, nbCell))
      allocate (divDensiteMolaireEnergieInterneSatFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
      allocate (divDensiteMolaireEnergieInterneSatNode(NbIncTotalPrimMax, NbPhase, nbNode))

      allocate (SmDensiteMolaireEnergieInterneSatCell(NbPhase, nbCell))
      allocate (SmDensiteMolaireEnergieInterneSatFrac(NbPhase, nbFrac))
      allocate (SmDensiteMolaireEnergieInterneSatNode(NbPhase, nbNode))

#endif

      ! Work arrays
      allocate (wa_UnsurViscosite(NbPhase, max_nb_control_volumes))
      allocate (wa_divUnsurViscosite(NbIncTotalPrimMax, NbPhase, max_nb_control_volumes))
      allocate (wa_SmUnsurViscosite(NbPhase, max_nb_control_volumes))

   end subroutine LoisThermoHydro_allocate

   subroutine LoisThermoHydro_free

      ! Work arrays
      deallocate (wa_UnsurViscosite)
      deallocate (wa_divUnsurViscosite)
      deallocate (wa_SmUnsurViscosite)

      ! densite massique
      deallocate (DensiteMassiqueCell)
      deallocate (DensiteMassiqueFrac)
      deallocate (DensiteMassiqueNode)
      deallocate (divDensiteMassiqueCell)
      deallocate (divDensiteMassiqueFrac)
      deallocate (divDensiteMassiqueNode)
      deallocate (SmDensiteMassiqueCell)
      deallocate (SmDensiteMassiqueFrac)
      deallocate (SmDensiteMassiqueNode)

      ! pression
      deallocate (divPressionCell)
      deallocate (divPressionFrac)
      deallocate (divPressionNode)
      deallocate (SmPressionCell)
      deallocate (SmPressionFrac)
      deallocate (SmPressionNode)

      ! pression capillaire
      deallocate (PressionCapCell)
      deallocate (PressionCapFrac)
      deallocate (PressionCapNode)
      deallocate (divPressionCapCell)
      deallocate (divPressionCapFrac)
      deallocate (divPressionCapNode)

      ! saturation
      deallocate (divSaturationCell)
      deallocate (divSaturationFrac)
      deallocate (divSaturationNode)

#ifdef _WIP_FREEFLOW_STRUCTURES_
      ! Freeflow phase molar flowrates
      deallocate (divFreeFlowMolarFlowrateNode)
      deallocate (SmFreeFlowMolarFlowrateNode)
      ! Freeflow phase molar flowrates * Comp
      deallocate (FreeFlowMolarFlowrateCompNode)
      deallocate (divFreeFlowMolarFlowrateCompNode)
      deallocate (SmFreeFlowMolarFlowrateCompNode)
      ! Hm * (Comp - atm_comp) if gas ; 0. otherwise
      deallocate (FreeFlowHmCompNode)
      deallocate (divFreeFlowHmCompNode)
      deallocate (SmFreeFlowHmCompNode)
      ! if gas-> FreeFlowMolarFlowrate(gas) * SpecificEnthalpy(water, gas)
      !    liquid-> FreeFlowMolarFlowrate(liquid) * Enthalpie(liquid)
      deallocate (FreeFlowMolarFlowrateEnthalpieNode)
      deallocate (divFreeFlowMolarFlowrateEnthalpieNode)
      deallocate (SmFreeFlowMolarFlowrateEnthalpieNode)
      ! HT * (T-atm_temperature) + net Radiation (which is a factor of T**4)
      deallocate (FreeFlowHTTemperatureNetRadiationNode)
      deallocate (divFreeFlowHTTemperatureNetRadiationNode)
      deallocate (SmFreeFlowHTTemperatureNetRadiationNode)
      ! Atmospheric enthalpie
      deallocate (AtmEnthalpieNode)
#endif

      ! densitemolaire * Permrel / viscosite * Comp
      deallocate (DensiteMolaireKrViscoCompCell)
      deallocate (DensiteMolaireKrViscoCompFrac)
      deallocate (DensiteMolaireKrViscoCompNode)
      deallocate (divDensiteMolaireKrViscoCompCell)
      deallocate (divDensiteMolaireKrViscoCompFrac)
      deallocate (divDensiteMolaireKrViscoCompNode)
      deallocate (SmDensiteMolaireKrViscoCompCell)
      deallocate (SmDensiteMolaireKrViscoCompFrac)
      deallocate (SmDensiteMolaireKrViscoCompNode)

      ! densitemolaire * Sat * Comp
      call MeshSchema_free_CompPhaseDOFFamilyArray(DensiteMolaireSatComp)
      deallocate (divDensiteMolaireSatCompCell)
      deallocate (divDensiteMolaireSatCompFrac)
      deallocate (divDensiteMolaireSatCompNode)
      call MeshSchema_free_CompPhaseDOFFamilyArray(SmDensiteMolaireSatComp)

      ! well inj
      deallocate (DensiteMolaireKrViscoCompWellInj)
      deallocate (divDensiteMolaireKrViscoCompWellInj)

      deallocate (DensiteMolaireKrViscoEnthalpieWellInj)
      deallocate (divDensiteMolaireKrViscoEnthalpieWellInj)

      deallocate (SmTemperatureCell)
      deallocate (SmTemperatureFrac)
      deallocate (SmTemperatureNode)

#ifdef _THERMIQUE_
      ! temperature
      deallocate (divTemperatureCell)
      deallocate (divTemperatureFrac)
      deallocate (divTemperatureNode)

      ! densitemolaire * Permrel / viscosite * Enthalpie
      deallocate (DensiteMolaireKrViscoEnthalpieCell)
      deallocate (DensiteMolaireKrViscoEnthalpieFrac)
      deallocate (DensiteMolaireKrViscoEnthalpieNode)
      deallocate (divDensiteMolaireKrViscoEnthalpieCell)
      deallocate (divDensiteMolaireKrViscoEnthalpieFrac)
      deallocate (divDensiteMolaireKrViscoEnthalpieNode)
      deallocate (SmDensiteMolaireKrViscoEnthalpieCell)
      deallocate (SmDensiteMolaireKrViscoEnthalpieFrac)
      deallocate (SmDensiteMolaireKrViscoEnthalpieNode)

      ! densitemolaire * energieinterne * Saturation
      call MeshSchema_free_PhaseDOFFamilyArray(DensiteMolaireEnergieInterneSat)
      deallocate (divDensiteMolaireEnergieInterneSatCell)
      deallocate (divDensiteMolaireEnergieInterneSatFrac)
      deallocate (divDensiteMolaireEnergieInterneSatNode)
      deallocate (SmDensiteMolaireEnergieInterneSatCell)
      deallocate (SmDensiteMolaireEnergieInterneSatFrac)
      deallocate (SmDensiteMolaireEnergieInterneSatNode)

      ! ! densitemolaire * energieinterne
      ! deallocate( DensiteMolaireEnergieInterneCell)
      ! deallocate( DensiteMolaireEnergieInterneFrac)
      ! deallocate( DensiteMolaireEnergieInterneNode)
      ! deallocate( divDensiteMolaireEnergieInterneCell)
      ! deallocate( divDensiteMolaireEnergieInterneFrac)
      ! deallocate( divDensiteMolaireEnergieInterneNode)
      ! deallocate( SmDensiteMolaireEnergieInterneCell)
      ! deallocate( SmDensiteMolaireEnergieInterneFrac)
      ! deallocate( SmDensiteMolaireEnergieInterneNode)
#endif

   end subroutine LoisThermoHydro_free

end module LoisThermoHydro
