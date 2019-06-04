!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LoisThermoHydro

  use CommonMPI, only: commRank, CommonMPI_abort
  use Thermodynamics, only: &
    f_EnergieInterne, f_Enthalpie, f_Viscosite, f_DensiteMolaire, &
    f_PermRel, f_PressionCapillaire, f_DensiteMassique
  use DefModel, only: &
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
     NodeDatabyWellInjLocal, NbWellProdLocal_Ncpus, &
     CellRocktypeLocal, FracRocktypeLocal, NodeRocktypeLocal

  implicit none

  ! densite massique
  ! Rq important: it contains values for all phases, not only phase present
  double precision, allocatable, dimension(:,:), protected :: &
       DensiteMassiqueCell, &
       DensiteMassiqueFrac, &
       DensiteMassiqueNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divDensiteMassiqueCell, &
       divDensiteMassiqueFrac, &
       divDensiteMassiqueNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmDensiteMassiqueCell, &
       SmDensiteMassiqueFrac, &
       SmDensiteMassiqueNode

  ! pression
  double precision, allocatable, dimension(:,:), protected :: &
       divPressionCell, &
       divPressionFrac, &
       divPressionNode
  double precision, allocatable, dimension(:), protected :: &
       SmPressionCell, &
       SmPressionFrac, &
       SmPressionNode

  ! pression capillaire
  double precision, allocatable, dimension(:,:), protected :: &
       PressionCapCell, &
       PressionCapFrac, &
       PressionCapNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divPressionCapCell, &
       divPressionCapFrac, &
       divPressionCapNode

  ! Saturation
  double precision, allocatable, dimension(:,:,:), protected :: &
       divSaturationCell, &
       divSaturationFrac, &
       divSaturationNode

  ! Freeflow phase molar flowrates
  double precision, allocatable, dimension(:,:,:), protected :: &
       divFreeFlowMolarFlowrateNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmFreeFlowMolarFlowrateNode

  ! DensiteMolaire*Kr/Viscosite*Comp
  double precision, allocatable, dimension(:,:,:), protected :: &
       DensiteMolaireKrViscoCompCell, &
       DensiteMolaireKrViscoCompFrac, &
       DensiteMolaireKrViscoCompNode
  double precision, allocatable, dimension(:,:,:,:), protected :: &
       divDensiteMolaireKrViscoCompCell, &
       divDensiteMolaireKrViscoCompFrac, &
       divDensiteMolaireKrViscoCompNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       SmDensiteMolaireKrViscoCompCell, &
       SmDensiteMolaireKrViscoCompFrac, &
       SmDensiteMolaireKrViscoCompNode

  ! DensiteMolaire*Kr/Viscosite*Comp for wells (injection and production)
  double precision, allocatable, dimension(:,:), protected :: &
       DensiteMolaireKrViscoCompWellInj
  double precision, allocatable, dimension(:,:), protected :: &
       divDensiteMolaireKrViscoCompWellInj

  ! DensiteMolaire * Sat * Comp
  double precision, allocatable, dimension(:,:,:), protected :: &
       DensiteMolaireSatCompCell, &
       DensiteMolaireSatCompFrac, &
       DensiteMolaireSatCompNode
  double precision, allocatable, dimension(:,:,:,:), protected :: &
       divDensiteMolaireSatCompCell, &
       divDensiteMolaireSatCompFrac, &
       divDensiteMolaireSatCompNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       SmDensiteMolaireSatCompCell, &
       SmDensiteMolaireSatCompFrac, &
       SmDensiteMolaireSatCompNode

  ! temperature
  double precision, allocatable, dimension(:,:), protected :: &
       divTemperatureCell, &
       divTemperatureFrac, &
       divTemperatureNode
  double precision, allocatable, dimension(:), protected :: &
       SmTemperatureCell, &
       SmTemperatureFrac, &
       SmTemperatureNode

  ! DensiteMolaire * PermRel / Viscosite * Enthalpie
  double precision, allocatable, dimension(:,:), protected :: &
       DensiteMolaireKrViscoEnthalpieCell, &
       DensiteMolaireKrViscoEnthalpieFrac, &
       DensiteMolaireKrViscoEnthalpieNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divDensiteMolaireKrViscoEnthalpieCell, &
       divDensiteMolaireKrViscoEnthalpieFrac, &
       divDensiteMolaireKrViscoEnthalpieNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmDensiteMolaireKrViscoEnthalpieCell, &
       SmDensiteMolaireKrViscoEnthalpieFrac, &
       SmDensiteMolaireKrViscoEnthalpieNode

  ! DensiteMolaire * PermRel / Viscosite * Enthalpie for injection wells
  double precision, allocatable, dimension(:), protected :: &
       DensiteMolaireKrViscoEnthalpieWellInj
  double precision, allocatable, dimension(:), protected :: &
       divDensiteMolaireKrViscoEnthalpieWellInj


  ! densitemolaire * energieinterne * Saturation
  double precision, allocatable, dimension(:,:), protected :: &
       DensiteMolaireEnergieInterneSatCell, &
       DensiteMolaireEnergieInterneSatFrac, &
       DensiteMolaireEnergieInterneSatNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divDensiteMolaireEnergieInterneSatCell, &
       divDensiteMolaireEnergieInterneSatFrac, &
       divDensiteMolaireEnergieInterneSatNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmDensiteMolaireEnergieInterneSatCell, &
       SmDensiteMolaireEnergieInterneSatFrac, &
       SmDensiteMolaireEnergieInterneSatNode


  ! tmp values to simpfy notations of numerotation
  ! ex. NbPhasePresente = NbPhasePresente_ctx(inc%ic)
  integer, private :: &
       NbPhasePresente, NbCompCtilde, &
       NbEqFermeture, NbEqEquilibre,  &
       NbIncPTC, NbIncPTCPrim, &
       NbIncTotalPrim, & 
                                !
       NumPhasePresente(NbPhase),               &
       NumCompCtilde(NbComp),                   &
       NumCompEqEquilibre(NbEqEquilibreMax),    &
       NumIncPTC2NumIncComp_comp(NbIncPTCMax),  &
       NumIncPTC2NumIncComp_phase(NbIncPTCMax), &
                                !
       Num2PhasesEqEquilibre(2, NbEqEquilibreMax), &
       NumIncComp2NumIncPTC(NbComp, NbPhase)

  public :: &
       LoisThermoHydro_allocate,   &
       LoisThermoHydro_free,       &
       LoisThermoHydro_compute, &
       LoisThermoHydro_divP_wellinj, &
       LoisThermoHydro_divPrim_nodes


  private :: &
       LoisThermoHydro_divPrim_cv,             & ! main function for prim divs for each control volume (cv)
       LoisThermoHydro_divPrim_FreeFlow_cv, & ! prim divs for Molar flowrates in Freeflow dof
       LoisThermoHydro_init_cv,                & ! init infos according to ic (context) for each control volume (cv)
                                !
       LoisThermoHydro_fill_gradient_dfdX,     & ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       LoisThermoHydro_dfdX_ps,                & ! fill dfdX_prim/dfdX_secd with the derivatives w.r.t. the primary/secondary unknowns 
       LoisThermoHydro_DensiteMolaire_cv,      & ! prim divs: densitemolaire
       LoisThermoHydro_Viscosite_cv,           & !            1/viscosite
       LoisThermoHydro_PermRel_cv,             & !            Permrel
       LoisThermoHydro_Inc_cv,                 & !            called with Pression and Temperature
       LoisThermoHydro_PressionCapillaire_cv,  & !            Pressioncapillaire
       LoisThermoHydro_Saturation_cv,          & !            Saturation
                                !
       LoisThermoHydro_DensiteMassique_cv,       & !          densitemassique
       LoisThermoHydro_DensiteMolaireKrViscoComp_cv,     & !  densitemolaire * Permrel / viscosite * Comp
       LoisThermoHydro_DensiteMolaireSatComp_cv,         & !  densitemolaire * Saturation * Comp
       LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv, & !  densitemolaire * Permrel / viscosite * Enthalpie
       LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv  !  densitemolaire * energieinterne * Saturation

#ifdef _THERMIQUE_

  private :: &
       LoisThermoHydro_EnergieInterne_cv,  & !  Enthalpie
       LoisThermoHydro_Enthalpie_cv
#endif

contains

  ! main subroutine of this module
  ! compute all prim div for future use
  ! loop of cell/frac/node
  ! subroutine Loisthermohydro_divPrim_cv compute all prim div for one control volume
  subroutine LoisThermoHydro_compute() &
        bind(C, name="LoisThermoHydro_compute")

    ! cell
     call LoisThermoHydro_divPrim_cv(NbCellLocal_Ncpus(commRank+1), IncCell,&
          CellRocktypeLocal, &
                                !
          dXssurdXpCell, &
          SmdXsCell, &
          SmFCell,   &
                                !
          NumIncTotalPrimCell,  &
          NumIncTotalSecondCell, &
                                !
          DensiteMassiqueCell,       &
          divDensiteMassiqueCell,  &
          SmDensiteMassiqueCell,     &
                                !
          divPressionCell,         &
          SmPressionCell,            &
                                !
          divTemperatureCell,         &
          SmTemperatureCell,            &
                                !
          divSaturationCell,       &
                                !
          PressionCapCell,           &
          divPressionCapCell,      &
                                !
          DensiteMolaireSatCompCell,      &
          divDensiteMolaireSatCompCell, &
          SmDensiteMolaireSatCompCell,    &
                                !
          DensiteMolaireKrViscoCompCell,      &
          divDensiteMolaireKrViscoCompCell, &
          SmDensiteMolaireKrViscoCompCell,    &
                                !
          DensiteMolaireEnergieInterneSatCell,      &
          divDensiteMolaireEnergieInterneSatCell, &
          SmDensiteMolaireEnergieInterneSatCell,    &
                                !
          DensiteMolaireKrViscoEnthalpieCell,      &
          divDensiteMolaireKrViscoEnthalpieCell, &
          SmDensiteMolaireKrViscoEnthalpieCell )

    ! frac
     call LoisThermoHydro_divPrim_cv(NbFracLocal_Ncpus(commRank+1), IncFrac,&
          FracRocktypeLocal, &
                                !
          dXssurdXpFrac, &
          SmdXsFrac, &
          SmFFrac,   &
                                !
          NumIncTotalPrimFrac,  &
          NumIncTotalSecondFrac, &
                                !
          DensiteMassiqueFrac,       &
          divDensiteMassiqueFrac,  &
          SmDensiteMassiqueFrac,     &
                                !
          divPressionFrac,         &
          SmPressionFrac,            &
                                !
          divTemperatureFrac,         &
          SmTemperatureFrac,            &
                                !
          divSaturationFrac,       &
                                !
          PressionCapFrac,           &
          divPressionCapFrac,      &
                                !
          DensiteMolaireSatCompFrac,      &
          divDensiteMolaireSatCompFrac, &
          SmDensiteMolaireSatCompFrac,    &
                                !
          DensiteMolaireKrViscoCompFrac,      &
          divDensiteMolaireKrViscoCompFrac, &
          SmDensiteMolaireKrViscoCompFrac,    &
                                !
          DensiteMolaireEnergieInterneSatFrac,      &
          divDensiteMolaireEnergieInterneSatFrac, &
          SmDensiteMolaireEnergieInterneSatFrac,    &
                                !
          DensiteMolaireKrViscoEnthalpieFrac,      &
          divDensiteMolaireKrViscoEnthalpieFrac, &
          SmDensiteMolaireKrViscoEnthalpieFrac )

    ! node
     call LoisThermoHydro_divPrim_cv(NbNodeLocal_Ncpus(commRank+1), IncNode,&
          NodeRocktypeLocal, &
                                !
          dXssurdXpNode, &
          SmdXsNode, &
          SmFNode,   &
                                !
          NumIncTotalPrimNode,  &
          NumIncTotalSecondNode, &
                                !
          DensiteMassiqueNode,       &
          divDensiteMassiqueNode,  &
          SmDensiteMassiqueNode,     &
                                !
          divPressionNode,         &
          SmPressionNode,            &
                                !
          divTemperatureNode,         &
          SmTemperatureNode,            &
                                !
          divSaturationNode,       &
                                !
          PressionCapNode,           &
          divPressionCapNode,      &
                                !
          DensiteMolaireSatCompNode,      &
          divDensiteMolaireSatCompNode, &
          SmDensiteMolaireSatCompNode,    &
                                !
          DensiteMolaireKrViscoCompNode,      &
          divDensiteMolaireKrViscoCompNode, &
          SmDensiteMolaireKrViscoCompNode,    &
                                !
          DensiteMolaireEnergieInterneSatNode,      &
          divDensiteMolaireEnergieInterneSatNode, &
          SmDensiteMolaireEnergieInterneSatNode,    &
                                !
          DensiteMolaireKrViscoEnthalpieNode,      &
          divDensiteMolaireKrViscoEnthalpieNode, &
          SmDensiteMolaireKrViscoEnthalpieNode )

#ifdef _WIP_FREEFLOW_STRUCTURES_
     ! FreeFlow nodes
     call LoisThermoHydro_divPrim_FreeFlow_cv(NbNodeLocal_Ncpus(commRank+1), IncNode,&
          NodeRocktypeLocal, &
                                !
          dXssurdXpNode, &
          SmdXsNode, &
          SmFNode,   &
                                !
          NumIncTotalPrimNode,  &
          NumIncTotalSecondNode, &
                                !
          divFreeFlowMolarFlowrateNode, &
          SmFreeFlowMolarFlowrateNode )
#endif

    ! well injection
    !   compute q_{w,s,i} (and directive of pression)
    !   for each local well and its nodes
     call LoisThermoHydro_divP_wellinj(NbWellInjLocal_Ncpus(commRank+1))

    ! SmDensiteMolaireKrViscoEnthalpieNode = 0.d0
    ! SmDensiteMolaireEnergieInterneSatNode = 0.d0
    ! SmDensiteMolaireKrViscoCompNode = 0.d0
    ! SmDensiteMolaireSatCompNode = 0.d0

    ! SmDensiteMolaireKrViscoEnthalpieFrac = 0.d0
    ! SmDensiteMolaireEnergieInterneSatFrac = 0.d0
    ! SmDensiteMolaireKrViscoCompFrac = 0.d0
    ! SmDensiteMolaireSatCompFrac = 0.d0

    ! SmDensiteMolaireKrViscoEnthalpieCell = 0.d0
    ! SmDensiteMolaireEnergieInterneSatCell = 0.d0
    ! SmDensiteMolaireKrViscoCompCell = 0.d0
    ! SmDensiteMolaireSatCompCell = 0.d0

    ! SmDensiteMassiqueCell = 0.d0
    ! SmDensiteMassiqueNode = 0.d0
    ! SmDensiteMassiqueFrac = 0.d0

  end subroutine LoisThermoHydro_compute


  ! all operations for one cv
  subroutine LoisThermoHydro_divPrim_cv(NbIncLocal, inc,&
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
       DensiteMolaireKrViscoEnthalpie,  divDensiteMolaireKrViscoEnthalpie,  SmDensiteMolaireKrViscoEnthalpie)

    ! input
    integer, intent(in) :: NbIncLocal

    type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

    integer, intent(in) :: rt (IndThermique+1, NbIncLocal)

    integer, intent(in) ::  &
         NumIncTotalPrimCV (NbIncTotalPrimMax, NbIncLocal),  &
         NumIncTotalSecondCV (NbEqFermetureMax, NbIncLocal)

    double precision, intent(in) :: &
         dXssurdXp (NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs (NbEqFermetureMax, NbIncLocal), &
         SmF (NbEqFermetureMax, NbIncLocal)

    ! output

    double precision, intent(out) :: &
                                !
         DensiteMassique (NbPhase, NbIncLocal), &
         divDensiteMassique (NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMassique (NbPhase, NbIncLocal), &
                                !
         divPression( NbIncTotalPrimMax, NbIncLocal), &
         SmPression(NbIncLocal), &
                                !
         divTemperature( NbIncTotalPrimMax, NbIncLocal), &
         SmTemperature(NbIncLocal), &
                                !
         divSaturation ( NbIncTotalPrimMax, NbPhase, NbIncLocal), &
                                !
         PressionCap (NbPhase, NbIncLocal), &
         divPressionCap (NbIncTotalPrimMax, NbPhase, NbIncLocal),&
                                !
         DensiteMolaireSatComp (NbComp, NbPhase, NbIncLocal), &
         divDensiteMolaireSatComp (NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmDensiteMolaireSatComp (NbComp, NbPhase, NbIncLocal), &
                                !
         DensiteMolaireKrViscoComp (NbComp, NbPhase, NbIncLocal), &
         divDensiteMolaireKrViscoComp (NbIncTotalPrimMax, NbComp, NbPhase, NbIncLocal), &
         SmDensiteMolaireKrViscoComp (NbComp, NbPhase, NbIncLocal), &
                                !
         DensiteMolaireEnergieInterneSat (NbPhase, NbIncLocal), &
         divDensiteMolaireEnergieInterneSat (NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMolaireEnergieInterneSat (NbPhase, NbIncLocal), &
                                !
         DensiteMolaireKrViscoEnthalpie (NbPhase, NbIncLocal), &
         divDensiteMolaireKrViscoEnthalpie (NbIncTotalPrimMax, NbPhase, NbIncLocal), &
         SmDensiteMolaireKrViscoEnthalpie (NbPhase, NbIncLocal)

    ! tmp
    double precision :: &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase), &
                                !
         UnsurViscosite(NbPhase), &
         divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
         SmUnsurViscosite(NbPhase), &
                                !
         PermRel(NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase)


    double precision :: &
         EnergieInterne(NbPhase), &
         divEnergieInterne(NbIncTotalPrimMax, NbPhase), &
         SmEnergieInterne(NbPhase), &
                                !
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncTotalPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)

    integer :: k

  ! loop over each local element
  do k=1, NbIncLocal

    ! init tmp values for each cv
    call LoisThermoHydro_init_cv(inc(k))

    ! viscosite
    call LoisThermoHydro_viscosite_cv(inc(k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite)

    ! densite massique
    call LoisThermoHydro_densitemassique_cv(inc(k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         DensiteMassique(:,k), divDensiteMassique(:,:,k), SmDensiteMassique(:,k))

    ! deniste molaire
    call LoisThermoHydro_densitemolaire_cv(inc(k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire)

    ! PermRel
    call LoisThermoHydro_PermRel_cv(inc(k), rt(:,k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         PermRel, divPermRel)

    ! Pression (unknown index is 1)
    call LoisThermoHydro_Inc_cv(1, inc(k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         dXssurdXp(:,:,k), SmdXs(:,k),  &
         divPression(:,k), SmPression(k))

    ! Pression Capillaire
    call LoisThermoHydro_PressionCapillaire_cv(rt(:,k), inc(k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         dXssurdXp(:,:,k), &
         PressionCap(:,k), divPressionCap(:,:,k))

    ! Saturation div
    call LoisThermoHydro_Saturation_cv(inc(k), &
      NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
      dXssurdXp(:,:,k), &
      divSaturation(:,:,k))

    ! term: DensiteMolaire * PermRel / Viscosite * Comp
    call LoisThermoHydro_DensiteMolaireKrViscoComp_cv( inc(k), &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
         PermRel, divPermRel, &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k),  &
         dXssurdXp(:,:,k), SmdXs(:,k), &
         DensiteMolaireKrViscoComp(:,:,k), divDensiteMolaireKrViscoComp(:,:,:,k), SmDensiteMolaireKrViscoComp(:,:,k))

    ! term: DensiteMolaire * Saturation * Comp
    call LoisThermoHydro_DensiteMolaireSatComp_cv( &
         inc(k), divSaturation(:,:,k), &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k),  &
         dXssurdXp(:,:,k), SmdXs(:,k), &
         DensiteMolaireSatComp(:,:,k), divDensiteMolaireSatComp(:,:,:,k), SmDensiteMolaireSatComp(:,:,k))


#ifdef _THERMIQUE_

    ! FIXME: Temperature (unknown index is 2)
    call LoisThermoHydro_Inc_cv(2, inc(k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         dXssurdXp(:,:,k), SmdXs(:,k),  &
         divTemperature(:,k), SmTemperature(k))

    ! energie interne
    call LoisThermoHydro_EnergieInterne_cv(inc(k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         EnergieInterne, divEnergieInterne, SmEnergieInterne)

    ! Enthalpie
    call LoisThermoHydro_Enthalpie_cv(inc(k), dXssurdXp(:,:,k), SmdXs(:,k), &
         NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
         Enthalpie, divEnthalpie, SmEnthalpie)

    ! term: DensiteMolaire * Energieinterne * Saturation
    call LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv(  &
         inc(k), divSaturation(:,:,k), &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
         EnergieInterne, divEnergieInterne, SmEnergieInterne, &
         DensiteMolaireEnergieInterneSat(:,k), divDensiteMolaireEnergieInterneSat(:,:,k), SmDensiteMolaireEnergieInterneSat(:,k))

    ! term: DensiteMolaire * PermRel / Viscosite * Enthalpie
    call LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv( &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
         PermRel, divPermRel, &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
         Enthalpie, divEnthalpie, SmEnthalpie, &
         DensiteMolaireKrViscoEnthalpie(:,k), divDensiteMolaireKrViscoEnthalpie(:,:,k), SmDensiteMolaireKrViscoEnthalpie(:,k))
#endif

 !       do i=1, NbPhasePresente_ctx(inc(k)%ic)
 !    write(*,*)"---------------------------------------------------------------------------------"
 !    write(*,*)"print laws for top k which context is ",inc(k)%ic
 !    write(*,*)"---------------------------------------------------------------------------------"
 !    write(*,*)"phase",NumPhasePresente(i)
 !    write(*,*)"---------------------------------------------------------------------------------"
 !    print*,"DensiteMolaire"
 !    print*, DensiteMolaire(i)
 !    print*, ""
!
 !    do j=1, NbIncTotalPrim
 !       print*, i, j, divDensiteMolaire(j,i)
 !    end do
 !    print*, ""
 !    print*, SmDensiteMolaire(i)
 !    print*, ""
!
 !    print*,"1.d0/UnsurViscosite(i)"
 !    print*, 1.d0/UnsurViscosite(i)
!
 !    print*, "Enthalpie"
 !    print*, Enthalpie(i)
 !    do j=1, NbIncTotalPrim
 !         print*, i, j, divEnthalpie(j,i)
 !    end do
 !    print*, ""
!
 !    print*, "PermRel"
 !    print*, PermRel(i)
 !    do j=1, NbIncTotalPrim
 !       print*, i, j, divPermRel(j,i)
 !    end do
 !    print*, ""
!
 !    print*, "DensiteMolaireKrViscoComp"
 !    print*, "all comp",DensiteMolaireKrViscoComp(:,i,k)
 !    print*, ""
!
 !    do j=1, NbIncTotalPrim
 !       print*, "all comp",i,j,divDensiteMolaireKrViscoComp(j,:,i,k)
 !       ! print*, divDensiteMolaire
 !    end do
 !    print*, ""
 !    print*, SmDensiteMolaireKrViscoComp(:,i,k)
 !    print*, ""
!
 !    write(*,*), "DensiteMolaireKrViscoEnthalpie"
 !    do j=1, NbIncTotalPrim
 !       write(*,*), i,j,divDensiteMolaireKrViscoEnthalpie(j,i,k)
 !       ! print*, divDensiteMolaire
 !    end do
 !    print*, ""
 !    write(*,*), SmDensiteMolaireKrViscoEnthalpie(i,k)
 !    print*, ""
!
 !    write(*,*), "DensiteMolaireEnergieInterneSat"
 !    do j=1, NbIncTotalPrim
 !       write(*,*), i,j,divDensiteMolaireEnergieInterneSat(j,i,k)
 !    end do
 !    print*, ""
 !    write(*,*), SmDensiteMolaireEnergieInterneSat(i,k)
 !    print*, ""
 !  end do
!
 !    write(*,*), "Temperature"
 !    do j=1, NbIncTotalPrim
 !       write(*,*), "same for all phase",j,divTemperature(j,k)
 !    end do
 !    print*, ""
 !    print*, SmTemperature(k)
 !    print*, ""
!
  end do

  end subroutine LoisThermoHydro_divPrim_cv

  ! Compute derivative of phase molar flowrates in the Freeflow dof
  subroutine LoisThermoHydro_divPrim_FreeFlow_cv(NbIncLocal, inc,&
       rt, &
       dXssurdXp, SmdXs, SmF, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, &
       divFreeFlowMolarFlowrate, SmFreeFlowMolarFlowrate)

    ! input
    integer, intent(in) :: NbIncLocal

    type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

    integer, intent(in) :: rt (IndThermique+1, NbIncLocal)

    integer, intent(in) ::  &
         NumIncTotalPrimCV (NbIncTotalPrimMax, NbIncLocal),  &
         NumIncTotalSecondCV (NbEqFermetureMax, NbIncLocal)

    double precision, intent(in) :: &
         dXssurdXp (NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs (NbEqFermetureMax, NbIncLocal), &
         SmF (NbEqFermetureMax, NbIncLocal)

    ! output
    double precision, intent(out) :: & 
       divFreeFlowMolarFlowrate(NbIncTotalPrimMax, NbPhase, NbIncLocal), &
       SmFreeFlowMolarFlowrate(NbPhase, NbIncLocal)

    integer :: k, i

  ! loop over each local element
  do k=1, NbIncLocal
     
     if(inc(k)%ic<2**NbPhase) cycle ! FIXME: loop over Freeflow dof only, avoid reservoir context

    ! init tmp values for each cv
     call LoisThermoHydro_init_cv(inc(k))

    do i=1, NbPhasePresente
       ! phase molar flowrate (unknown index is NbIncPTC+NbPhasePresente-1+i) (FIXME: -1 because one saturation is eliminated)
       call LoisThermoHydro_Inc_cv(NbIncPTC+NbPhasePresente-1+i, inc(k), &
            NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
            dXssurdXp(:,:,k), SmdXs(:,k),  &
            divFreeFlowMolarFlowrate(:,i,k), SmFreeFlowMolarFlowrate(i,k))
    enddo

  end do

  end subroutine LoisThermoHydro_divPrim_FreeFlow_cv


  !> Update thermo Laws of nodes
  subroutine LoisThermoHydro_divPrim_nodes

     call LoisThermoHydro_divPrim_cv(NbNodeLocal_Ncpus(commRank+1), IncNode,&
          NodeRocktypeLocal, &
                                !
          dXssurdXpNode, &
          SmdXsNode, &
          SmFNode,   &
                                !
          NumIncTotalPrimNode,  &
          NumIncTotalSecondNode, &
                                !
          DensiteMassiqueNode,       &
          divDensiteMassiqueNode,  &
          SmDensiteMassiqueNode,     &
                                !
          divPressionNode,         &
          SmPressionNode,            &
                                !
          divTemperatureNode,         &
          SmTemperatureNode,            &
                                !
          divSaturationNode,       &
                                !
          PressionCapNode,           &
          divPressionCapNode,      &
                                !
          DensiteMolaireSatCompNode,      &
          divDensiteMolaireSatCompNode, &
          SmDensiteMolaireSatCompNode,    &
                                !
          DensiteMolaireKrViscoCompNode,      &
          divDensiteMolaireKrViscoCompNode, &
          SmDensiteMolaireKrViscoCompNode,    &
                                !
          DensiteMolaireEnergieInterneSatNode,      &
          divDensiteMolaireEnergieInterneSatNode, &
          SmDensiteMolaireEnergieInterneSatNode,    &
                                !
          DensiteMolaireKrViscoEnthalpieNode,      &
          divDensiteMolaireKrViscoEnthalpieNode, &
          SmDensiteMolaireKrViscoEnthalpieNode )

  end subroutine LoisThermoHydro_divPrim_nodes


  subroutine LoisThermoHydro_init_cv(inc)

    type(TYPE_IncCVReservoir), intent(in) :: inc

    NbPhasePresente = NbPhasePresente_ctx(inc%ic)
    NbCompCtilde = NbCompCtilde_ctx(inc%ic)

    NbEqFermeture = NbEqFermeture_ctx(inc%ic)
    NbEqEquilibre = NbEqEquilibre_ctx(inc%ic)

    NbIncPTC = NbIncPTC_ctx(inc%ic)
    NbIncPTCPrim  = NbIncPTC - NbEqFermeture


    ! ps. if there is only one phase, phase is secd
    NbIncTotalPrim = NbIncTotalPrim_ctx(inc%ic)

    NumPhasePresente(:) = NumPhasePresente_ctx(:,inc%ic)
    NumCompCtilde(:) = NumCompCtilde_ctx(:,inc%ic)
    NumCompEqEquilibre(:) = NumCompEqEquilibre_ctx(:,inc%ic)
    NumIncPTC2NumIncComp_comp(:) = NumIncPTC2NumIncComp_comp_ctx(:,inc%ic)
    NumIncPTC2NumIncComp_phase(:) = NumIncPTC2NumIncComp_phase_ctx(:,inc%ic)

    Num2PhasesEqEquilibre(:,:) = Num2PhasesEqEquilibre_ctx(:,:,inc%ic)
    NumIncComp2NumIncPTC(:,:) = NumIncComp2NumIncPTC_ctx(:,:,inc%ic)

  end subroutine LoisThermoHydro_init_cv

  ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
  subroutine LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)
   integer, intent(in) :: iph ! num of phase
    double precision, intent(in) :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision, intent(out) :: dfdX(NbIncTotalMax)

    ! tmp
    integer :: j, jc, jph
    double precision :: dfS_elim

    dfdX(:) = 0.d0

    dfdX(1) = dPf  ! P

#ifdef _THERMIQUE_
    dfdX(2) = dTf
#endif

    do j=1, NbComp
         if(MCP(j,iph)==1) then
            jc = NumIncComp2NumIncPTC(j,iph)
            dfdX(jc) = dCf(j)
         end if
    enddo
    
!    Previous implementation:
!    do j=1, NbPhase ! S
!         jc = j + NbIncPTC
!         dfdX(jc) = dSf(j)
!    enddo
    dfS_elim = dSf(NumPhasePresente(NbPhasePresente))
    do j=1, NbPhasePresente-1
         jph = NumPhasePresente(j)

         jc = j + NbIncPTC
         dfdX(jc) = dSf(jph) - dfS_elim ! last saturation is eliminated
    enddo
    
  end subroutine LoisThermoHydro_fill_gradient_dfdX

  subroutine LoisThermoHydro_dfdX_ps(iph, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
       dfdX_prim, dfdX_secd)

    ! input 
    integer, intent(in) :: iph ! num of phase
    integer, intent(in) :: &
      NumIncTotalPrimCV(NbIncTotalPrimMax), &
      NumIncTotalSecondCV(NbEqFermetureMax)
    double precision, intent(in) :: dfdX(NbIncTotalMax)  ! dfdX = (df/dP, df/dT, df/dC, df/dS)
    ! output
    double precision, intent(out) :: &
      dfdX_prim(NbIncTotalPrimMax, NbPhase), &
      dfdX_secd(NbEqFermetureMax, NbPhase)

    ! tmp
      integer :: j
  
    ! prim and secd part of dfdX
    do j=1, NbIncTotalPrim
      dfdX_prim(j,iph) = dfdX(NumIncTotalPrimCV(j))
    end do

    ! secd unknowns
    do j=1, NbEqFermeture 
         dfdX_secd(j,iph) = dfdX(NumIncTotalSecondCV(j))
    end do

  end subroutine LoisThermoHydro_dfdX_ps

  subroutine LoisThermoHydro_densitemassique_cv(inc, dXssurdXp, SmdXs, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

#ifdef DEBUG_LOISTHEMOHYDRO
    do i = 1, NbPhasePresente
      iph = NumPhasePresente(i)
      write(*,*) 'Phase', iph, 'MCP row=', MCP(:,iph)
    end do
#endif

    dfdX_secd(:,:) = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_DensiteMassique(iph,inc%Pression,inc%Temperature, &
            inc%Comp(:,iph), inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(iph) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(iph, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=densitemassique
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhase, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

    ! - dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhase,  &
         -1.d0, dfdX_secd, NbEqFermetureMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_densitemassique_cv



  subroutine LoisThermoHydro_viscosite_cv(inc, dXssurdXp, SmdXs, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    dfdX_secd(:,:) = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_Viscosite(iph,inc%Pression,inc%Temperature, &
            inc%Comp(:,iph), inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

#ifndef NDEBUG
    if(f.le.0) call CommonMPI_abort('inconsistent viscosity value')
#endif

       val(i) = 1.d0/f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)
       dfdX(:) = -dfdX(:)/f**2

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(i, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    !print*, SmdXs

    ! dv/dXp - dv/dXs*dXs/dXp, v=viscosite
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order,
    ! consider all the mats as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbEqFermetureMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_viscosite_cv



  subroutine LoisThermoHydro_densitemolaire_cv(inc, dXssurdXp, SmdXs, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    dfdX_secd(:,:) = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_DensiteMolaire(iph,inc%Pression,inc%Temperature, &
            inc%Comp(:,iph), inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(i, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=densitemolaire
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbEqFermetureMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_densitemolaire_cv



  subroutine LoisThermoHydro_PermRel_cv(inc, rt, dXssurdXp, SmdXs, &
          NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
    INTEGER, INTENT(IN) :: rt(IndThermique+1)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)

    ! tmp
    double precision :: f, dSf(NbPhase), dCf(NbComp)
    double precision :: dfdX(NbIncTotalMax)
    double precision :: dfdX_secd(NbEqFermetureMax, NbPhase) !=NbEqFermetureMax
    integer :: i, iph

! previous implementation
 !  val(:) = 0.d0
 !  dval(:,:) = 0.d0
!
 !  ! if there is only one presente phase,
 !  ! this Saturation must be secd --> dval=0, Sm=0
 !  do i=1, NbPhasePresente
 !     iph = NumPhasePresente(i)
!
 !     call f_PermRel(rt, iph, inc%Saturation, f, dSf)
!
 !     val(i) = f
 !     dfS_secd = dSf( NumPhasePresente(NbPhasePresente)) ! the last is secd, FIXME: sum S=1 in hard, last saturation is eliminated
!
 !     ! alpha=1,2,...,NbPhasepresente-1
 !     ! Ps. NbIncPTCPrim+NbPhasePresente-1=NbIncTotalPrim
 !     do j=1, NbPhasePresente - 1
 !        jph = NumPhasePresente(j)
!
 !        dval(j+NbIncPTCPrim,i) = dSf(jph) - dfS_secd
 !     end do
 !  end do

    val(:) = 0.d0
    dval(:,:) = 0.d0
    dfdX_secd(:,:) = 0.d0

    dCf = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_PermRel(rt, iph, inc%Saturation, f, dSf)

       val(i) = f

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, 0.d0, 0.d0, dCf, dSf, dfdX)

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(iph, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=densitemassique
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhase, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

!
!    ! - dv/dXs*SmdXs  ! FIXME: add Sm ?
!    Smval(:) = 0.d0
!    call dgemv('T', NbEqFermeture, NbPhase,  &
!         -1.d0, dfdX_secd, NbEqFermetureMax, &
!         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_PermRel_cv

  subroutine LoisThermoHydro_Inc_cv(index_inc, inc, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, &
       dXssurdXp, SmdXs, &
       dval, Smval)

    ! input
    integer, intent(in) :: index_inc
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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

    IF(ANY(NumIncTotalPrimCV == index_inc))THEN ! index_inc (P or T) is prim
      do i=1,NbIncTotalPrimMax
        if (NumIncTotalPrimCV(i) == index_inc) then
          dval(i) = 1.d0
        endif
      enddo
    ELSE IF(ANY(NumIncTotalSecondCV == index_inc))THEN ! index_inc (P or T) is secd
      do i=1,NbEqFermeture
        if (NumIncTotalSecondCV(i) == index_inc) then
          dval(:) = - dXssurdXp(:,i)
          Smval = - SmdXs(i)
        endif
      enddo
    ELSE ! index_inc (P or T) not found
      if(index_inc == 1) call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, P not found ')
      if(index_inc == 2) call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, T not found ')
    ENDIF

  end subroutine LoisThermoHydro_Inc_cv



  subroutine LoisThermoHydro_Saturation_cv(inc, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, &
       dXssurdXp, &
       dval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
    integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV( NbEqFermetureMax)
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax)

    ! output
    double precision, intent(out) :: dval(NbIncTotalPrimMax, NbPhase)

    ! tmp
    integer :: i, k

    dval(:,:) = 0.d0

    ! FIXME: Elimination of the last present phase (sum S =1 forced in the code)
    DO i=1, NbPhasePresente - 1

      IF(ANY(NumIncTotalPrimCV == i+NbIncPTC))THEN ! S is prim
        do k=1,NbIncTotalPrimMax
          if (NumIncTotalPrimCV(k) == i+NbIncPTC) then
            dval(k,i) = 1.d0

            dval(k,NbPhasePresente) = -1.d0
          endif
        enddo
      ELSE IF(ANY(NumIncTotalSecondCV == i+NbIncPTC))THEN ! S is secd
        do k=1,NbEqFermeture
          if (NumIncTotalSecondCV(k) == i+NbIncPTC) then
            dval(:,i) = - dXssurdXp(:,k)

            dval(:,NbPhasePresente) = dXssurdXp(:,k)
          endif
        enddo
      ELSE IF(inc%ic<2**NbPhase)THEN! S not found for a reservoir dof
        call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, S not found ')
      ENDIF
    ENDDO

  end subroutine LoisThermoHydro_Saturation_cv



  subroutine LoisThermoHydro_PressionCapillaire_cv(rt, inc, &
          NumIncTotalPrimCV, NumIncTotalSecondCV, &
          dXssurdXp, &
          val, dval)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0

    do iph=1, NbPhase

       call f_PressionCapillaire(rt, iph, inc%Saturation, f, dSf)

       val(iph) = f

       dfS_secd = dSf(NumPhasePresente(NbPhasePresente))

       do j=1, NbPhasePresente - 1
          ! Look for S, is it primary or secondary unknowns ?
          ! FIXME: Elimination of the last present phase (sum S =1 forced in the code)
          IF(ANY(NumIncTotalPrimCV == j+NbIncPTC))THEN ! S is prim
               jph = NumPhasePresente(j)
               do k=1,NbIncTotalPrimMax
                 if (NumIncTotalPrimCV(k) == j+NbIncPTC) then
                   dval(k,iph) = dSf(jph) - dfS_secd
                 endif
               enddo
          ELSE IF(ANY(NumIncTotalSecondCV == j+NbIncPTC))THEN ! S is secd
               do k=1,NbEqFermeture
                 if (NumIncTotalSecondCV(k) == j+NbIncPTC) then
                    dval(:,iph) = - dXssurdXp(:,k)
call CommonMPI_abort("je ne sais pas ce qui se passe quand S secondaire")
                 endif
               enddo
          ELSE IF(inc%ic<2**NbPhase)THEN! S not found for a reservoir dof
               call CommonMPI_abort(' pb in NumIncTotal in LoisThermoHydro, S not found ')
          ENDIF
       enddo

          
    end do ! iph

  end subroutine LoisThermoHydro_PressionCapillaire_cv


#ifdef _THERMIQUE_

  subroutine LoisThermoHydro_EnergieInterne_cv(inc, dXssurdXp, SmdXs, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    dfdX_secd(:,:) = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_EnergieInterne(iph,inc%Pression,inc%Temperature, &
            inc%Comp(:,iph), inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(i, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=energieinterne
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhasePresente , NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbEqFermetureMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_EnergieInterne_cv

  subroutine LoisThermoHydro_Enthalpie_cv(inc, dXssurdXp, SmdXs, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    dfdX_secd(:,:) = 0.d0

    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_Enthalpie(iph,inc%Pression,inc%Temperature, &
            inc%Comp(:,iph), inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(i, NumIncTotalPrimCV, NumIncTotalSecondCV, dfdX, &
            dval, dfdX_secd)
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=viscosite
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncTotalPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncTotalPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncTotalPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbEqFermetureMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_Enthalpie_cv

#endif


  ! term: desitemolaire * Saturation
  subroutine LoisThermoHydro_DensiteMolaireSat_cv( &
       inc, divSaturation, &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         divSaturation(NbIncTotalPrimMax, NbPhase), &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

    ! tmp
    integer :: i, iph, k

    val(:) = 0.d0
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       val(i) = DensiteMolaire(i) * inc%Saturation(iph)
    end do

    ! 2. dval: d(DensiteMolaire * Saturation)
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       do k=1, NbIncTotalPrim

          dval(k,i) = &
               divDensiteMolaire(k,i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       Smval(i) = SmDensiteMolaire(i)*inc%Saturation(iph)
    end do

  end subroutine LoisThermoHydro_DensiteMolaireSat_cv


  ! term: desitemolaire * Saturation * Comp
  subroutine LoisThermoHydro_DensiteMolaireSatComp_cv( &
       inc, divSaturation, &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       NumIncTotalPrimCV, NumIncTotalSecondCV,  &
       dXssurdXp, SmdXs, &
       val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         divSaturation(NbIncTotalPrimMax, NbPhase), &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase)

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
    double precision :: dv

    double precision :: dvi(NbIncTotalPrimMax)

    val(:,:) = 0.d0
    dval(:,:,:) = 0.d0
    Smval(:,:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       do icp=1, NbComp

          ! only {alpha | alpha \in Q_k \cap P_i} is useful
          ! To understant better, change the order of the loop do i=.. and the loop do icp=..
          if(MCP(icp,iph)==1) then ! P_i
             val(icp,i) = DensiteMolaire(i)*inc%Saturation(iph) &
                  * inc%Comp(icp,iph)
          end if
       end do
    end do

    ! 2.1 div(DensiteMolaire * Saturation) * C_i^alpha
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       ! 2.1.1 compute dvi, tmp vector, used in 2.1.2
       do k=1, NbIncTotalPrim

          dvi(k) = &
               divDensiteMolaire(k,i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)
       end do

       ! 2.1.2
       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             do k=1, NbIncTotalPrim
                dval(k,icp,i) = dvi(k) * inc%Comp(icp,iph)
             end do

          end if
       end do ! end of P_i

    end do ! end of 2.1

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       dv = SmDensiteMolaire(i)*inc%Saturation(iph)

       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             Smval(icp,i) = dv * inc%Comp(icp,iph)
          end if
       end do
    end do ! end of 2.2

    ! 3.1 + DensiteMolaire*Saturation * div(C_i^alpha)
    ! 3.2 + DensiteMolaire*Saturation * Sm(C_i^alpha)

    ! loop of inc prim
    do j=1, NbIncTotalPrim
       numj = NumIncTotalPrimCV(j)

       ! numi is C_i^alpha in (P,T,C,S)
       if( (numj>(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          ! for j, only div(C_jcp^jph) is not zero
          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then
                dval(j,jcp,i) = dval(j,jcp,i) &
                     + DensiteMolaire(i)*inc%Saturation(iph)
             end if
          end do

       end if
    end do ! end of inc prim

    ! loop of secd unklnowns
    do j=1, NbEqFermeture 
       numj = NumIncTotalSecondCV(j)

       ! numi is C in (P,T,C)
       if( (numj>(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then

                dv = DensiteMolaire(i)*inc%Saturation(iph)

                do s=1, NbIncTotalPrim
                   dval(s,jcp,i) = dval(s,jcp,i) - dXssurdXp(s,j) * dv
                end do

                Smval(jcp,i) = Smval(jcp,i) - SmdXs(j) * dv
             end if
          end do

       end if
    end do ! end of inc secd

  end subroutine LoisThermoHydro_DensiteMolaireSatComp_cv


  ! div and Sm of term: DensiteMolaire*Kr/Viscosite
  subroutine LoisThermoHydro_DensiteMolaireKrVisco_cv( &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel,             &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       val, dval, Smval)

    ! input
    double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase), &
                                !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase), &
         divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
                                !
         SmDensiteMolaire(NbPhase),&
         SmUnSurViscosite(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

    ! tmp
    integer :: i, k

    val(:) = 0.d0
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       val(i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
    end do

    ! 2. dval: d(DensiteMolaire*Kr/Viscosite)
    do i=1, NbPhasePresente

       do k=1, NbIncTotalPrim

          dval(k,i) = &
               divDensiteMolaire(k,i)*PermRel(i)*UnsurViscosite(i)    &
               + divPermRel(k,i)*DensiteMolaire(i)*UnsurViscosite(i)  &
               + divUnsurViscosite(k,i)*DensiteMolaire(i)*PermRel(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente

       Smval(i) = &
            SmDensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
            + SmUnSurViscosite(i)*DensiteMolaire(i)*PermRel(i)
    end do

  end subroutine LoisThermoHydro_DensiteMolaireKrVisco_cv


  ! div and Sm of term: DensiteMolaire*Kr/Viscosite*Comp
  subroutine LoisThermoHydro_DensiteMolaireKrViscoComp_cv(inc, &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel,             &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       NumIncTotalPrimCV, NumIncTotalSecondCV,  &
       dXssurdXp, SmdXs, &
       val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase),   &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase),   &
                                !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divPermRel(NbIncTotalPrimMax, NbPhase), &
         divUnsurViscosite(NbIncTotalPrimMax, NbPhase), &
                                !
         SmDensiteMolaire(NbPhase), &
         SmUnSurViscosite(NbPhase)

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
    double precision :: dv

    double precision :: dvi(NbIncTotalPrimMax)

    val(:,:) = 0.d0
    dval(:,:,:) = 0.d0
    Smval(:,:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       do icp=1, NbComp

          ! only {alpha | alpha \in Q_k \cap P_i} is useful
          ! To understand better, change the order of the loop do i=.. and the loop do icp=..
          if(MCP(icp,iph)==1) then ! P_i
             val(icp,i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
                  * inc%Comp(icp,iph)
          end if
       end do
    end do

    ! 2.1 div(DensiteMolaire*Kr/Viscosite)*C_i^alpha
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       ! 2.1.1 compute dvi, tmp vector, used in 2.1.2
       do k=1, NbIncTotalPrim
          dvi(k) = &
               divDensiteMolaire(k,i)*PermRel(i)*UnsurViscosite(i)    &
               + divPermRel(k,i)*DensiteMolaire(i)*UnsurViscosite(i)  &
               + divUnsurViscosite(k,i)*DensiteMolaire(i)*PermRel(i)
       end do

       ! 2.1.2
       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             do k=1, NbIncTotalPrim
                dval(k,icp,i) = dvi(k) * inc%Comp(icp,iph)
             end do

          end if
       end do ! end of P_i

    end do ! end of 2.1

    ! 2.2. Sm(DensiteMolaire*Kr/Viscosite)*C_i^alpha
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       dv = &
            SmDensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
            + SmUnSurViscosite(i)*DensiteMolaire(i)*PermRel(i)

       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             Smval(icp,i) = dv * inc%Comp(icp,iph)
          end if
       end do
    end do ! end of 2.2

    ! 3.1 + DensiteMolaire*Kr/Viscosite * div(C_i^alpha)
    ! 3.2 + DensiteMolaire*Kr/Viscosite * Sm(C_i^alpha)

    ! loop of inc prim
    do j=1, NbIncTotalPrim
       numj = NumIncTotalPrimCV(j)

       ! numi is C_i^alpha in (P,T,C,S)
       if( (numj>(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          ! for j, only div(C_jcp^jph) is not zero
          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then
                dval(j,jcp,i) = dval(j,jcp,i) &
                     + DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
             end if
          end do

       end if
    end do ! end of inc prim

    ! secd unknowns
    do j=1, NbEqFermeture 
       numj = NumIncTotalSecondCV(j)

       ! numi is C in (P,T,C)
       if( (numj>(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then

                dv = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)

                do s=1, NbIncTotalPrim
                   dval(s,jcp,i) = dval(s,jcp,i) - dXssurdXp(s,j) * dv
                end do

                Smval(jcp,i) = Smval(jcp,i) - SmdXs(j) * dv
             end if
          end do

       end if
    end do ! end of inc secd

  end subroutine LoisThermoHydro_DensiteMolaireKrViscoComp_cv


  subroutine LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv( &
       inc, divSaturation, &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       EnergieInterne, divEnergieInterne, SmEnergieInterne, &
       val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         EnergieInterne(NbPhase), &
                                !
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         divEnergieInterne(NbIncTotalPrimMax, NbPhase), &
         divSaturation(NbIncTotalPrimMax, NbPhase),     &
                                !
         SmDensiteMolaire(NbPhase),&
         SmEnergieInterne(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncTotalPrimMax, NbPhase), Smval(NbPhase)

    ! tmp
    integer :: i, iph, k

    val(:) = 0.d0
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       val(i) = DensiteMolaire(i)*EnergieInterne(i)*inc%Saturation(iph)
    end do

    ! 2. dval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       do k=1, NbIncTotalPrim

          dval(k,i) = &
               divDensiteMolaire(k,i)*EnergieInterne(i)*inc%Saturation(iph) &
               + divEnergieInterne(k,i)*DensiteMolaire(i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)*EnergieInterne(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       Smval(i) = &
            SmDensiteMolaire(i)*EnergieInterne(i)*inc%Saturation(iph) &
            + SmEnergieInterne(i)*DensiteMolaire(i)*inc%Saturation(iph)
    end do

  end subroutine LoisThermoHydro_DensiteMolaireEnergieInterneSat_cv


  subroutine LoisThermoHydro_DensiteMolaireKrViscoEnthalpie_cv( &
       DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel, &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       Enthalpie, divEnthalpie, SmEnthalpie, &
       val, dval, Smval)

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase),   &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase),   &
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
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       val(i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)*Enthalpie(i)
    end do

    ! 2. dval
    do i=1, NbPhasePresente

       do k=1, NbIncTotalPrim

          dval(k,i) = &
               divDensiteMolaire(k,i)*PermRel(i)*UnsurViscosite(i)*Enthalpie(i)    &
               + divPermRel(k,i)*DensiteMolaire(i)*UnsurViscosite(i)*Enthalpie(i)  &
               + divUnsurViscosite(k,i)*DensiteMolaire(i)*PermRel(i)*Enthalpie(i)  &
               + divEnthalpie(k,i)*DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente

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
    integer :: rt(IndThermique+1)
    integer :: s, i, k

do k=1, NbIncLocal
    Sw(:) = 0.d0
    Sw(LIQUID_PHASE) = 1.d0

    ! node of well k
    do s=NodeDatabyWellInjLocal%Pt(k)+1, NodeDatabyWellInjLocal%Pt(k+1)

       Pws = PerfoWellInj(s)%Pression ! P_{w,s}

       Tw = DataWellInjLocal(k)%Temperature     ! T_w
       Cw(:) = DataWellInjLocal(k)%CompTotal(:) ! C_w

       rt = NodeRocktypeLocal(:,s)

       ! Permrel
       call f_PermRel(rt,LIQUID_PHASE, Sw, PermRel, dSf)

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

       do i=1, NbComp

          ! value
          DensiteMolaireKrViscoCompWellInj(i,s) = Cw(i) * PermRel * DensiteMolaire / Viscosite

          ! div of pression
          divDensiteMolaireKrViscoCompWellInj(i,s) = Cw(i) * PermRel  &
               * (dP_DensiteMolaire / Viscosite - DensiteMolaire * dP_Viscosite / (Viscosite**2) )
       end do

#ifdef _THERMIQUE_
       DensiteMolaireKrViscoEnthalpieWellInj(s) = Enthalpie * PermRel * DensiteMolaire / Viscosite

       divDensiteMolaireKrViscoEnthalpieWellInj(s) = &
            + dP_DensiteMolaire / Viscosite * Enthalpie &
            + dP_Enthalpie * DensiteMolaire / Viscosite &
            - dP_Viscosite / (Viscosite**2) * DensiteMolaire * Enthalpie
#endif
    end do ! nodes of well k
end do ! wells

  end subroutine LoisThermoHydro_divP_wellinj


  ! allocate
  subroutine LoisThermoHydro_allocate

    integer :: nbCell, nbFrac, nbNode, nbNodeInj

    nbCell = NbCellLocal_Ncpus(commRank+1)
    nbFrac = NbFracLocal_Ncpus(commRank+1)
    nbNode = NbNodeLocal_Ncpus(commRank+1)
    nbNodeInj = NodeByWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)
    ! print*, 'LoisThermoHydro_allocate', nbCell, nbFrac, nbNode, nbNodeInj

    ! densite massique
    allocate( DensiteMassiqueCell(NbPhase, nbCell))
    allocate( DensiteMassiqueFrac(NbPhase, nbFrac))
    allocate( DensiteMassiqueNode(NbPhase, nbNode))

    allocate( divDensiteMassiqueCell(NbIncTotalPrimMax, NbPhase, nbCell))
    allocate( divDensiteMassiqueFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
    allocate( divDensiteMassiqueNode(NbIncTotalPrimMax, NbPhase, nbNode))

    allocate( SmDensiteMassiqueCell(NbPhase, nbCell))
    allocate( SmDensiteMassiqueFrac(NbPhase, nbFrac))
    allocate( SmDensiteMassiqueNode(NbPhase, nbNode))

    ! pression
    allocate( divPressionCell(NbIncTotalPrimMax, nbCell) )
    allocate( divPressionFrac(NbIncTotalPrimMax, nbFrac) )
    allocate( divPressionNode(NbIncTotalPrimMax, nbNode) )

    allocate( SmPressionCell( nbCell))
    allocate( SmPressionFrac( nbFrac))
    allocate( SmPressionNode( nbNode))

    ! pression capillaire
    allocate( PressionCapCell(NbPhase, nbCell))
    allocate( PressionCapFrac(NbPhase, nbFrac))
    allocate( PressionCapNode(NbPhase, nbNode))

    allocate( divPressionCapCell(NbIncTotalPrimMax, NbPhase, nbCell))
    allocate( divPressionCapFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
    allocate( divPressionCapNode(NbIncTotalPrimMax, NbPhase, nbNode))

    ! Saturation
    allocate( divSaturationCell(NbIncTotalPrimMax, NbPhase, nbCell))
    allocate( divSaturationFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
    allocate( divSaturationNode(NbIncTotalPrimMax, NbPhase, nbNode))
#ifdef _WIP_FREEFLOW_STRUCTURES_
    ! Freeflow phase molar flowrates
    allocate( divFreeFlowMolarFlowrateNode(NbIncTotalPrimMax, NbPhase, nbNode))
    allocate( SmFreeFlowMolarFlowrateNode(NbPhase, nbNode))
#endif
    ! DensiteMolaire*Kr/Viscosite * Comp
    allocate( DensiteMolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
    allocate( DensiteMolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
    allocate( DensiteMolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

    allocate( divDensiteMolaireKrViscoCompCell(NbIncTotalPrimMax, NbComp, NbPhase, nbCell))
    allocate( divDensiteMolaireKrViscoCompFrac(NbIncTotalPrimMax, NbComp, NbPhase, nbFrac))
    allocate( divDensiteMolaireKrViscoCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))

    allocate( SmDensiteMolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
    allocate( SmDensiteMolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
    allocate( SmDensiteMolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

    ! DensiteMolaire * Saturation * Comp
    allocate( DensiteMolaireSatCompCell(NbComp, NbPhase, nbCell))
    allocate( DensiteMolaireSatCompFrac(NbComp, NbPhase, nbFrac))
    allocate( DensiteMolaireSatCompNode(NbComp, NbPhase, nbNode))

    allocate( divDensiteMolaireSatCompCell(NbIncTotalPrimMax, NbComp, NbPhase, nbCell))
    allocate( divDensiteMolaireSatCompFrac(NbIncTotalPrimMax, NbComp, NbPhase, nbFrac))
    allocate( divDensiteMolaireSatCompNode(NbIncTotalPrimMax, NbComp, NbPhase, nbNode))

    allocate( SmDensiteMolaireSatCompCell(NbComp, NbPhase, nbCell))
    allocate( SmDensiteMolaireSatCompFrac(NbComp, NbPhase, nbFrac))
    allocate( SmDensiteMolaireSatCompNode(NbComp, NbPhase, nbNode))

    ! well inj
    allocate(DensiteMolaireKrViscoCompWellInj(NbComp, nbNodeInj))
    allocate(divDensiteMolaireKrViscoCompWellInj(NbComp, nbNodeInj))

    allocate(DensiteMolaireKrViscoEnthalpieWellInj(nbNodeInj))
    allocate(divDensiteMolaireKrViscoEnthalpieWellInj(nbNodeInj))

    ! the following arrays must be allocated even if there is no energy transfer
    allocate( SmTemperatureCell( nbCell))
    allocate( SmTemperatureFrac( nbFrac))
    allocate( SmTemperatureNode( nbNode))

#ifdef _THERMIQUE_

    ! temperature
    allocate( divTemperatureCell(NbIncTotalPrimMax, nbCell) )
    allocate( divTemperatureFrac(NbIncTotalPrimMax, nbFrac) )
    allocate( divTemperatureNode(NbIncTotalPrimMax, nbNode) )

    ! DensiteMolaire * PermRel / Viscosite * Enthalpie
    allocate( DensiteMolaireKrViscoEnthalpieCell(NbPhase, nbCell))
    allocate( DensiteMolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
    allocate( DensiteMolaireKrViscoEnthalpieNode(NbPhase, nbNode))

    allocate( divDensiteMolaireKrViscoEnthalpieCell(NbIncTotalPrimMax, NbPhase, nbCell))
    allocate( divDensiteMolaireKrViscoEnthalpieFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
    allocate( divDensiteMolaireKrViscoEnthalpieNode(NbIncTotalPrimMax, NbPhase, nbNode))

    allocate( SmDensiteMolaireKrViscoEnthalpieCell(NbPhase, nbCell))
    allocate( SmDensiteMolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
    allocate( SmDensiteMolaireKrViscoEnthalpieNode(NbPhase, nbNode))

    ! densitemolaire * energieinterne * Saturation
    allocate( DensiteMolaireEnergieInterneSatCell(NbPhase, nbCell))
    allocate( DensiteMolaireEnergieInterneSatFrac(NbPhase, nbFrac))
    allocate( DensiteMolaireEnergieInterneSatNode(NbPhase, nbNode))

    allocate( divDensiteMolaireEnergieInterneSatCell(NbIncTotalPrimMax, NbPhase, nbCell))
    allocate( divDensiteMolaireEnergieInterneSatFrac(NbIncTotalPrimMax, NbPhase, nbFrac))
    allocate( divDensiteMolaireEnergieInterneSatNode(NbIncTotalPrimMax, NbPhase, nbNode))

    allocate( SmDensiteMolaireEnergieInterneSatCell(NbPhase, nbCell))
    allocate( SmDensiteMolaireEnergieInterneSatFrac(NbPhase, nbFrac))
    allocate( SmDensiteMolaireEnergieInterneSatNode(NbPhase, nbNode))

#endif

  end subroutine LoisThermoHydro_allocate


  ! free
  subroutine LoisThermoHydro_free

   ! densite massique
    deallocate( DensiteMassiqueCell)
    deallocate( DensiteMassiqueFrac)
    deallocate( DensiteMassiqueNode)
    deallocate( divDensiteMassiqueCell)
    deallocate( divDensiteMassiqueFrac)
    deallocate( divDensiteMassiqueNode)
    deallocate( SmDensiteMassiqueCell)
    deallocate( SmDensiteMassiqueFrac)
    deallocate( SmDensiteMassiqueNode)

    ! pression
    deallocate( divPressionCell)
    deallocate( divPressionFrac)
    deallocate( divPressionNode)
    deallocate( SmPressionCell)
    deallocate( SmPressionFrac)
    deallocate( SmPressionNode)

    ! pression capillaire
    deallocate( PressionCapCell)
    deallocate( PressionCapFrac)
    deallocate( PressionCapNode)
    deallocate( divPressionCapCell)
    deallocate( divPressionCapFrac)
    deallocate( divPressionCapNode)

    ! saturation
    deallocate( divSaturationCell)
    deallocate( divSaturationFrac)
    deallocate( divSaturationNode)

#ifdef _WIP_FREEFLOW_STRUCTURES_
    ! Freeflow phase molar flowrates
    deallocate( divFreeFlowMolarFlowrateNode)
    deallocate( SmFreeFlowMolarFlowrateNode)
#endif

    ! densitemolaire * Permrel / viscosite * Comp
    deallocate( DensiteMolaireKrViscoCompCell)
    deallocate( DensiteMolaireKrViscoCompFrac)
    deallocate( DensiteMolaireKrViscoCompNode)
    deallocate( divDensiteMolaireKrViscoCompCell)
    deallocate( divDensiteMolaireKrViscoCompFrac)
    deallocate( divDensiteMolaireKrViscoCompNode)
    deallocate( SmDensiteMolaireKrViscoCompCell)
    deallocate( SmDensiteMolaireKrViscoCompFrac)
    deallocate( SmDensiteMolaireKrViscoCompNode)

    ! densitemolaire * Sat * Comp
    deallocate( DensiteMolaireSatCompCell)
    deallocate( DensiteMolaireSatCompFrac)
    deallocate( DensiteMolaireSatCompNode)
    deallocate( divDensiteMolaireSatCompCell)
    deallocate( divDensiteMolaireSatCompFrac)
    deallocate( divDensiteMolaireSatCompNode)
    deallocate( SmDensiteMolaireSatCompCell)
    deallocate( SmDensiteMolaireSatCompFrac)
    deallocate( SmDensiteMolaireSatCompNode)

    ! well inj
    deallocate(DensiteMolaireKrViscoCompWellInj)
    deallocate(divDensiteMolaireKrViscoCompWellInj)

    deallocate(DensiteMolaireKrViscoEnthalpieWellInj)
    deallocate(divDensiteMolaireKrViscoEnthalpieWellInj)

    deallocate( SmTemperatureCell)
    deallocate( SmTemperatureFrac)
    deallocate( SmTemperatureNode)

#ifdef _THERMIQUE_
    ! temperature
    deallocate( divTemperatureCell)
    deallocate( divTemperatureFrac)
    deallocate( divTemperatureNode)

    ! densitemolaire * Permrel / viscosite * Enthalpie
    deallocate( DensiteMolaireKrViscoEnthalpieCell)
    deallocate( DensiteMolaireKrViscoEnthalpieFrac)
    deallocate( DensiteMolaireKrViscoEnthalpieNode)
    deallocate( divDensiteMolaireKrViscoEnthalpieCell)
    deallocate( divDensiteMolaireKrViscoEnthalpieFrac)
    deallocate( divDensiteMolaireKrViscoEnthalpieNode)
    deallocate( SmDensiteMolaireKrViscoEnthalpieCell)
    deallocate( SmDensiteMolaireKrViscoEnthalpieFrac)
    deallocate( SmDensiteMolaireKrViscoEnthalpieNode)

    ! densitemolaire * energieinterne * Saturation
    deallocate( DensiteMolaireEnergieInterneSatCell)
    deallocate( DensiteMolaireEnergieInterneSatFrac)
    deallocate( DensiteMolaireEnergieInterneSatNode)
    deallocate( divDensiteMolaireEnergieInterneSatCell)
    deallocate( divDensiteMolaireEnergieInterneSatFrac)
    deallocate( divDensiteMolaireEnergieInterneSatNode)
    deallocate( SmDensiteMolaireEnergieInterneSatCell)
    deallocate( SmDensiteMolaireEnergieInterneSatFrac)
    deallocate( SmDensiteMolaireEnergieInterneSatNode)

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
