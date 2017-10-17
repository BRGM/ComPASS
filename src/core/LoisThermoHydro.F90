!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LoisThermoHydro

  use DefModel
  use NumbyContext
  use IncCV

  implicit none

  ! dXssurdXp and SmdXs
  double precision, allocatable, dimension(:,:,:), protected :: &
       dXssurdXpCell, &
       dXssurdXpFrac, &
       dXssurdXpNode

  double precision, allocatable, dimension(:,:), protected :: &
       SmdXsCell, &
       SmdXsFrac, &
       SmdXsNode

  double precision, allocatable, dimension(:,:), protected :: &
       SmFCell, &
       SmFFrac, &
       SmFNode

  ! num inc prim secd
  integer, allocatable, dimension(:,:), protected :: &
       NumIncPTCSPrimCell,  &
       NumIncPTCSPrimFrac,  &
       NumIncPTCSPrimNode,  &
       NumIncPTCSecondCell, &
       NumIncPTCSecondFrac, &
       NumIncPTCSecondNode

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

  ! DensiteMolaire*Kr/Viscosite*Comp
  double precision, allocatable, dimension(:,:,:), protected :: &
       DensitemolaireKrViscoCompCell, &
       DensitemolaireKrViscoCompFrac, &
       DensitemolaireKrViscoCompNode
  double precision, allocatable, dimension(:,:,:,:), protected :: &
       divDensitemolaireKrViscoCompCell, &
       divDensitemolaireKrViscoCompFrac, &
       divDensitemolaireKrViscoCompNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       SmDensitemolaireKrViscoCompCell, &
       SmDensitemolaireKrViscoCompFrac, &
       SmDensitemolaireKrViscoCompNode

  ! DensiteMolaire*Kr/Viscosite*Comp for wells (injection and production)
  double precision, allocatable, dimension(:,:), protected :: &
       DensitemolaireKrViscoCompWellInj
  double precision, allocatable, dimension(:,:), protected :: &
       divDensitemolaireKrViscoCompWellInj

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
       DensitemolaireKrViscoEnthalpieCell, &
       DensitemolaireKrViscoEnthalpieFrac, &
       DensitemolaireKrViscoEnthalpieNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divDensitemolaireKrViscoEnthalpieCell, &
       divDensitemolaireKrViscoEnthalpieFrac, &
       divDensitemolaireKrViscoEnthalpieNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmDensitemolaireKrViscoEnthalpieCell, &
       SmDensitemolaireKrViscoEnthalpieFrac, &
       SmDensitemolaireKrViscoEnthalpieNode

  ! DensiteMolaire * PermRel / Viscosite * Enthalpie for injection wells
  double precision, allocatable, dimension(:), protected :: &
       DensitemolaireKrViscoEnthalpieWellInj
  double precision, allocatable, dimension(:), protected :: &
       divDensitemolaireKrViscoEnthalpieWellInj


  ! densitemolaire * energieinterne * Saturation
  double precision, allocatable, dimension(:,:), protected :: &
       DensitemolaireEnergieInterneSatCell, &
       DensitemolaireEnergieInterneSatFrac, &
       DensitemolaireEnergieInterneSatNode
  double precision, allocatable, dimension(:,:,:), protected :: &
       divDensitemolaireEnergieInterneSatCell, &
       divDensitemolaireEnergieInterneSatFrac, &
       divDensitemolaireEnergieInterneSatNode
  double precision, allocatable, dimension(:,:), protected :: &
       SmDensitemolaireEnergieInterneSatCell, &
       SmDensitemolaireEnergieInterneSatFrac, &
       SmDensitemolaireEnergieInterneSatNode


  ! tmp values to simpfy notations of numerotation
  ! ex. NbPhasePresente = NbPhasePresente_ctx(inc%ic)
  integer, private :: &
       NbPhasePresente, NbCompCtilde, &
       NbEqFermeture, NbEqEquilibre,  &
       NbIncPTC, NbIncPTCPrim, NbIncPTCS, &
       NbIncPTCSPrim, NbIncPTCSecond, &
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
       LoisThermoHydro_init_cv,                & ! init infos according to ic (context) for each control volume (cv)
       LoisThermoHydro_dFsurdX_cv,             & ! compute dF/dX for each cv
       LoisThermoHydro_ps_cv,                  & ! choose first and second variables for each cv
       LoisThermoHydro_dXssurdXp_cv,           & ! compute dFs/dXp
                                !
       LoisThermoHydro_Densitemolaire_cv,      & ! prim divs: densitemolaire
       LoisThermoHydro_Viscosite_cv,           & !            1/viscosite
       LoisThermoHydro_PermRel_cv,             & !            Permrel
       LoisThermoHydro_Pression_cv,            & !            Pression
       LoisThermoHydro_PressionCapillaire_cv,  & !            Pressioncapillaire
       LoisThermoHydro_Saturation_cv,          & !            Saturation
                                !
       LoisThermoHydro_Densitemassique_cv,       & !          densitemassique
       LoisThermoHydro_DensitemolaireKrViscoComp_cv,     & !  densitemolaire * Permrel / viscosite * Comp
       LoisThermoHydro_DensitemolaireSatComp_cv,         & !  densitemolaire * Saturation * Comp
       LoisThermoHydro_DensitemolaireKrViscoEnthalpie_cv, & !  densitemolaire * Permrel / viscosite * Enthalpie
       LoisThermoHydro_DensitemolaireEnergieInterneSat_cv  !  densitemolaire * energieinterne * Saturation

#ifdef _THERMIQUE_

  private :: &
       LoisThermoHydro_Temperature_cv,     &
       LoisThermoHydro_EnergieInterne_cv,  & !  Enthalpie
       LoisThermoHydro_Enthalpie_cv
#endif

contains

  ! main subroutine of this module
  ! compute all prim div for future use
  ! loop of cell/frac/node
  ! subroutine Loisthermohydro_divPrim_cv compute all prim div for one control volume
  subroutine LoisThermoHydro_compute

    integer :: k
    integer :: rocktypeinc

    ! cell
    do k=1, NbCellLocal_Ncpus(commRank+1)

       rocktypeinc = 1

       call LoisThermoHydro_divPrim_cv(rocktypeinc, IncCell(k), &
                                !
            dXssurdXpCell(:,:,k), &
            SmdXsCell(:,k), &
            SmFCell(:,k),   &
                                !
            NumIncPTCSPrimCell(:,k),  &
            NumIncPTCSecondCell(:,k), &
                                !
            DensitemassiqueCell(:,k),       &
            divDensitemassiqueCell(:,:,k),  &
            SmDensitemassiqueCell(:,k),     &
                                !
            divPressionCell(:,k),         &
            SmPressionCell(k),            &
                                !
            divTemperatureCell(:,k),         &
            SmTemperatureCell(k),            &
                                !
            divSaturationCell(:,:,k),       &
                                !
            PressionCapCell(:,k),           &
            divPressionCapCell(:,:,k),      &
                                !
            DensitemolaireSatCompCell(:,:,k),      &
            divDensitemolaireSatCompCell(:,:,:,k), &
            SmDensitemolaireSatCompCell(:,:,k),    &
                                !
            DensitemolaireKrViscoCompCell(:,:,k),      &
            divDensitemolaireKrViscoCompCell(:,:,:,k), &
            SmDensitemolaireKrViscoCompCell(:,:,k),    &
                                !
            DensitemolaireEnergieInterneSatCell(:,k),      &
            divDensitemolaireEnergieInterneSatCell(:,:,k), &
            SmDensitemolaireEnergieInterneSatCell(:,k),    &
                                !
            DensitemolaireKrViscoEnthalpieCell(:,k),      &
            divDensitemolaireKrViscoEnthalpieCell(:,:,k), &
            SmDensitemolaireKrViscoEnthalpieCell(:,k) )
    end do

    ! frac
    do k=1, NbFracLocal_Ncpus(commRank+1)

       rocktypeinc = 2

       call LoisThermoHydro_divPrim_cv(rocktypeinc, IncFrac(k), &
                                !
            dXssurdXpFrac(:,:,k), &
            SmdXsFrac(:,k), &
            SmFFrac(:,k),   &
                                !
            NumIncPTCSPrimFrac(:,k),  &
            NumIncPTCSecondFrac(:,k), &
                                !
            DensitemassiqueFrac(:,k),       &
            divDensitemassiqueFrac(:,:,k),  &
            SmDensitemassiqueFrac(:,k),     &
                                !
            divPressionFrac(:,k),         &
            SmPressionFrac(k),            &
                                !
            divTemperatureFrac(:,k),         &
            SmTemperatureFrac(k),            &
                                !
            divSaturationFrac(:,:,k),       &
                                !
            PressionCapFrac(:,k),           &
            divPressionCapFrac(:,:,k),      &
                                !
            DensitemolaireSatCompFrac(:,:,k),      &
            divDensitemolaireSatCompFrac(:,:,:,k), &
            SmDensitemolaireSatCompFrac(:,:,k),    &
                                !
            DensitemolaireKrViscoCompFrac(:,:,k),      &
            divDensitemolaireKrViscoCompFrac(:,:,:,k), &
            SmDensitemolaireKrViscoCompFrac(:,:,k),    &
                                !
            DensitemolaireEnergieInterneSatFrac(:,k),      &
            divDensitemolaireEnergieInterneSatFrac(:,:,k), &
            SmDensitemolaireEnergieInterneSatFrac(:,k),    &
                                !
            DensitemolaireKrViscoEnthalpieFrac(:,k),      &
            divDensitemolaireKrViscoEnthalpieFrac(:,:,k), &
            SmDensitemolaireKrViscoEnthalpieFrac(:,k) )
    end do

    ! node
    do k=1, NbNodeLocal_Ncpus(commRank+1)

       if(IdNodeLocal(k)%Frac=="y") then
          rocktypeinc = 2
       else if(IdNodeLocal(k)%Frac=="n") then
          rocktypeinc = 1
       end if

       call LoisThermoHydro_divPrim_cv(rocktypeinc, IncNode(k), &
                                !
            dXssurdXpNode(:,:,k), &
            SmdXsNode(:,k), &
            SmFNode(:,k),   &
                                !
            NumIncPTCSPrimNode(:,k),  &
            NumIncPTCSecondNode(:,k), &
                                !
            DensitemassiqueNode(:,k),       &
            divDensitemassiqueNode(:,:,k),  &
            SmDensitemassiqueNode(:,k),     &
                                !
            divPressionNode(:,k),         &
            SmPressionNode(k),            &
                                !
            divTemperatureNode(:,k),         &
            SmTemperatureNode(k),            &
                                !
            divSaturationNode(:,:,k),       &
                                !
            PressionCapNode(:,k),           &
            divPressionCapNode(:,:,k),      &
                                !
            DensitemolaireSatCompNode(:,:,k),      &
            divDensitemolaireSatCompNode(:,:,:,k), &
            SmDensitemolaireSatCompNode(:,:,k),    &
                                !
            DensitemolaireKrViscoCompNode(:,:,k),      &
            divDensitemolaireKrViscoCompNode(:,:,:,k), &
            SmDensitemolaireKrViscoCompNode(:,:,k),    &
                                !
            DensitemolaireEnergieInterneSatNode(:,k),      &
            divDensitemolaireEnergieInterneSatNode(:,:,k), &
            SmDensitemolaireEnergieInterneSatNode(:,k),    &
                                !
            DensitemolaireKrViscoEnthalpieNode(:,k),      &
            divDensitemolaireKrViscoEnthalpieNode(:,:,k), &
            SmDensitemolaireKrViscoEnthalpieNode(:,k) )
    end do


    ! well injection
    !   compute q_{w,s,i} (and directive of pression)
    !   for well k and all nodes of well k
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       call LoisThermoHydro_divP_wellinj(k)
    end do

    ! SmDensitemolaireKrViscoEnthalpieNode(:,:) = 0.d0
    ! SmDensitemolaireEnergieInterneSatNode(:,:) = 0.d0
    ! SmDensitemolaireKrViscoCompNode(:,:,:) = 0.d0
    ! SmDensitemolaireSatCompNode(:,:,:) = 0.d0

    ! SmDensitemolaireKrViscoEnthalpieFrac(:,:) = 0.d0
    ! SmDensitemolaireEnergieInterneSatFrac(:,:) = 0.d0
    ! SmDensitemolaireKrViscoCompFrac(:,:,:) = 0.d0
    ! SmDensitemolaireSatCompFrac(:,:,:) = 0.d0

    ! SmDensitemolaireKrViscoEnthalpieCell(:,:) = 0.d0
    ! SmDensitemolaireEnergieInterneSatCell(:,:) = 0.d0
    ! SmDensitemolaireKrViscoCompCell(:,:,:) = 0.d0
    ! SmDensitemolaireSatCompCell(:,:,:) = 0.d0

    ! SmDensitemassiqueCell(:,:) = 0.d0
    ! SmDensitemassiqueNode(:,:) = 0.d0
    ! SmDensitemassiqueFrac(:,:) = 0.d0

  end subroutine LoisThermoHydro_compute


  ! all operations for one cv
  subroutine LoisThermoHydro_divPrim_cv(rocktypeinc, inc, &
       dXssurdXp, SmdXs, SmF, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       Densitemassique, divDensitemassique, SmDensitemassique, &
       divPression, SmPression, &
       divTemperature, SmTemperature, &
       divSaturation, &
       PressionCap, divPressionCap, &
       DensitemolaireSatComp, divDensitemolaireSatComp, SmDensitemolaireSatComp, &
       DensitemolaireKrViscoComp, divDensitemolaireKrViscoComp, SmDensitemolaireKrViscoComp, &
       DensitemolaireEnergieInterneSat, divDensitemolaireEnergieInterneSat, SmDensitemolaireEnergieInterneSat, &
       DensitemolaireKrViscoEnthalpie,  divDensitemolaireKrViscoEnthalpie,  SmDensitemolaireKrViscoEnthalpie)

    ! input
    type(Type_IncCV), intent(in) :: inc

    integer, intent(in) :: rocktypeinc

    ! output

    integer, intent(out) ::  &
         NumIncPTCSPrimCV (NbIncPTCSPrimMax),  &
         NumIncPTCSecondCV (NbEqFermetureMax)

    double precision, intent(out) :: &
                                !
         dXssurdXp (NbIncPTCSPrimMax, NbEqFermetureMax), & ! (col,row) index order
         SmdXs (NbEqFermetureMax), &
         SmF (NbEqFermetureMax),   &
                                !
         DensiteMassique (NbPhase), &
         divDensiteMassique (NbIncPTCSPrimMax, NbPhase), &
         SmDensiteMassique (NbPhase), &
                                !
         divPression( NbIncPTCSPrimMax), &
         SmPression, &
                                !
         divTemperature( NbIncPTCSPrimMax), &
         SmTemperature, &
                                !
         divSaturation ( NbIncPTCSPrimMax, NbPhase), &
                                !
         PressionCap (NbPhase), &
         divPressionCap (NbIncPTCSPrimMax, NbPhase),&
                                !
         DensiteMolaireSatComp (NbComp, NbPhase), &
         divDensiteMolaireSatComp (NbIncPTCSPrimMax, NbComp, NbPhase), &
         SmDensiteMolaireSatComp (NbComp, NbPhase), &
                                !
         DensitemolaireKrViscoComp (NbComp, NbPhase), &
         divDensitemolaireKrViscoComp (NbIncPTCSPrimMax, NbComp, NbPhase), &
         SmDensitemolaireKrViscoComp (NbComp, NbPhase), &
                                !
         DensitemolaireEnergieInterneSat (NbPhase), &
         divDensitemolaireEnergieInterneSat (NbIncPTCSPrimMax, NbPhase), &
         SmDensitemolaireEnergieInterneSat (NbPhase), &
                                !
         DensitemolaireKrViscoEnthalpie (NbPhase), &
         divDensitemolaireKrViscoEnthalpie (NbIncPTCSPrimMax, NbPhase), &
         SmDensitemolaireKrViscoEnthalpie (NbPhase)

    ! tmp
    double precision :: &
         dFsurdX( NbIncPTCSMax, NbEqFermetureMax) ! (col,row) index order

    double precision :: &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase), &
                                !
         UnsurViscosite(NbPhase), &
         divUnsurViscosite(NbIncPTCSPrimMax, NbPhase), &
         SmUnsurViscosite(NbPhase), &
                                !
         PermRel(NbPhase), &
         divPermRel(NbIncPTCSPrimMax, NbPhase)


    double precision :: &
         EnergieInterne(NbPhase), &
         divEnergieInterne(NbIncPTCSPrimMax, NbPhase), &
         SmEnergieInterne(NbPhase), &
                                !
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncPTCSPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)


    ! init tmp values for each cv
    call LoisThermoHydro_init_cv(inc)

    ! compute dF/dX
    ! dFsurdX: (col, row) index order
    call LoisThermoHydro_dFsurdX_cv(inc, dFsurdX, SmF)

    ! choose inconnues prim and secd
    call LoisThermoHydro_ps_cv(inc, dFsurdX, pschoice, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! compute dXssurdxp
    call LoisThermoHydro_dXssurdXp_cv(inc, dFsurdX, SmF, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, dXssurdXp, SmdXs)

    ! viscosite
    call LoisThermoHydro_viscosite_cv(inc, dXssurdXp, SmdXs, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite)

    ! densite massique
    call LoisThermoHydro_densitemassique_cv(inc, dXssurdXp, SmdXs, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         DensiteMassique, divDensiteMassique, SmDensiteMassique)

    ! deniste molaire
    call LoisThermoHydro_densitemolaire_cv(inc, dXssurdXp, SmdXs, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire)

    ! PermRel
    call LoisThermoHydro_PermRel_cv(inc, PermRel, divPermRel)

    ! Pression
    call LoisThermoHydro_Pression_cv(inc, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         dXssurdXp, SmdXs,  &
         divPression, SmPression)

    ! Pression Capillaire
    call LoisThermoHydro_PressionCapillaire_cv(rocktypeinc, inc, PressionCap, divPressionCap)

    ! Saturation div
    call LoisThermoHydro_Saturation_cv(inc, divSaturation)

    ! term: DensiteMolaire * PermRel / Viscosite * Comp
    call LoisThermoHydro_DensitemolaireKrViscoComp_cv( inc, &
         Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
         PermRel, divPermRel, &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV,  &
         dXssurdXp, SmdXs, &
         DensitemolaireKrViscoComp, divDensitemolaireKrViscoComp, SmDensitemolaireKrViscoComp)

    ! term: DensiteMolaire * Saturation * Comp
    call LoisThermoHydro_DensitemolaireSatComp_cv( &
         inc, divSaturation, &
         Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV,  &
         dXssurdXp, SmdXs, &
         DensitemolaireSatComp, divDensiteMolaireSatComp, SmDensiteMolaireSatComp)


#ifdef _THERMIQUE_

    ! Temperature
    call LoisThermoHydro_Temperature_cv(inc, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         dXssurdXp, SmdXs,  &
         divTemperature, SmTemperature)

    ! energie interne
    call LoisThermoHydro_EnergieInterne_cv(inc, dXssurdXp, SmdXs, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         EnergieInterne, divEnergieInterne, SmEnergieInterne)

    ! Enthalpie
    call LoisThermoHydro_Enthalpie_cv(inc, dXssurdXp, SmdXs, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         Enthalpie, divEnthalpie, SmEnthalpie)

    ! term: DensiteMolaire * Energieinterne * Saturation
    call LoisThermoHydro_DensitemolaireEnergieInterneSat_cv(  &
         inc, divSaturation, &
         DensiteMolaire, divDensiteMolaire, SmDensiteMolaire, &
         EnergieInterne, divEnergieInterne, SmEnergieInterne, &
         DensitemolaireEnergieInterneSat, divDensiteMolaireEnergieInterneSat, SmDensitemolaireEnergieInterneSat)

    ! term: DensiteMolaire * PermRel / Viscosite * Enthalpie
    call LoisThermoHydro_DensitemolaireKrViscoEnthalpie_cv( &
         Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
         PermRel, divPermRel, &
         UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
         Enthalpie, divEnthalpie, SmEnthalpie, &
         DensitemolaireKrViscoEnthalpie, divDensitemolaireKrViscoEnthalpie, SmDensitemolaireKrViscoEnthalpie)
#endif

    ! if(commRank==1) then
    !    do i=1, NbPhase

    ! print*, Densitemolaire(i)
    ! print*, ""

    ! do j=1, NbIncPTCSPrim
    !    print*, i, j, divDensitemolaire(j,i)
    ! end do
    ! print*, ""
    ! print*, SmDensitemolaire(i)
    ! print*, ""

    ! print*, 1.d0/UnsurViscosite(i)

    ! print*, Enthalpie(i)
    ! print*, divEnthalpie(:,i)
    ! print*, ""

    ! print*, PermRel(i)
    ! do j=1, NbIncPTCSPrim
    !    print*, i, j, divPermRel(j,i)
    ! end do
    ! print*, ""
    ! print*, SmDensitemolaire(i)
    ! print*, ""

    ! print*, DensitemolaireKrViscoComp(i,1)
    ! print*, ""

    ! do j=1, NbIncPTCSPrim
    !    print*, divDensitemolaireKrViscoComp(j,i,1)
    !    ! print*, divDensitemolaire
    ! end do
    ! print*, ""
    ! print*, SmDensitemolaireKrViscoComp(:,i)
    ! print*, ""

    ! do j=1, NbIncPTCSPrim
    !    write(*,'(ES22.14)'), divDensitemolaireKrViscoEnthalpie(j,i)
    !    ! print*, divDensitemolaire
    ! end do
    ! print*, ""
    ! write(*,*), SmDensitemolaireKrViscoEnthalpie(i)
    ! print*, ""

    ! do j=1, NbIncPTCSPrim
    !    write(*,'(ES22.14)'), divDensitemolaireEnergieInterneSat(j,i)
    ! end do
    ! print*, ""
    ! write(*,'(ES22.14)'), SmDensitemolaireEnergieInterneSat(i)
    ! print*, ""

    ! do j=1, NbIncPTCSPrim
    !    write(*,'(ES22.14)'), divSaturation(j,i)
    ! end do

    ! do j=1, NbIncPTCSPrim
    !    write(*,'(ES22.14)'), divTemperature(j)
    ! end do
    ! print*, ""
    ! print*, SmTemperature
    ! print*, ""

    !    end do
    ! end if

  end subroutine LoisThermoHydro_divPrim_cv

  !> Update thermo Laws of nodes
  subroutine LoisThermoHydro_divPrim_nodes

    integer :: k, rocktypeinc

    do k=1, NbNodeLocal_Ncpus(commRank+1)

       if(IdNodeLocal(k)%Frac=="y") then
          rocktypeinc = 2
       else if(IdNodeLocal(k)%Frac=="n") then
          rocktypeinc = 1
       end if

       call LoisThermoHydro_divPrim_cv(rocktypeinc, IncNode(k), &
                                !
            dXssurdXpNode(:,:,k), &
            SmdXsNode(:,k), &
            SmFNode(:,k),   &
                                !
            NumIncPTCSPrimNode(:,k),  &
            NumIncPTCSecondNode(:,k), &
                                !
            DensitemassiqueNode(:,k),       &
            divDensitemassiqueNode(:,:,k),  &
            SmDensitemassiqueNode(:,k),     &
                                !
            divPressionNode(:,k),         &
            SmPressionNode(k),            &
                                !
            divTemperatureNode(:,k),         &
            SmTemperatureNode(k),            &
                                !
            divSaturationNode(:,:,k),       &
                                !
            PressionCapNode(:,k),           &
            divPressionCapNode(:,:,k),      &
                                !
            DensitemolaireSatCompNode(:,:,k),      &
            divDensitemolaireSatCompNode(:,:,:,k), &
            SmDensitemolaireSatCompNode(:,:,k),    &
                                !
            DensitemolaireKrViscoCompNode(:,:,k),      &
            divDensitemolaireKrViscoCompNode(:,:,:,k), &
            SmDensitemolaireKrViscoCompNode(:,:,k),    &
                                !
            DensitemolaireEnergieInterneSatNode(:,k),      &
            divDensitemolaireEnergieInterneSatNode(:,:,k), &
            SmDensitemolaireEnergieInterneSatNode(:,k),    &
                                !
            DensitemolaireKrViscoEnthalpieNode(:,k),      &
            divDensitemolaireKrViscoEnthalpieNode(:,:,k), &
            SmDensitemolaireKrViscoEnthalpieNode(:,k) )
    end do

  end subroutine LoisThermoHydro_divPrim_nodes


  subroutine LoisThermoHydro_init_cv(inc)

    type(Type_IncCV), intent(in) :: inc

    NbPhasePresente = NbPhasePresente_ctx(inc%ic)
    NbCompCtilde = NbCompCtilde_ctx(inc%ic)

    NbEqFermeture = NbEqFermeture_ctx(inc%ic)
    NbEqEquilibre = NbEqEquilibre_ctx(inc%ic)

    NbIncPTC = NbIncPTC_ctx(inc%ic)
    NbIncPTCS = NbIncPTC + NbPhasePresente

    NbIncPTCPrim  = NbIncPTC - NbEqFermeture

    ! ps. if there is only one phase, phase is secd
    NbIncPTCSPrim = NbIncPTCS - NbEqFermeture - 1
    NbIncPTCSecond = NbEqFermeture

    NumPhasePresente(:) = NumPhasePresente_ctx(:,inc%ic)
    NumCompCtilde(:) = NumCompCtilde_ctx(:,inc%ic)
    NumCompEqEquilibre(:) = NumCompEqEquilibre_ctx(:,inc%ic)
    NumIncPTC2NumIncComp_comp(:) = NumIncPTC2NumIncComp_comp_ctx(:,inc%ic)
    NumIncPTC2NumIncComp_phase(:) = NumIncPTC2NumIncComp_phase_ctx(:,inc%ic)

    Num2PhasesEqEquilibre(:,:) = Num2PhasesEqEquilibre_ctx(:,:,inc%ic)
    NumIncComp2NumIncPTC(:,:) = NumIncComp2NumIncPTC_ctx(:,:,inc%ic)

  end subroutine LoisThermoHydro_init_cv


  ! compute dFsurdX for each control volume
  subroutine LoisThermoHydro_dFsurdX_cv(inc, dFsurdX, SmF)

    type(Type_IncCV), intent(in) :: inc
    double precision, intent(out) :: &  ! (col, row) index order
         dFsurdX(NbIncPTCSMax, NbEqFermetureMax)

    double precision, intent(out) :: &
         SmF(NbEqFermetureMax)

    integer :: i, mi, iph, iph1, iph2, icp, j, numj, numc1, numc2
    double precision :: &
         f1, dPf1, dTf1, dCf1(NbComp), &
         f2, dPf2, dTf2, dCf2(NbComp)

    dFsurdX(:,:) = 0.d0
    SmF(:) = 0.d0

    ! 1. sum_alpha C_i^alpha = 1, for i

    ! loop for rows associate with C_i^alpha in dFsurdX

    j = 1 + IndThermique

    do i=1, NbPhasePresente  ! row is i, col is j
       iph = NumPhasePresente(i)

       do icp=1, NbComp ! loop for cols
          if(MCP(icp,iph)==1) then
             j = j + 1
             dFsurdX(j,i) = 1.d0

             SmF(i) = SmF(i) + inc%Comp(icp,iph)
          end if
       end do
    enddo

    do i=1, NbPhasePresente
       SmF(i) = SmF(i) - 1.d0
    end do

    mi = NbPhasePresente ! row is mi+i

    ! 2. f_i^alpha * C_i^alpha = f_i^beta * C_i^beta
    do i=1, NbEqEquilibre !

       icp = NumCompEqEquilibre(i) ! component

       iph1 = Num2PhasesEqEquilibre(1,i) ! phase alpha
       iph2 = Num2PhasesEqEquilibre(2,i) ! phase beta

       numc1 = NumIncComp2NumIncPTC(icp,iph1) ! num of C_i^alpha in IncPTC
       numc2 = NumIncComp2NumIncPTC(icp,iph2) ! num of C_i^beta in IncPTC

       ! fugacity and div
       call f_Fugacity(iph1, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph1), inc%Saturation, &
            f1, dPf1, dTf1, dCf1)
       call f_Fugacity(iph2, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph2), inc%Saturation, &
            f2, dPf2, dTf2, dCf2)

       ! div Pression
       dFsurdX(1,i+mi) = dPf1*inc%Comp(icp,iph1) - dPf2*inc%Comp(icp,iph2)

#ifdef _THERMIQUE_

       ! div Temperature
       dFsurdX(2,i+mi) = dTf1*inc%Comp(icp,iph1) - dTf2*inc%Comp(icp,iph2)
#endif

       ! d (f(P,T,C)*C_i)/dC_i = f + df/dC_i*C_i
       ! d (f(P,T,C)*C_i)/dC_j = df/dC_i*C_j, j!=i
       dFsurdX(numc1, i+mi) = f1

       do j=1, NbComp
          if(MCP(j,iph1)==1) then ! phase iph1 contains component j

             numj = NumIncComp2NumIncPTC(j,iph1) ! num of C_j^iph1 in Inc
             dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
                  + inc%Comp(j,iph1)*dCf1(j)
          end if
       end do

       dFsurdX(numc2,i+mi) = -f2

       do j=1, NbComp
          if(MCP(j,iph2)==1) then ! phase iph1 contains component j

             numj = NumIncComp2NumIncPTC(j,iph2) ! num of C_j^iph2 in Inc
             dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
                  - inc%Comp(j,iph2)*dCf2(j)
          end if
       end do

       ! SmF
       SmF(i+mi) = f1*inc%Comp(icp,iph1) - f2*inc%Comp(icp,iph2)
    end do

  end subroutine LoisThermoHydro_dFsurdX_cv


  ! choose prim and sec inconnus for each CV
  ! fill inc%Nb/NumIncPrim/Secd
  subroutine LoisThermoHydro_ps_cv(inc, dFsurdX, pschoicecv, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    type(Type_IncCV), intent(in) :: inc
    double precision, intent(in) :: dFsurdX(NbIncPTCSMax, NbEqFermetureMax)
    integer, intent(in) :: pschoicecv
    integer, intent(out) :: NumIncPTCSPrimCV( NbIncPTCSPrimMax)
    integer, intent(out) :: NumIncPTCSecondCV( NbEqFermetureMax)

    integer :: i, ic

    NumIncPTCSPrimCV(:) = 0
    NumIncPTCSecondCV(:) = 0

    if(pschoicecv==1) then ! manually

       ic = inc%ic

       ! prim variable
       do i=1, NbIncPTCSPrim_ctx(ic)
          NumIncPTCSPrimCV(i) = psprim(i,ic)
       end do

       ! secd variable
       do i=1, NbEqFermeture_ctx(ic)
          NumIncPTCSecondCV(i) = pssecd(i,ic)
       end do

    else if(pschoicecv==2) then ! Glouton method
       ! call LoisThermoHydro_IncSecondGluton(inc, dFsurdX, &
       ! NumIncPTCSPrimCV, NumIncPTCSecondCV)

    else if (pschoicecv==3) then ! Gauss method
       call LoisThermoHydro_IncSecondGauss(inc, dFsurdX, &
            NumIncPTCSPrimCV, NumIncPTCSecondCV)
    end if

  end subroutine LoisThermoHydro_ps_cv


  ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
  subroutine LoisThermoHydro_dXssurdXp_cv(inc, dFsurdX, SmF, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, dXssurdXp, SmdXs)

    ! inputs
    type(Type_IncCV), intent(in) :: inc
    double precision, intent(in) ::  & ! (col, row) index order
         dFsurdX(NbIncPTCSMax, NbEqFermetureMax), &
         SmF(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! tmp
    double precision :: & ! (row,col) index order, lapack
         dFsurdX_prim(NbEqFermetureMax, NbIncPTCSPrimMax), &
         dFsurdX_secd(NbEqFermetureMax, NbEqFermetureMax)

    ! parameters for lapack
    integer :: ipiv(NbEqFermetureMax), info, Ierr, errcode
    integer :: i, j

    ! from dFsurdX, take out the cols of prim and secd variable
    do j=1, NbIncPTCSPrim
       do i=1, NbEqFermeture
          dFsurdX_prim(i,j) = dFsurdX(NumIncPTCSPrimCV(j),i)
       enddo
    enddo

    do j=1, NbEqFermeture ! = NbIncPTCSecond
       do i=1, NbEqFermeture
          dFsurdX_secd(i,j) = dFsurdX(NumIncPTCSecondCV(j),i)
       enddo
    enddo

    ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
    call dgetrf(NbEqFermeture, NbEqFermeture, &
         dFsurdX_secd, NbEqFermetureMax, ipiv, info)

    if(info /=0) then
       write(0,'(A,I0)') "dgetrf error in dXssurdxp, info = ", info
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    call dgetrs('N', NbEqFermeture, NbIncPTCSPrim, &
         dFsurdX_secd, NbEqFermetureMax, &
         ipiv, dFsurdX_prim, NbEqFermetureMax, info)
    if(info /=0) then
       write(0,'(A,I0)') "dgetrs error in dXssurdxp, info = ", info
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    do j=1, NbIncPTCSPrim
       do i=1, NbEqFermeture
          dXssurdXp(j,i) = dFsurdX_prim(i,j)
       enddo
    enddo

    ! SmdXs = dFsurdXs**(-1) * SmF
    call dgetrs('N', NbEqFermeture, 1, &
         dFsurdX_secd, NbEqFermetureMax, &
         ipiv, SmF, NbEqFermetureMax, info)
    if(info /=0) then
       write(0,'(A,I0)') "dgetrs error in SmdXs, info = ", info
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    do i=1, NbEqFermeture
       SmdXs(i) = SmF(i)
    end do

  end subroutine LoisThermoHydro_dXssurdXp_cv


  subroutine LoisThermoHydro_densitemassique_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    double precision, intent(in) :: &    ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)
    double precision, intent(out) :: Smval(NbPhase)

    ! tmp
    double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax, NbPhase) !=NbEqFermetureMax

    integer :: iph, i, j, jc

    ! 1. val
    ! 2. dval
    ! 3. Smval

    val(:) = 0.d0
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

#ifdef DEBUG_LOISTHEMOHYDRO
do iph = 1, NbPhase
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
       dfdX(1) = dPf  ! P

#ifdef _THERMIQUE_
       dfdX(2) = dTf
#endif

       do j= 2+IndThermique, NbIncPTC ! C
          dfdX(j) = 0.d0
       enddo
       do j=1, NbComp
          if(MCP(j,iph)==1) then
             jc = NumIncComp2NumIncPTC(j,iph)
             dfdX(jc) = dCf(j)
          end if
       enddo

       do j=1, NbPhase ! S
          jc = j + NbIncPTC
          dfdX(jc) = dSf(j)
       enddo

       ! prim and secd part of dfdX
       do j=1, NbIncPTCSPrim
          dval(j,iph) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
       end do

       do j=1, NbIncPTCSecond ! = NbEqFermeture
          dfdX_secd(j,iph) = dfdX( NumIncPTCSecondCV(j))
       end do

    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=densitemassique
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncPTCSPrim, NbPhase, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! - dv/dXs*SmdXs
    Smval(:) = 0.d0
    call dgemv('T', NbEqFermeture, NbPhase,  &
         -1.d0, dfdX_secd, NbIncPTCSecondMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_densitemassique_cv



  subroutine LoisThermoHydro_viscosite_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)
    double precision, intent(out) :: Smval(NbPhase)

    ! tmp
    double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax, NbPhase)

    double precision :: f2 ! =f**2
    integer :: i, iph, j, jc

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
            inc%Comp, inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = 1.d0/f ! val
       f2 = f**2

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       dfdX(1) = -dPf/f2  ! P

#ifdef _THERMIQUE_
       dfdX(2) = -dTf/f2
#endif

       do j= 2+IndThermique, NbIncPTC ! C
          dfdX(j) = 0.d0
       enddo
       do j=1, NbComp
          if(MCP(j,iph)==1) then
             jc = NumIncComp2NumIncPTC(j,iph)
             dfdX(jc) = -dCf(j)/f2
          end if
       enddo

       do j=1, NbPhase ! S
          jc = j + NbIncPTC
          dfdX(jc) = -dSf(j)/f2
       enddo

       ! prim and secd part of dfdX
       do j=1, NbIncPTCSPrim
          dval(j,i) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
       end do

       do j=1, NbIncPTCSecond ! = NbEqFermeture
          dfdX_secd(j,i) = dfdX( NumIncPTCSecondCV(j))
       end do
    end do

    !print*, SmdXs

    ! dv/dXp - dv/dXs*dXs/dXp, v=viscosite
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order,
    ! consider all the mats as transpose
    call dgemm('N','N', NbIncPTCSPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbIncPTCSecondMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_viscosite_cv



  subroutine LoisThermoHydro_densitemolaire_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)
    double precision, intent(out) :: Smval(NbPhase)

    ! tmp
    double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax, NbPhase)

    integer :: i, iph, j, jc

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
       dfdX(1) = dPf  ! P

#ifdef _THERMIQUE_
       dfdX(2) = dTf
#endif

       do j= 2+IndThermique, NbIncPTC ! C
          dfdX(j) = 0.d0
       enddo
       do j=1, NbComp
          if(MCP(j,iph)==1) then
             jc = NumIncComp2NumIncPTC(j,iph)
             dfdX(jc) = dCf(j)
          end if
       enddo

       do j=1, NbPhase ! S
          jc = j + NbIncPTC
          dfdX(jc) = dSf(j)
       enddo

       ! prim and secd part of dfdX
       do j=1, NbIncPTCSPrim
          dval(j,i) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
       end do

       do j=1, NbIncPTCSecond ! = NbEqFermeture
          dfdX_secd(j,i) = dfdX( NumIncPTCSecondCV(j))
       end do

    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=densitemolaire
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncPTCSPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbIncPTCSecondMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_densitemolaire_cv



  subroutine LoisThermoHydro_PermRel_cv(inc, val, dval)

    ! input
    type(Type_IncCV), intent(in)  :: inc

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)

    ! tmp
    double precision :: f, dSf(NbPhase), dfS_secd
    integer :: i, iph, j, jph

    val(:) = 0.d0
    dval(:,:) = 0.d0

    ! if there is only one presente phase,
    ! this Saturation mush be secd --> dval=0, Sm=0
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       call f_PermRel(iph, inc%Saturation, f, dSf)

       val(i) = f
       dfS_secd = dSf( NumPhasePresente(NbPhasePresente)) ! the last is secd

       ! alpha=1,2,...,NbPhasepresente-1
       ! Ps. NbIncPTCPrim+NbPhasePresente-1=NbIncPTCSPrim
       do j=1, NbPhasePresente - 1
          jph = NumPhasePresente(j)

          dval(j+NbIncPTCPrim,i) = dSf(jph) - dfS_secd
       end do
    end do

  end subroutine LoisThermoHydro_PermRel_cv


  subroutine LoisThermoHydro_Pression_cv( inc, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       dXssurdXp, SmdXs, &
       dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV( NbEqFermetureMax)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: dval(NbIncPTCSPrimMax)
    double precision, intent(out) :: Smval

    integer :: j

    dval(:) = 0.d0

    ! P is prim, first in NumIncPTCSPrimCV
    if( NumIncPTCSPrimCV(1)==1) then
       dval(1) = 1.d0
       Smval = 0.d0
    else
       do j=1, NbIncPTCSecond
          if( NumIncPTCSecondCV(j)==1) then

             ! row j of dXssurdXp (stored in (col, row) index order)
             dval(1:NbIncPTCSPrim) = -dXssurdXp(1:NbIncPTCSPrim,j)
             Smval = -SmdXs(j)
             exit
          end if

       end do
    end if

  end subroutine LoisThermoHydro_Pression_cv



  subroutine LoisThermoHydro_Saturation_cv(inc, dval)

    ! input
    type(Type_IncCV), intent(in)  :: inc

    ! output
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)

    ! tmp
    integer :: i

    dval(:,:) = 0.d0

    ! alpha=1,2,...,NbPhasePresente-1
    do i=1, NbPhasePresente - 1
       dval(i+NbIncPTCPrim,i) = 1.d0
    end do

    ! alpha = NbPhasePresente
    do i=1, NbPhasePresente - 1
       dval(i+NbIncPTCPrim, NbPhasePresente) = -1.d0
    end do

  end subroutine LoisThermoHydro_Saturation_cv



  subroutine LoisThermoHydro_PressionCapillaire_cv(rocktypeinc, inc, val, dval)

    ! input
    integer, intent(in) :: rocktypeinc
    type(Type_IncCV), intent(in)  :: inc

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)

    ! tmp
    double precision :: f, dSf(NbPhase), dfS_secd
    integer :: i, iph, j, jph

    val(:) = 0.d0
    dval(:,:) = 0.d0

    do iph=1, NbPhase

       call f_PressionCapillaire(rocktypeinc, iph, inc%Saturation, f, dSf)

       val(iph) = f

       dfS_secd = dSf( NumPhasePresente( NbPhasePresente))

       do j=1, NbPhasePresente - 1
          jph = NumPhasePresente(j)

          dval(j+NbIncPTCPrim,iph) = dSf(jph) - dfS_secd
       end do

    end do

  end subroutine LoisThermoHydro_PressionCapillaire_cv


#ifdef _THERMIQUE_

  subroutine LoisThermoHydro_EnergieInterne_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)
    double precision, intent(out) :: Smval(NbPhase)

    ! tmp
    double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax, NbPhase)

    integer :: i, iph, j, jc

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
            inc%Comp, inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       dfdX(1) = dPf  ! P

#ifdef _THERMIQUE_
       dfdX(2) = dTf
#endif

       do j= 2+IndThermique, NbIncPTC ! C
          dfdX(j) = 0.d0
       enddo
       do j=1, NbComp
          if(MCP(j,iph)==1) then
             jc = NumIncComp2NumIncPTC(j,iph)
             dfdX(jc) = dCf(j)
          end if
       enddo

       do j=1, NbPhase ! S
          jc = j + NbIncPTC
          dfdX(jc) = dSf(j)
       enddo

       ! prim and secd part of dfdX
       do j=1, NbIncPTCSPrim
          dval(j,i) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
       end do

       do j=1, NbIncPTCSecond ! = NbEqFermeture
          dfdX_secd(j,i) = dfdX( NumIncPTCSecondCV(j))
       end do
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=energieinterne
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncPTCSPrim, NbPhasePresente , NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbIncPTCSecondMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_EnergieInterne_cv

#endif


#ifdef _THERMIQUE_

  subroutine LoisThermoHydro_Temperature_cv( inc, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       dXssurdXp, SmdXs, &
       dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    integer, intent(in) :: &
         NumIncPTCSPrimCV( NbIncPTCSPrimMax), &
         NumIncPTCSecondCV( NbEqFermetureMax)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: dval(NbIncPTCSPrimMax)
    double precision, intent(out) :: Smval

    ! tmp
    double precision :: f, dPf
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax)

    integer :: j,jT

    dval(:) = 0.d0
    Smval = 0.d0

    ! P is prim, T is prim: NumIncPTCSPrimCV(1)=1, NumIncPTCSPrimCV(2)=2
    ! P is prim, T is secd: NumIncPTCSPrimCV(1)=1, NumIncPTCSPrimCV(2)>2
    ! TODO

    ! TODO
    ! Rewrite using ANY(NumIncPTCSPrimCV==1) and get index
    if((NumIncPTCSPrimCV(1)==1) .and. &
         (NumIncPTCSPrimCV(2)==2)) then ! P is prim, T is prim

       dval(2) = 1.d0
       Smval = 0.d0

    else if((NumIncPTCSPrimCV(1)==1) .and. &
         (NumIncPTCSPrimCV(2)>2)) then ! P is prim, T is secd

       ! T = Tsat(P)
    !   call DefModel_Tsat(inc%Pression, f, dPf)

       ! Fill dfdX = (df/dP, df/dT, df/dC, df/dS)
     !  dfdX(:) = 0.d0

     !  dfdX(1) = dPf
     !  dfdX(2) = 0.d0

       jT = 0
       do j=1,NbEqFermeture
          if (NumIncPTCSecondCV(j).eq.2) then
             jT = j
             dval(:) =  - dXssurdXp(:,j)
             Smval = - SmdXs(j)
          endif
       enddo
       if (jT.eq.0) then
          write(*,*)' pb dans derprim temp T second non trouvee ',jT
          stop
       endif

  ! autre calcul equivalent (mais plus couteux)
   !    dfdX(:) = 0.d0
   !    dfdX(2) = 1.d0
       ! prim and secd part of dfdX
    !   do j=1, NbIncPTCSPrim
    !      dval(j) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
    !   end do

    !   dfdX_secd(:) = 0.d0
    !   do j=1, NbIncPTCSecond ! = NbEqFermeture
    !      dfdX_secd(j) = dfdX( NumIncPTCSecondCV(j))
    !   end do

       ! dT/dXp - dT/dXs*dXs/dXp
       ! dval = dfdX_prim - dXssurdXp*dfdX_secd
       ! all the mats is in (col, row) index order, only need to consider as transpose
   !    call dgemm('N','N', NbIncPTCSPrim, 1, NbEqFermeture, &
   !         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
       !     dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! -dT/dXs*SmdXs
  !  call dgemv('T', NbEqFermeture, 1,  &
  !       -1.d0, dfdX_secd, NbIncPTCSecondMax, &
  !       SmdXs, 1, 0.d0, Smval, 1)

   ! debug: on compare les deux methodes : OK
   ! write(*,*)' - dXssurdXp 1 ',- dXssurdXp(:,jT)
   ! write(*,*)' - SmdXs 1 ',- SmdXs(jT)
   ! write(*,*)' dT ',dval
   ! write(*,*)' SmT ',Smval
   ! stop

    else if( NumIncPTCSPrimCV(1)==2) then

       ! TODO: P is secd, T is prim
       print*, "TODO error in div Temperature"

    else

       ! TODO: P and T are both secd ?
       print*, "TODO error in div Temperature"
    end if

  end subroutine LoisThermoHydro_Temperature_cv


  subroutine LoisThermoHydro_Enthalpie_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in)  :: inc
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)
    double precision, intent(out) :: Smval(NbPhase)

    ! tmp
    double precision :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision :: dfdX(NbIncPTCSMax)
    double precision :: dfdX_secd(NbIncPTCSecondMax, NbPhase)

    integer :: i, iph, j, jc

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
            inc%Comp, inc%Saturation, &
            f, dPf, dTf, dCf, dSf)

       val(i) = f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       dfdX(1) = dPf  ! P

#ifdef _THERMIQUE_
       dfdX(2) = dTf
#endif

       do j= 2+IndThermique, NbIncPTC ! C
          dfdX(j) = 0.d0
       enddo
       do j=1, NbComp
          if(MCP(j,iph)==1) then
             jc = NumIncComp2NumIncPTC(j,iph)
             dfdX(jc) = dCf(j)
          end if
       enddo

       do j=1, NbPhase ! S
          jc = j + NbIncPTC
          dfdX(jc) = dSf(j)
       enddo

       ! prim and secd part of dfdX
       do j=1, NbIncPTCSPrim
          dval(j,i) = dfdX( NumIncPTCSPrimCV(j)) ! dval=dfdX_prim
       end do

       do j=1, NbIncPTCSecond ! = NbEqFermeture
          dfdX_secd(j,i) = dfdX( NumIncPTCSecondCV(j))
       end do
    end do

    ! dv/dXp - dv/dXs*dXs/dXp, v=viscosite
    ! dval = dfdX_prim - dXssurdXp*dfdX_secd
    ! all the mats is in (col, row) index order, only need to consider as transpose
    call dgemm('N','N', NbIncPTCSPrim, NbPhasePresente, NbEqFermeture, &
         -1.d0, dXssurdXp, NbIncPTCSPrimMax, &
         dfdX_secd, NbEqFermetureMax, 1.d0, dval, NbIncPTCSPrimMax)

    ! -dv/dXs*SmdXs
    call dgemv('T', NbEqFermeture, NbPhasePresente,  &
         -1.d0, dfdX_secd, NbIncPTCSecondMax, &
         SmdXs, 1, 0.d0, Smval, 1)

  end subroutine LoisThermoHydro_Enthalpie_cv

#endif


  ! term: desitemolaire * Saturation
  subroutine LoisThermoHydro_DensitemolaireSat_cv( &
       inc, divSaturation, &
       Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
       val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         divSaturation(NbIncPTCSPrimMax, NbPhase), &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncPTCSPrimMax, NbPhase), Smval(NbPhase)

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

       do k=1, NbIncPTCSPrim

          dval(k,i) = &
               divDensitemolaire(k,i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       Smval(i) = SmDensitemolaire(i)*inc%Saturation(iph)
    end do

  end subroutine LoisThermoHydro_DensitemolaireSat_cv


  ! term: desitemolaire * Saturation * Comp
  subroutine LoisThermoHydro_DensitemolaireSatComp_cv( &
       inc, divSaturation, &
       Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV,  &
       dXssurdXp, SmdXs, &
       val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         divSaturation(NbIncPTCSPrimMax, NbPhase), &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncPTCSPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

    ! tmp
    integer :: i, iph, icp, k, j, jcp, jph, numj, s
    double precision :: dv

    double precision :: dvi(NbIncPTCSPrimMax)

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
       do k=1, NbIncPTCSPrim

          dvi(k) = &
               divDensitemolaire(k,i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)
       end do

       ! 2.1.2
       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             do k=1, NbIncPTCSPrim
                dval(k,icp,i) = dvi(k) * inc%Comp(icp,iph)
             end do

          end if
       end do ! end of P_i

    end do ! end of 2.1

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       dv = SmDensitemolaire(i)*inc%Saturation(iph)

       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             Smval(icp,i) = dv * inc%Comp(icp,iph)
          end if
       end do
    end do ! end of 2.2

    ! 3.1 + DensiteMolaire*Saturation * div(C_i^alpha)
    ! 3.2 + DensiteMolaire*Saturation * Sm(C_i^alpha)

    ! loop of inc prim
    do j=1, NbIncPTCSPrim
       numj = NumIncPTCSPrimCV(j)

       ! numi is C_i^alpha in (P,T,C,S)
       if( (numj>=(1+IndThermique)) .and. (numj<=NbIncPTC)) then

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

    ! loop of inc secd
    do j=1, NbIncPTCSecond ! = NbEqFermeture
       numj = NumIncPTCSecondCV(j)

       ! numi is C in (P,T,C)
       if( (numj>=(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then

                dv = DensiteMolaire(i)*inc%Saturation(iph)

                do s=1, NbIncPTCSPrim
                   dval(s,jcp,i) = dval(s,jcp,i) - dXssurdXp(s,j) * dv
                end do

                Smval(jcp,i) = Smval(jcp,i) - SmdXs(j) * dv
             end if
          end do

       end if
    end do ! end of inc secd

  end subroutine LoisThermoHydro_DensitemolaireSatComp_cv


  ! div and Sm of term: DensiteMolaire*Kr/Viscosite
  subroutine LoisThermoHydro_DensitemolaireKrVisco_cv( &
       Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel,             &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       val, dval, Smval)

    ! input
    double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase), &
                                !
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         divPermRel(NbIncPTCSPrimMax, NbPhase), &
         divUnsurViscosite(NbIncPTCSPrimMax, NbPhase), &
                                !
         SmDensiteMolaire(NbPhase),&
         SmUnSurViscosite(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncPTCSPrimMax, NbPhase), Smval(NbPhase)

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

       do k=1, NbIncPTCSPrim

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

  end subroutine LoisThermoHydro_DensitemolaireKrVisco_cv


  ! div and Sm of term: DensiteMolaire*Kr/Viscosite*Comp
  subroutine LoisThermoHydro_DensitemolaireKrViscoComp_cv(inc, &
       Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel,             &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV,  &
       dXssurdXp, SmdXs, &
       val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in) :: inc

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase),   &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase),   &
                                !
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         divPermRel(NbIncPTCSPrimMax, NbPhase), &
         divUnsurViscosite(NbIncPTCSPrimMax, NbPhase), &
                                !
         SmDensiteMolaire(NbPhase), &
         SmUnSurViscosite(NbPhase)

    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV(NbEqFermetureMax)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: &
         val(NbComp, NbPhase), &
         dval(NbIncPTCSPrimMax, NbComp, NbPhase), &
         Smval(NbComp, NbPhase)

    ! tmp
    integer :: i, iph, icp, k, j, jcp, jph, numj, s
    double precision :: dv

    double precision :: dvi(NbIncPTCSPrimMax)

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
             val(icp,i) = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i) &
                  * inc%Comp(icp,iph)
          end if
       end do
    end do

    ! 2.1 div(DensiteMolaire*Kr/Viscosite)*C_i^alpha
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       ! 2.1.1 compute dvi, tmp vector, used in 2.1.2
       do k=1, NbIncPTCSPrim
          dvi(k) = &
               divDensiteMolaire(k,i)*PermRel(i)*UnsurViscosite(i)    &
               + divPermRel(k,i)*DensiteMolaire(i)*UnsurViscosite(i)  &
               + divUnsurViscosite(k,i)*DensiteMolaire(i)*PermRel(i)
       end do

       ! 2.1.2
       do icp=1, NbComp ! P_i
          if(MCP(icp,iph)==1) then

             do k=1, NbIncPTCSPrim
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
    do j=1, NbIncPTCSPrim
       numj = NumIncPTCSPrimCV(j)

       ! numi is C_i^alpha in (P,T,C,S)
       if( (numj>=(1+IndThermique)) .and. (numj<=NbIncPTC)) then

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

    ! loop of inc secd
    do j=1, NbIncPTCSecond ! = NbEqFermeture
       numj = NumIncPTCSecondCV(j)

       ! numi is C in (P,T,C)
       if( (numj>=(1+IndThermique)) .and. (numj<=NbIncPTC)) then

          jcp = NumIncPTC2NumIncComp_comp(numj)
          jph = NumIncPTC2NumIncComp_phase(numj)

          ! find i that iph=jph
          do i=1, NbPhasePresente
             iph = NumPhasePresente(i)

             if(iph==jph) then

                dv = DensiteMolaire(i)*PermRel(i)*UnsurViscosite(i)

                do s=1, NbIncPTCSPrim
                   dval(s,jcp,i) = dval(s,jcp,i) - dXssurdXp(s,j) * dv
                end do

                Smval(jcp,i) = Smval(jcp,i) - SmdXs(j) * dv
             end if
          end do

       end if
    end do ! end of inc secd

  end subroutine LoisThermoHydro_DensitemolaireKrViscoComp_cv


  subroutine LoisThermoHydro_DensitemolaireEnergieInterneSat_cv( &
       inc, divSaturation, &
       Densitemolaire, divDensitemolaire, SmDensitemolaire, &
       EnergieInterne, divEnergieInterne, SmEnergieInterne, &
       val, dval, Smval)

    ! input
    type(Type_IncCV), intent(in) :: inc ! contains Saturation

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase), &
         EnergieInterne(NbPhase), &
                                !
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         divEnergieInterne(NbIncPTCSPrimMax, NbPhase), &
         divSaturation(NbIncPTCSPrimMax, NbPhase),     &
                                !
         SmDensiteMolaire(NbPhase),&
         SmEnergieInterne(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncPTCSPrimMax, NbPhase), Smval(NbPhase)

    ! tmp
    integer :: i, iph, k

    val(:) = 0.d0
    dval(:,:) = 0.d0
    Smval(:) = 0.d0

    ! 1. val
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       val(i) = Densitemolaire(i)*EnergieInterne(i)*inc%Saturation(iph)
    end do

    ! 2. dval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       do k=1, NbIncPTCSPrim

          dval(k,i) = &
               divDensitemolaire(k,i)*EnergieInterne(i)*inc%Saturation(iph) &
               + divEnergieInterne(k,i)*DensiteMolaire(i)*inc%Saturation(iph) &
               + divSaturation(k,i)*DensiteMolaire(i)*EnergieInterne(i)
       end do
    end do

    ! 3. Smval
    do i=1, NbPhasePresente
       iph = NumPhasePresente(i)

       Smval(i) = &
            SmDensitemolaire(i)*EnergieInterne(i)*inc%Saturation(iph) &
            + SmEnergieInterne(i)*DensiteMolaire(i)*inc%Saturation(iph)
    end do

  end subroutine LoisThermoHydro_DensitemolaireEnergieInterneSat_cv


  subroutine LoisThermoHydro_DensitemolaireKrViscoEnthalpie_cv( &
       Densitemolaire, divDensiteMolaire, SmDensiteMolaire, &
       PermRel, divPermRel, &
       UnsurViscosite, divUnsurViscosite, SmUnSurViscosite, &
       Enthalpie, divEnthalpie, SmEnthalpie, &
       val, dval, Smval)

    double precision, intent(in) :: &
         DensiteMolaire(NbPhase),   &
         PermRel(NbPhase), &
         UnsurViscosite(NbPhase),   &
                                !
         divDensiteMolaire(NbIncPTCSPrimMax, NbPhase), &
         divPermRel(NbIncPTCSPrimMax, NbPhase), &
         divUnsurViscosite(NbIncPTCSPrimMax, NbPhase), &
                                !
         SmDensiteMolaire(NbPhase), &
         SmUnSurViscosite(NbPhase), &
                                !
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncPTCSPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)

    ! output
    double precision, intent(out) :: &
         val(NbPhase), dval(NbIncPTCSPrimMax, NbPhase), Smval(NbPhase)

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

       do k=1, NbIncPTCSPrim

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

  end subroutine LoisThermoHydro_DensitemolaireKrViscoEnthalpie_cv


  ! Choix des inconnues primaires et secondaires
  !  a partir de la matrice dFsurdX par algorithme glouton
  !  minimisant les angles successifs
  subroutine LoisThermoHydro_IncSecondGlouton(inc, dFsurdX, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! input
    type(Type_IncCV), intent(in) :: inc
    double precision, intent(in) :: dFsurdX(NbIncPTCSMax, NbEqFermetureMax)

    ! output
    integer, intent(out) :: NumIncPTCSPrimCV(NbEqFermetureMax)
    integer, intent(out) :: NumIncPTCSecondCV(NbIncPTCSMax)

    ! tmp
    double precision :: &
         dFsurdXnrm2(NbIncPTCSMax), &
         ctrit(NbIncPTCSMax), ctrit_max

    integer :: ctrit_maxidxs(NbIncPTCSMax)

    double precision :: rnormProjVinc, ss

    double precision :: &
         BaseOrthonormale( NbEqFermetureMax, NbIncPTCSMax)

    integer :: is, i, j, nj, j1

    ! two steps for NumIncPTCSecondcv
    ! 1. first
    ! 2. others

    ! dFsurdXnrm2(i): norm 2 of line i of dFsurdX
    do j=1, NbIncPTCS
       do i=1, NbEqFermeture
          dFsurdXnrm2(j) = dFsurdXnrm2(j) + dFsurdX(j,i)**2
       end do
       dFsurdXnrm2(j) = dsqrt(dFsurdXnrm2(j))
    end do

    ! loop for choosing secd inconnues (size is NbEqFermeture), index: is
    do is=1, NbEqFermeture

       ! compute ctrit
       if(is==1) then

          ctrit(:) = dFsurdXnrm2(:)
       else !

          do j=1, NbIncPTCS ! i: loop index of P T C S

             ss = 0.d0
             do i=1, NbEqFermeture
                ss = ss + BaseOrthonormale(i,j)*dFsurdX(j,i)
             end do

             rnormProjVinc = rnormProjVinc + ss**2
          end do

          rnormProjVinc = sqrt(rnormProjVinc)

          ! maximise distance = rnormeVinc - rnormeProjVinc
          !   where rnormVinc = dFsurdXnrm2(j)
          ctrit(is) = dFsurdXnrm2(is) - rnormProjVinc
       end if

       ! max of ctrit
       ctrit_max = -100.d0
       do j=1, NbIncPTC
          if(ctrit(j)>ctrit_max) then
             ctrit_max = ctrit(j)
          end if
       end do

       ! ctrit_maxidxs: all elements that takes the max value
       nj = 0
       do j=1, NbIncPTCS
          if(abs(ctrit(j)-ctrit_max)<eps) then
             ctrit_maxidxs(nj) = j
             nj = nj + 1
          end if
       end do

       ! which one is second ? (i1, j1)
       NumIncPTCSecondCV(is) = j1

       ! update BaseOrthonomale
       do j=1, NbEqFermeture
          BaseOrthonormale(j,is) = dFsurdX(is,j)
       end do

       do i=1, is-1

          ss = 0.d0
          do j=1, NbEqFermeture
             ss = ss + BaseOrthonormale(j,i)*dFsurdX(is,j)
          end do

          BaseOrthonormale(:,is) = &
               BaseOrthonormale(:,is) - ss*BaseOrthonormale(:,i)
       end do

       ! normalisation
       ss = 0.d0
       do j=1, NbEqFermeture
          ss = ss + BaseOrthonormale(j,is)**2
       enddo

       ss = dsqrt(ss)
       do j =1, NbEqFermeture
          BaseOrthonormale(j,is) = BaseOrthonormale(j,is)/ss
       enddo

    end do ! end loop of is for choosing second inconnues

    ! TODO

  end subroutine LoisThermoHydro_IncSecondGlouton


  !
  subroutine LoisThermoHydro_IncSecondGauss(inc, dFsurdX, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! input
    type(Type_IncCV), intent(in) :: inc
    double precision, intent(in) :: dFsurdX(NbIncPTCSMax, NbEqFermetureMax)

    ! output
    integer, intent(out) :: NumIncPTCSPrimCV(NbIncPTCSPrimMax)
    integer, intent(out) :: NumIncPTCSecondCV(NbEqFermetureMax)

    ! tmps
    double precision :: &
         BB(NbIncPTCSMax, NbEqFermetureMax) ! copy of dFsurdX

    double precision :: pivotmax
    integer :: npivot, pivot(2, NbEqFermetureMax*NbIncPTCSMax) ! ??? size
    logical :: pivot_P, pivot_T

    integer :: & ! E^{eq}, E^{inc}
         NbSetInc, &
         NbSetEq,  &
         NumSetInc(NbIncPTCSMax),   &
         NumSetEq(NbEqFermetureMax)

    logical :: &
         NumIncPTCSPrim_idx(NbIncPTCSMax)

    integer :: is, i, j, numi, numj, n
    integer :: i1, j1, icp, iph, icp1, iph1, k

    ! 1. NumIncPTCSecondcv
    !    1.1 choose secd in (P,T,C), Gauss
    !    1.2 choose secd in S, the first is secd
    ! 2. NumIncPTCSPrimcv

    ! 1.1 choose secd in (P,T,C), Gauss

    ! if(commRank==0) then
    !    do i=1, NbEqFermeture
    !       do j=1, NbIncPTC
    !          print*, dFsurdX(j,i)
    !       end do
    !       print*, ""
    !    end do
    ! end if

    ! init set of inconnus and equations
    NbSetInc = NbIncPTC
    NbSetEq  = NbEqFermeture

    do j=1, NbSetInc
       NumSetInc(j) = j
    end do

    do i=1, NbSetEq
       NumSetEq(i) = i
    end do

    BB(:,:) = dFsurdX(:,:) ! (col, row) index order

    do is=1, NbEqFermeture ! = nb of secd inconnues

       ! max element of abs(BB(:,:))
       pivotmax = -1.d0
       do numi=1, NbSetEq
          i = NumSetEq(numi)

          do numj=1, NbSetInc
             j = NumSetInc(numj)

             if(abs(BB(j,i))>pivotmax) then
                pivotmax = abs(BB(j,i))
             end if
          end do
       end do


       ! set of element (BB) that takes the maximum value
       npivot = 0
       pivot_P = .false.
       pivot_T = .false.

       do numi=1, NbSetEq
          i = NumSetEq(numi)

          do numj=1, NbSetInc
             j = NumSetInc(numj)

             if(abs(abs(BB(j,i))-pivotmax)<eps) then

                npivot = npivot + 1

                pivot(1,npivot) = i ! num eq in BB, not in subset of BB
                pivot(2,npivot) = j ! num inc in BB, not in subset of BB

                ! check if j is P or T
                if(j==1) then
                   pivot_P = .true.
                else if(j==2) then
                   pivot_T = .true.
                end if
             end if

          end do
       end do

       ! if(commRank==1 .and. is==2) then
       !    print*, pivotmax, npivot
       !    do j=1, npivot
       !       print*, pivot(1,j), pivot(2,j)
       !    end do
       !    !print*, ""
       ! end if


       ! choose (i1,j1) is in set pivot(:,),
       ! (i1,j1) is num of BB, not subset of BB

       ! sinon on prend T si elle est dans l'ensemble des pivots max
       if(pivot_T .eqv. .true.) then

          do k=1, npivot
             if(pivot(2,k)==2) then
                i1 = pivot(1,k)
                j1 = 2
             end if
          end do

          NumIncPTCSecondCV(is) = 2 ! j1=2
       else
          ! sinon on prend la plus grande composition dans la phase
          ! avec la plus petite saturation

          i1 = pivot(1,1) ! num of Eq
          j1 = pivot(2,1) ! num of Inc
          icp1 = NumIncPTC2NumIncComp_comp(j1)
          iph1 = NumIncPTC2NumIncComp_phase(j1)

          ! if(is==2 .and. commRank==1) then
          !    print*,"init", i1, j1
          ! end if

          do k=2, npivot

             i = pivot(1,k)
             j = pivot(2,k)
             icp = NumIncPTC2NumIncComp_comp(j)  ! icp in C_{icp}^iph
             iph = NumIncPTC2NumIncComp_phase(j) ! iph in C_{icp}^iph

             !              if(commRank==1 .and. is==2) then
             ! !                print*, i1,j1
             !                 print*, & !i, j ,inc%Saturation(iph), inc%Saturation(iph1), &
             !                      icp, iph, j!, inc%Comp(icp,iph), inc%Comp(icp1, iph1)
             !                 print*, ""
             !              end if

             ! update (i1, j1)
             if(inc%Saturation(iph)<inc%Saturation(iph1)) then ! if < (Saturation)
                i1 = i
                j1 = j
                icp1 = icp
                iph1 = iph

             else if ( ( abs(inc%Saturation(iph)-inc%Saturation(iph1))<eps)) then ! if = (Saturation) and ...
                if (inc%Comp(icp,iph) > inc%Comp(icp1,iph1) )  then !
                   i1 = i
                   j1 = j
                   icp1 = icp
                   iph1 = iph
                end if
             end if

          end do

          NumIncPTCSecondCV(is) = j1
       end if ! end for choosing (i1, j1), NumIncPTCSecondcv(is)=j1

       ! update sets: NumSetInc and NumSetEq, remove (i1, j1)
       n = 0
       do i=1, NbSetEq
          numi = NumSetEq(i)

          if(numi .ne. i1) then
             n = n + 1
             NumSetEq(n) = numi
          end if
       end do
       NbSetEq = n

       n = 0
       do j=1, NbSetInc
          numj = NumSetInc(j)

          if(numj .ne. j1) then
             n = n + 1
             NumSetInc(n) = numj
          end if
       end do
       NbSetInc = n

       ! schur complement
       do numi=1, NbSetEq
          i = NumSetEq(numi)

          do numj=1, NbSetInc
             j = NumSetInc(numj)

             !print*, i1, j1, numi, numj

             BB(j,i) = BB(j,i) &
                  - BB(j1,i)*BB(j,i1)/BB(j1,i1)
          end do
       end do

    end do ! end loop of is

    ! 2. NumIncPTCSPrimcv = {1,2,...,NbIncPTC}/NumIncPTCSecondcv
    !                       + S^alpha, alpha=1:Nphase-1
    NumIncPTCSPrim_idx(:) = .true.
    do j=1, NbEqFermeture
       NumIncPTCSPrim_idx( NumIncPTCSecondCV(j)) = .false. ! not prim
    end do

    n = 0
    do j=1, NbIncPTC
       if(NumIncPTCSPrim_idx(j) .eqv. .true.) then ! prim
          n = n + 1
          NumIncPTCSPrimCV(n) = j
       end if
    end do

    ! last S is secd
    ! if there is only one phase, phase is secd
    do j=NbIncPTC+1, NbIncPTCS-1
       n = n + 1
       NumIncPTCSPrimCV(n) = j
    end do

  end subroutine LoisThermoHydro_IncSecondGauss


  subroutine LoisThermoHydro_divP_wellinj(k)

    integer, intent(in) :: k

    double precision :: Pws, Tw, Sw(NbPhase), Cw(NbComp)
    double precision :: &
         DensiteMolaire, dP_DensiteMolaire, &
         Viscosite, dP_Viscosite, &
         Enthalpie, dP_Enthalpie

    double precision :: dSf(NbPhase), dTf, dCf(NbComp), PermRel
    integer :: s, i

    Sw(:) = 0
    Sw(PHASE_WATER) = 1.d0

    ! node of well k
    do s=NodeDatabyWellInjLocal%Pt(k)+1, NodeDatabyWellInjLocal%Pt(k+1)

       Pws = PerfoWellInj(s)%Pression ! P_{w,s}

       Tw = DataWellInjLocal(k)%Temperature     ! T_w
       Cw(:) = DataWellInjLocal(k)%CompTotal(:) ! C_w

       ! Permrel
       call f_PermRel(PHASE_WATER, Sw, PermRel, dSf)

       ! Molar density
       call f_DensiteMolaire(PHASE_WATER, Pws, Tw, Cw, Sw, &
            DensiteMolaire, dP_DensiteMolaire, dTf, dCf, dSf)

       ! Viscosite
       call f_Viscosite(PHASE_WATER, Pws, Tw, Cw, Sw, &
            Viscosite, dP_Viscosite, dTf, dCf, dSf)

#ifdef _THERMIQUE_
       ! Enthalpie
       call f_Enthalpie(PHASE_WATER, Pws, Tw, Cw, Sw, &
            Enthalpie, dP_Enthalpie, dTf, dCf, dSf)
#endif

       do i=1, NbComp

          ! value
          DensitemolaireKrViscoCompWellInj(i,s) = Cw(i) * PermRel * DensiteMolaire / Viscosite

          ! div of pression
          divDensitemolaireKrViscoCompWellInj(i,s) = Cw(i) * PermRel  &
               * (dP_DensiteMolaire / Viscosite - DensiteMolaire * dP_Viscosite / (Viscosite**2) )
       end do

#ifdef _THERMIQUE_
       DensitemolaireKrViscoEnthalpieWellInj(s) = Enthalpie * PermRel * DensiteMolaire / Viscosite

       divDensitemolaireKrViscoEnthalpieWellInj(s) = &
            + dP_DensiteMolaire / Viscosite * Enthalpie &
            + dP_Enthalpie * DensiteMolaire / Viscosite &
            - dP_Viscosite / (Viscosite**2) * DensiteMolaire * Enthalpie
#endif
    end do

  end subroutine LoisThermoHydro_divP_wellinj



  ! allocate
  subroutine LoisThermoHydro_allocate

    integer :: nbCell, nbFrac, nbNode, nbNodeInj

    nbCell = NbCellLocal_Ncpus(commRank+1)
    nbFrac = NbFracLocal_Ncpus(commRank+1)
    nbNode = NbNodeLocal_Ncpus(commRank+1)
    nbNodeInj = NodeByWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)
    ! print*, 'LoisThermoHydro_allocate', nbCell, nbFrac, nbNode, nbNodeInj

    ! dXssurdXp and SmdXs
    allocate( dXssurdXpCell(NbIncPTCSPrimMax, NbEqFermetureMax, nbCell))
    allocate( dXssurdXpFrac(NbIncPTCSPrimMax, NbEqFermetureMax, nbFrac))
    allocate( dXssurdXpNode(NbIncPTCSPrimMax, NbEqFermetureMax, nbNode))

    allocate( SmdXsCell(NbEqFermetureMax, nbCell))
    allocate( SmdXsFrac(NbEqFermetureMax, nbFrac))
    allocate( SmdXsNode(NbEqFermetureMax, nbNode))

    allocate( SmFCell(NbEqFermetureMax, nbCell))
    allocate( SmFFrac(NbEqFermetureMax, nbFrac))
    allocate( SmFNode(NbEqFermetureMax, nbNode))

    ! Num IncPTCSPrim and IncPTCSecond
    allocate( NumIncPTCSPrimCell(NbIncPTCSPrimMax, nbCell))
    allocate( NumIncPTCSPrimFrac(NbIncPTCSPrimMax, nbFrac))
    allocate( NumIncPTCSPrimNode(NbIncPTCSPrimMax, nbNode))
    allocate( NumIncPTCSecondCell(NbEqFermetureMax,nbCell ))
    allocate( NumIncPTCSecondFrac(NbEqFermetureMax,nbFrac ))
    allocate( NumIncPTCSecondNode(NbEqFermetureMax,nbNode ))

    ! deniste massique
    allocate( DensiteMassiqueCell(NbPhase, nbCell))
    allocate( DensiteMassiqueFrac(NbPhase, nbFrac))
    allocate( DensiteMassiqueNode(NbPhase, nbNode))

    allocate( divDensiteMassiqueCell(NbIncPTCSPrimMax, NbPhase, nbCell))
    allocate( divDensiteMassiqueFrac(NbIncPTCSPrimMax, NbPhase, nbFrac))
    allocate( divDensiteMassiqueNode(NbIncPTCSPrimMax, NbPhase, nbNode))

    allocate( SmDensiteMassiqueCell(NbPhase, nbCell))
    allocate( SmDensiteMassiqueFrac(NbPhase, nbFrac))
    allocate( SmDensiteMassiqueNode(NbPhase, nbNode))

    ! pression
    allocate( divPressionCell(NbIncPTCSPrimMax, nbCell) )
    allocate( divPressionFrac(NbIncPTCSPrimMax, nbFrac) )
    allocate( divPressionNode(NbIncPTCSPrimMax, nbNode) )

    allocate( SmPressionCell( nbCell))
    allocate( SmPressionFrac( nbFrac))
    allocate( SmPressionNode( nbNode))

    ! pression capillaire
    allocate( PressionCapCell(NbPhase, nbCell))
    allocate( PressionCapFrac(NbPhase, nbFrac))
    allocate( PressionCapNode(NbPhase, nbNode))

    allocate( divPressionCapCell(NbIncPTCSPrimMax, NbPhase, nbCell))
    allocate( divPressionCapFrac(NbIncPTCSPrimMax, NbPhase, nbFrac))
    allocate( divPressionCapNode(NbIncPTCSPrimMax, NbPhase, nbNode))

    ! Saturation
    allocate( divSaturationCell(NbIncPTCSPrimMax, NbPhase, nbCell))
    allocate( divSaturationFrac(NbIncPTCSPrimMax, NbPhase, nbFrac))
    allocate( divSaturationNode(NbIncPTCSPrimMax, NbPhase, nbNode))

    ! DensiteMolaire*Kr/Viscosite * Comp
    allocate( DensitemolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
    allocate( DensitemolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
    allocate( DensitemolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

    allocate( divDensitemolaireKrViscoCompCell(NbIncPTCSPrimMax, NbComp, NbPhase, nbCell))
    allocate( divDensitemolaireKrViscoCompFrac(NbIncPTCSPrimMax, NbComp, NbPhase, nbFrac))
    allocate( divDensitemolaireKrViscoCompNode(NbIncPTCSPrimMax, NbComp, NbPhase, nbNode))

    allocate( SmDensitemolaireKrViscoCompCell(NbComp, NbPhase, nbCell))
    allocate( SmDensitemolaireKrViscoCompFrac(NbComp, NbPhase, nbFrac))
    allocate( SmDensitemolaireKrViscoCompNode(NbComp, NbPhase, nbNode))

    ! DensiteMolaire * Saturation * Comp
    allocate( DensitemolaireSatCompCell(NbComp, NbPhase, nbCell))
    allocate( DensitemolaireSatCompFrac(NbComp, NbPhase, nbFrac))
    allocate( DensitemolaireSatCompNode(NbComp, NbPhase, nbNode))

    allocate( divDensitemolaireSatCompCell(NbIncPTCSPrimMax, NbComp, NbPhase, nbCell))
    allocate( divDensitemolaireSatCompFrac(NbIncPTCSPrimMax, NbComp, NbPhase, nbFrac))
    allocate( divDensitemolaireSatCompNode(NbIncPTCSPrimMax, NbComp, NbPhase, nbNode))

    allocate( SmDensitemolaireSatCompCell(NbComp, NbPhase, nbCell))
    allocate( SmDensitemolaireSatCompFrac(NbComp, NbPhase, nbFrac))
    allocate( SmDensitemolaireSatCompNode(NbComp, NbPhase, nbNode))

    ! well inj
    allocate(DensitemolaireKrViscoCompWellInj(NbComp, nbNodeInj))
    allocate(divDensitemolaireKrViscoCompWellInj(NbComp, nbNodeInj))

    allocate(DensitemolaireKrViscoEnthalpieWellInj(nbNodeInj))
    allocate(divDensitemolaireKrViscoEnthalpieWellInj(nbNodeInj))

#ifdef _THERMIQUE_

    ! temperature
    allocate( divTemperatureCell(NbIncPTCSPrimMax, nbCell) )
    allocate( divTemperatureFrac(NbIncPTCSPrimMax, nbFrac) )
    allocate( divTemperatureNode(NbIncPTCSPrimMax, nbNode) )

    allocate( SmTemperatureCell( nbCell))
    allocate( SmTemperatureFrac( nbFrac))
    allocate( SmTemperatureNode( nbNode))


    ! DensiteMolaire * PermRel / Viscosite * Enthalpie
    allocate( DensitemolaireKrViscoEnthalpieCell(NbPhase, nbCell))
    allocate( DensitemolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
    allocate( DensitemolaireKrViscoEnthalpieNode(NbPhase, nbNode))

    allocate( divDensitemolaireKrViscoEnthalpieCell(NbIncPTCSPrimMax, NbPhase, nbCell))
    allocate( divDensitemolaireKrViscoEnthalpieFrac(NbIncPTCSPrimMax, NbPhase, nbFrac))
    allocate( divDensitemolaireKrViscoEnthalpieNode(NbIncPTCSPrimMax, NbPhase, nbNode))

    allocate( SmDensitemolaireKrViscoEnthalpieCell(NbPhase, nbCell))
    allocate( SmDensitemolaireKrViscoEnthalpieFrac(NbPhase, nbFrac))
    allocate( SmDensitemolaireKrViscoEnthalpieNode(NbPhase, nbNode))

    ! densitemolaire * energieinterne * Saturation
    allocate( DensitemolaireEnergieInterneSatCell(NbPhase, nbCell))
    allocate( DensitemolaireEnergieInterneSatFrac(NbPhase, nbFrac))
    allocate( DensitemolaireEnergieInterneSatNode(NbPhase, nbNode))

    allocate( divDensitemolaireEnergieInterneSatCell(NbIncPTCSPrimMax, NbPhase, nbCell))
    allocate( divDensitemolaireEnergieInterneSatFrac(NbIncPTCSPrimMax, NbPhase, nbFrac))
    allocate( divDensitemolaireEnergieInterneSatNode(NbIncPTCSPrimMax, NbPhase, nbNode))

    allocate( SmDensitemolaireEnergieInterneSatCell(NbPhase, nbCell))
    allocate( SmDensitemolaireEnergieInterneSatFrac(NbPhase, nbFrac))
    allocate( SmDensitemolaireEnergieInterneSatNode(NbPhase, nbNode))

#endif

  end subroutine LoisThermoHydro_allocate


  ! free
  subroutine LoisThermoHydro_free

    ! dXssurdXp
    deallocate( dXssurdXpCell)
    deallocate( dXssurdXpFrac)
    deallocate( dXssurdXpNode)

    ! Smdxs
    deallocate( SmdXsCell)
    deallocate( SmdXsFrac)
    deallocate( SmdXsNode)

    ! SmF
    deallocate( SmFCell)
    deallocate( SmFFrac)
    deallocate( SmFNode)


    ! Num IncPTCSPrim and IncPTCSecond
    deallocate( NumIncPTCSPrimCell)
    deallocate( NumIncPTCSPrimFrac)
    deallocate( NumIncPTCSPrimNode)
    deallocate( NumIncPTCSecondCell)
    deallocate( NumIncPTCSecondFrac)
    deallocate( NumIncPTCSecondNode)

    ! densitemassique
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

    ! densitemolaire * Permrel / viscosite * Comp
    deallocate( DensitemolaireKrViscoCompCell)
    deallocate( DensitemolaireKrViscoCompFrac)
    deallocate( DensitemolaireKrViscoCompNode)
    deallocate( divDensitemolaireKrViscoCompCell)
    deallocate( divDensitemolaireKrViscoCompFrac)
    deallocate( divDensitemolaireKrViscoCompNode)
    deallocate( SmDensitemolaireKrViscoCompCell)
    deallocate( SmDensitemolaireKrViscoCompFrac)
    deallocate( SmDensitemolaireKrViscoCompNode)

    ! densitemolaire * Sat * Comp
    deallocate( DensitemolaireSatCompCell)
    deallocate( DensitemolaireSatCompFrac)
    deallocate( DensitemolaireSatCompNode)
    deallocate( divDensitemolaireSatCompCell)
    deallocate( divDensitemolaireSatCompFrac)
    deallocate( divDensitemolaireSatCompNode)
    deallocate( SmDensitemolaireSatCompCell)
    deallocate( SmDensitemolaireSatCompFrac)
    deallocate( SmDensitemolaireSatCompNode)

    ! well inj
    deallocate(DensitemolaireKrViscoCompWellInj)
    deallocate(divDensitemolaireKrViscoCompWellInj)

    deallocate(DensitemolaireKrViscoEnthalpieWellInj)
    deallocate(divDensitemolaireKrViscoEnthalpieWellInj)

#ifdef _THERMIQUE_
    ! temperature
    deallocate( divTemperatureCell)
    deallocate( divTemperatureFrac)
    deallocate( divTemperatureNode)
    deallocate( SmTemperatureCell)
    deallocate( SmTemperatureFrac)
    deallocate( SmTemperatureNode)

    ! densitemolaire * Permrel / viscosite * Enthalpie
    deallocate( DensitemolaireKrViscoEnthalpieCell)
    deallocate( DensitemolaireKrViscoEnthalpieFrac)
    deallocate( DensitemolaireKrViscoEnthalpieNode)
    deallocate( divDensitemolaireKrViscoEnthalpieCell)
    deallocate( divDensitemolaireKrViscoEnthalpieFrac)
    deallocate( divDensitemolaireKrViscoEnthalpieNode)
    deallocate( SmDensitemolaireKrViscoEnthalpieCell)
    deallocate( SmDensitemolaireKrViscoEnthalpieFrac)
    deallocate( SmDensitemolaireKrViscoEnthalpieNode)

    ! densitemolaire * energieinterne * Saturation
    deallocate( DensitemolaireEnergieInterneSatCell)
    deallocate( DensitemolaireEnergieInterneSatFrac)
    deallocate( DensitemolaireEnergieInterneSatNode)
    deallocate( divDensitemolaireEnergieInterneSatCell)
    deallocate( divDensitemolaireEnergieInterneSatFrac)
    deallocate( divDensitemolaireEnergieInterneSatNode)
    deallocate( SmDensitemolaireEnergieInterneSatCell)
    deallocate( SmDensitemolaireEnergieInterneSatFrac)
    deallocate( SmDensitemolaireEnergieInterneSatNode)

    ! ! densitemolaire * energieinterne
    ! deallocate( DensitemolaireEnergieInterneCell)
    ! deallocate( DensitemolaireEnergieInterneFrac)
    ! deallocate( DensitemolaireEnergieInterneNode)
    ! deallocate( divDensitemolaireEnergieInterneCell)
    ! deallocate( divDensitemolaireEnergieInterneFrac)
    ! deallocate( divDensitemolaireEnergieInterneNode)
    ! deallocate( SmDensitemolaireEnergieInterneCell)
    ! deallocate( SmDensitemolaireEnergieInterneFrac)
    ! deallocate( SmDensitemolaireEnergieInterneNode)
#endif

  end subroutine LoisThermoHydro_free


  ! compute prim values using secd values
  ! secd = SmdX - dXssurdXp * prim
  subroutine LoisThermoHydro_PrimToSecd( &
       vnode, vfrac, vcell)    ! prim and secd

    double precision, dimension(:,:), intent(inout) :: &
         vnode, vfrac, vcell

    double precision :: &
         xp(NbCompThermique), &
         xs(NbEqFermetureMax)

    integer :: k, ic, i, iph
    integer :: &
         NbPhasePresente, NbEqFermeture, &
         NbIncPTC, NbIncPTCS, NbIncPTCPrim, NbIncPTCSPrim

    ! node
    do k=1, NbNodeLocal_Ncpus(commRank+1)

       ic = IncNode(k)%ic
       NbPhasePresente = NbPhasePresente_ctx(ic)
       NbEqFermeture = NbEqFermeture_ctx(ic)
       NbIncPTC  = NbIncPTC_ctx(ic)
       NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente
       NbIncPTCSPrim = NbIncPTCSPrim_ctx(ic)
       NbIncPTCPrim = NbIncPTC - NbEqFermeture

       xp(1:NbCompThermique) = vnode(1:NbCompThermique,k)
       xs(1:NbEqFermeture) = SmdXsNode(1:NbEqFermeture,k)

       call dgemv('T', NbIncPTCSPrim, NbEqFermeture, &
            -1.d0, dXssurdXpNode(:,:,k), NbIncPTCSPrimMax, &
            xp(:), 1, -1.d0, xs(:), 1)

       vnode(:,k) = 0.d0

       ! copy prim P,T,C,S
       do i=1, NbIncPTCSPrim
          vnode( NumIncPTCSPrimNode(i,k),k) = xp(i)
       end do

       ! copy secd P,T,C
       do i=1, NbEqFermeture
          vnode( NumIncPTCSecondNode(i,k),k) = xs(i)
       end do

       ! if NbPhasePresente=1,
       !    then this phase must be secd
       ! else last saturation is secd, the others are prim
       !    secd S = - sum_{S^alpha is prim} S^alpha
       if( NbPhasePresente==1) then

          ! iph is this phase in (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(1,ic)
          vnode(iph,k) = 0.d0
       else

          ! iph is last present phase in vector (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente,ic)
          vnode(iph,k) = 0.d0
          do i=1, NbPhasePresente-1
             vnode(iph,k) = vnode(iph,k) - xp(NbIncPTCPrim+i)
          end do
       end if

       ! term prim n_k(X_j^n)
       do i=1, NbCompCtilde_ctx(ic)
          vnode(NbIncPTCS+i,k) = xp(NbIncPTCSPrim+i)
       end do
    end do

    ! frac
    do k=1, NbFracLocal_Ncpus(commRank+1)

       ic = IncFrac(k)%ic
       NbPhasePresente = NbPhasePresente_ctx(ic)
       NbEqFermeture = NbEqFermeture_ctx(ic)
       NbIncPTC  = NbIncPTC_ctx(ic)
       NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente
       NbIncPTCPrim = NbIncPTC - NbEqFermeture
       NbIncPTCSPrim = NbIncPTCSPrim_ctx(ic)

       xp(1:NbCompThermique) = vfrac(1:NbCompThermique,k)
       xs(1:NbEqFermeture) = SmdXsFrac(1:NbEqFermeture,k)

       call dgemv('T', NbIncPTCSPrim, NbEqFermeture, &
            -1.d0, dXssurdXpFrac(:,:,k), NbIncPTCSPrimMax, &
            xp(:), 1, -1.d0, xs(:), 1)

       vfrac(:,k) = 0.d0

       ! copy prim P,T,C,S
       do i=1, NbIncPTCSPrim
          vfrac( NumIncPTCSPrimFrac(i,k),k) = xp(i)
       end do

       ! copy secd P,T,C
       do i=1, NbEqFermeture
          vfrac( NumIncPTCSecondFrac(i,k),k) = xs(i)
       end do

       ! if NbPhasePresente=1,
       !    then this phase must be secd
       ! else last saturation is secd, the others are prim
       !    secd S = - sum_{S^alpha is prim} S^alpha
       if( NbPhasePresente==1) then

          ! iph is this phase in (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(1,ic)
          vfrac(iph,k) = 0.d0
       else

          ! iph is last present phase in vector (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente,ic)
          vfrac(iph,k) = 0.d0
          do i=1, NbPhasePresente-1
             vfrac(iph,k) = vfrac(iph,k) - xp(NbIncPTCPrim+i)
          end do
       end if

       ! term n_k(X_j^n)
       do i=1, NbCompCtilde_ctx(ic) ! =NbCompThermique-NbIncPTCSPrim
          vfrac(NbIncPTCS+i,k) = xp(NbIncPTCSPrim+i)
       end do
    end do

    ! cell
    do k=1, NbCellLocal_Ncpus(commRank+1)

       ic = IncCell(k)%ic
       NbPhasePresente = NbPhasePresente_ctx(ic)
       NbEqFermeture = NbEqFermeture_ctx(ic)
       NbIncPTC  = NbIncPTC_ctx(ic)
       NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente
       NbIncPTCPrim = NbIncPTC - NbEqFermeture
       NbIncPTCSPrim = NbIncPTCSPrim_ctx(ic)

       xp(1:NbCompThermique) = vcell(1:NbCompThermique,k)
       xs(1:NbEqFermeture) = SmdXsCell(1:NbEqFermeture,k)

       call dgemv('T', NbIncPTCSPrim, NbEqFermeture, &
            -1.d0, dXssurdXpCell(:,:,k), NbIncPTCSPrimMax, &
            xp(:), 1, -1.d0, xs(:), 1)

       vcell(:,k) = 0.d0

       ! copy prim P,T,C,S
       do i=1, NbIncPTCSPrim
          vcell( NumIncPTCSPrimCell(i,k),k) = xp(i)
       end do

       ! copy secd P,T,C
       do i=1, NbEqFermeture
          vcell( NumIncPTCSecondCell(i,k),k) = xs(i)
       end do

       ! if NbPhasePresente=1,
       !    then this phase must be secd
       ! else last saturation is secd, the others are prim
       !    secd S = - sum_{S^alpha is prim} S^alpha
       if( NbPhasePresente==1) then

          ! iph is this phase in (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(1,ic)
          vcell(iph,k) = 0.d0

       else

          ! iph is last present phase in vector (P,T,C,S,n_i)
          iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente,ic)
          vcell(iph,k) = 0.d0
          do i=1, NbPhasePresente-1
             vcell(iph,k) = vcell(iph,k) - xp(NbIncPTCPrim+i)
          end do
       end if

       ! term n_k(X_j^n)
       do i=1, NbCompCtilde_ctx(ic)
          vcell(NbIncPTCS+i,k) = xp(NbIncPTCSPrim+i)
       end do
    end do

  end subroutine LoisThermoHydro_PrimToSecd


end module LoisThermoHydro
