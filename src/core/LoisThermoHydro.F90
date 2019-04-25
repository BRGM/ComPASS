!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LoisThermoHydro

  use DefModel
  use Thermodynamics
  use NumbyContext
  use IncCVReservoir
  use IncCVWells
  use IncPrimSecd

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
                                !
       LoisThermoHydro_fill_gradient_dfdX,     & ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       LoisThermoHydro_dfdX_ps,                & ! fill dfdX_prim/dfdX_secd with the derivatives w.r.t. the primary/secondary unknowns 
       LoisThermoHydro_Densitemolaire_cv,      & ! prim divs: densitemolaire
       LoisThermoHydro_Viscosite_cv,           & !            1/viscosite
       LoisThermoHydro_PermRel_cv,             & !            Permrel
       LoisThermoHydro_Inc_cv,                 & !            called with Pression and Temperature
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

    integer :: k

    ! cell
    do k=1, NbCellLocal_Ncpus(commRank+1)

       call LoisThermoHydro_divPrim_cv(IncCell(k), CellRocktypeLocal(:,k), &
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

       call LoisThermoHydro_divPrim_cv(IncFrac(k), FracRocktypeLocal(:,k), &
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



       call LoisThermoHydro_divPrim_cv(IncNode(k), NodeRocktypeLocal(:,k), &
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
  subroutine LoisThermoHydro_divPrim_cv(inc, rt, &
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
    type(TYPE_IncCVReservoir), intent(in) :: inc

    integer, intent(in) :: rt (IndThermique+1)

    integer, intent(in) ::  &
         NumIncPTCSPrimCV (NbIncPTCSPrimMax),  &
         NumIncPTCSecondCV (NbEqFermetureMax)

    double precision, intent(in) :: &
         dXssurdXp (NbIncPTCSPrimMax, NbEqFermetureMax), & ! (col,row) index order
         SmdXs (NbEqFermetureMax), &
         SmF (NbEqFermetureMax)

    ! output

    double precision, intent(out) :: &
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
    call LoisThermoHydro_PermRel_cv(inc, rt, PermRel, divPermRel)

    ! Pression (unknown index is 1)
    call LoisThermoHydro_Inc_cv(1, inc, &
         NumIncPTCSPrimCV, NumIncPTCSecondCV, &
         dXssurdXp, SmdXs,  &
         divPression, SmPression)

    ! Pression Capillaire
    call LoisThermoHydro_PressionCapillaire_cv(rt, inc, PressionCap, divPressionCap)

    ! Saturation div
    call LoisThermoHydro_Saturation_cv(inc, &
      NumIncPTCSPrimCV, NumIncPTCSecondCV, &
      dXssurdXp, &
      divSaturation)

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

    ! Temperature (unknown index is 2)
    call LoisThermoHydro_Inc_cv(2, inc, &
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

    integer :: k

    do k=1, NbNodeLocal_Ncpus(commRank+1)

       call LoisThermoHydro_divPrim_cv(IncNode(k), NodeRocktypeLocal(:,k), &
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

    type(TYPE_IncCVReservoir), intent(in) :: inc

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

  ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
  subroutine LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)
   integer, intent(in) :: iph ! num of phase
    double precision, intent(in) :: dPf, dTf, dCf(NbComp), dSf(NbPhase)
    double precision, intent(out) :: dfdX(NbIncPTCSMax)

    ! tmp
    integer :: j, jc

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
    
    do j=1, NbPhase ! S
         jc = j + NbIncPTC
         dfdX(jc) = dSf(j)
    enddo
    
  end subroutine LoisThermoHydro_fill_gradient_dfdX

  subroutine LoisThermoHydro_dfdX_ps(iph, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
       dfdX_prim, dfdX_secd)

    ! input 
    integer, intent(in) :: iph ! num of phase
    integer, intent(in) :: &
      NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
      NumIncPTCSecondCV(NbEqFermetureMax)
    double precision, intent(in) :: dfdX(NbIncPTCSMax)  ! dfdX = (df/dP, df/dT, df/dC, df/dS)
    ! output
    double precision, intent(out) :: &
      dfdX_prim(NbIncPTCSPrimMax, NbPhase), &
      dfdX_secd(NbIncPTCSecondMax, NbPhase)

    ! tmp
      integer :: j
  
    ! prim and secd part of dfdX
    do j=1, NbIncPTCSPrim
      dfdX_prim(j,iph) = dfdX(NumIncPTCSPrimCV(j))
    end do

    do j=1, NbIncPTCSecond ! = NbEqFermeture
         dfdX_secd(j,iph) = dfdX(NumIncPTCSecondCV(j))
    end do

  end subroutine LoisThermoHydro_dfdX_ps

  subroutine LoisThermoHydro_densitemassique_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
       call LoisThermoHydro_dfdX_ps(iph, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
            dval, dfdX_secd)
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
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
    if(f.leq.0) then
        call CommonMPI_abort('inconsistent viscosity value')
    endif
#endif

       val(i) = 1.d0/f ! val

       ! fill dfdX = (df/dP, df/dT, df/dC, df/dS)
       call LoisThermoHydro_fill_gradient_dfdX(iph, dPf, dTf, dCf, dSf, dfdX)
       dfdX(:) = -dfdX(:)/f**2

       ! fill dval with the derivatives w.r.t. the primary unknowns (dval=dfdX_prim)
       ! and dfdX_secd w.r.t. the secondary unknowns 
       call LoisThermoHydro_dfdX_ps(i, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
            dval, dfdX_secd)
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
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
       call LoisThermoHydro_dfdX_ps(i, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
            dval, dfdX_secd)
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



  subroutine LoisThermoHydro_PermRel_cv(inc, rt, val, dval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
    INTEGER, INTENT(IN) :: rt(IndThermique+1)

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

       call f_PermRel(rt, iph, inc%Saturation, f, dSf)

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

  subroutine LoisThermoHydro_Inc_cv(index_inc, inc, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       dXssurdXp, SmdXs, &
       dval, Smval)

    ! input
    integer, intent(in) :: index_inc
    type(TYPE_IncCVReservoir), intent(in)  :: inc
    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV( NbEqFermetureMax)

    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: dval(NbIncPTCSPrimMax)
    double precision, intent(out) :: Smval

    integer :: i

    dval(:) = 0.d0

    IF(ANY(NumIncPTCSPrimCV == index_inc))THEN ! index_inc (P or T) is prim
      do i=1,NbIncPTCSPrimMax
        if (NumIncPTCSPrimCV(i) == index_inc) then
          dval(i) = 1.d0
          Smval = 0.d0
        endif
      enddo
    ELSE IF(ANY(NumIncPTCSecondCV == index_inc))THEN ! index_inc (P or T) is secd
      do i=1,NbEqFermeture
        if (NumIncPTCSecondCV(i) == index_inc) then
          dval(:) = - dXssurdXp(:,i)
          Smval = - SmdXs(i)
        endif
      enddo
    ELSE ! index_inc (P or T) not found
      if(index_inc == 1) write(*,*)' pb in derprim, P not found '
      if(index_inc == 2) write(*,*)' pb in derprim, T not found '
      write(*,*)' primary unknown'
      write(*,*) NumIncPTCSPrimCV
      write(*,*)' secondary unknown'
      write(*,*) NumIncPTCSecondCV
      stop
    ENDIF

  end subroutine LoisThermoHydro_Inc_cv



  subroutine LoisThermoHydro_Saturation_cv(inc, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       dXssurdXp, &
       dval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
    integer, intent(in) :: &
         NumIncPTCSPrimCV(NbIncPTCSPrimMax), &
         NumIncPTCSecondCV( NbEqFermetureMax)
    double precision, intent(in) :: & ! (col, row) index order
         dXssurdXp(NbIncPTCSPrimMax, NbEqFermetureMax)

    ! output
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)

    ! tmp
    integer :: i, iph, k

    dval(:,:) = 0.d0

    ! secondary staturation is the first non present phase
    DO i=1, NbPhasePresente - 1
      iph = NumPhasePresente(i)

      IF(ANY(NumIncPTCSPrimCV == i+NbIncPTC))THEN ! S is prim
        do k=1,NbIncPTCSPrimMax
          if (NumIncPTCSPrimCV(k) == i+NbIncPTC) then
            dval(k,i) = 1.d0

            dval(k,NbPhasePresente) = -1.d0
          endif
        enddo
      ELSE IF(ANY(NumIncPTCSecondCV == i+NbIncPTC))THEN ! S is secd
        do k=1,NbEqFermeture
          if (NumIncPTCSecondCV(k) == i+NbIncPTC) then
            dval(:,i) = - dXssurdXp(:,k)

            dval(:,NbPhasePresente) = dXssurdXp(:,k)
          endif
        enddo
      ELSE ! S not found
        write(*,*)' pb dans derprim, S non trouvee '
        write(*,*)' primary unknown'
        write(*,*) NumIncPTCSPrimCV
        write(*,*)' secondary unknown'
        write(*,*) NumIncPTCSecondCV
        write(*,*)' saturation'
        write(*,*) i
        stop
      ENDIF
    ENDDO

  end subroutine LoisThermoHydro_Saturation_cv



  subroutine LoisThermoHydro_PressionCapillaire_cv(rt, inc, val, dval)

    ! input
    integer, intent(in) :: rt(IndThermique+1)
    type(TYPE_IncCVReservoir), intent(in)  :: inc

    ! output
    double precision, intent(out) :: val(NbPhase)
    double precision, intent(out) :: dval(NbIncPTCSPrimMax, NbPhase)

    ! tmp
    double precision :: f, dSf(NbPhase), dfS_secd
    integer :: i, iph, j, jph

    val(:) = 0.d0
    dval(:,:) = 0.d0

    do iph=1, NbPhase

       call f_PressionCapillaire(rt, iph, inc%Saturation, f, dSf)

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
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
       call LoisThermoHydro_dfdX_ps(i, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
            dval, dfdX_secd)
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

  subroutine LoisThermoHydro_Enthalpie_cv(inc, dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, val, dval, Smval)

    ! input
    type(TYPE_IncCVReservoir), intent(in)  :: inc
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
       call LoisThermoHydro_dfdX_ps(i, NumIncPTCSPrimCV, NumIncPTCSecondCV, dfdX, &
            dval, dfdX_secd)
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
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

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
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

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
    type(TYPE_IncCVReservoir), intent(in) :: inc

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
    type(TYPE_IncCVReservoir), intent(in) :: inc ! contains Saturation

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


  subroutine LoisThermoHydro_divP_wellinj(k)

    integer, intent(in) :: k

    double precision :: Pws, Tw, Sw(NbPhase), Cw(NbComp)
    double precision :: &
         DensiteMolaire, dP_DensiteMolaire, &
         Viscosite, dP_Viscosite, &
         Enthalpie, dP_Enthalpie

    double precision :: dSf(NbPhase), dTf, dCf(NbComp), PermRel
    integer :: rt(IndThermique+1)
    integer :: s, i

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

    ! densite massique
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

    ! the following arrays must be allocated even if there is no energy transfer
    allocate( SmTemperatureCell( nbCell))
    allocate( SmTemperatureFrac( nbFrac))
    allocate( SmTemperatureNode( nbNode))

#ifdef _THERMIQUE_

    ! temperature
    allocate( divTemperatureCell(NbIncPTCSPrimMax, nbCell) )
    allocate( divTemperatureFrac(NbIncPTCSPrimMax, nbFrac) )
    allocate( divTemperatureNode(NbIncPTCSPrimMax, nbNode) )

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

    deallocate( SmTemperatureCell)
    deallocate( SmTemperatureFrac)
    deallocate( SmTemperatureNode)

#ifdef _THERMIQUE_
    ! temperature
    deallocate( divTemperatureCell)
    deallocate( divTemperatureFrac)
    deallocate( divTemperatureNode)

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

end module LoisThermoHydro
