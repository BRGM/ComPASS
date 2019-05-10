!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecd

  use iso_c_binding, only: c_double
  use mpi, only: MPI_Abort
  use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

  use SchemeParameters, only: eps

  use DefModel, only: &
     NbComp, NbPhase, MCP, &
     NbIncPTCSPrimMax, NbEqFermetureMax, NbIncPTCSMax, &
     pschoice, psprim, pssecd, &
     NumPhasePresente_ctx, NbPhasePresente_ctx, &
     NbIncPTCMax, IndThermique, NbCompThermique, NbEqEquilibreMax

  use NumbyContext, only: &
     NbEqFermeture_ctx, NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, &
     NumIncComp2NumIncPTC_ctx, NbIncPTCSPrim_ctx, NumCompCtilde_ctx, &
     NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
     NbEqEquilibre_ctx, NbIncPTC_ctx, NbCompCtilde_ctx

  use IncCVReservoir, only: &
     TYPE_IncCVReservoir, &
     IncCell, IncFrac, IncNode
  use MeshSchema, only: &
     NbCellLocal_Ncpus, NbFracLocal_Ncpus, NbNodeLocal_Ncpus, &
     NodeRocktypeLocal, CellRocktypeLocal, FracRocktypeLocal, &
     NodeByWellInjLocal

  use Newton, only: Newton_increments_pointers, Newton_increments, Newton_pointers_to_values
  use Thermodynamics, only: f_Fugacity

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


  ! tmp values to simpfy notations of numerotation
  ! ex. NbPhasePresente = NbPhasePresente_ctx(inc%ic)
  integer, private :: &
       NbPhasePresente, NbCompCtilde, &
       NbEqFermeture, NbEqEquilibre,  & ! NbEqEquilibre -> Nombre d'Equation d'Equilibre thermodynamique fct du contexte (i.e. égalité des fugacités)
       NbIncPTC, NbIncPTCPrim, NbIncPTCS, &
       NbIncPTCSPrim, NbIncPTCSecond, &
                                !
       NumPhasePresente(NbPhase),               & ! Num -> identifiants de la phase présente
       NumCompCtilde(NbComp),                   & ! Num -> identifiants des composants absents
       NumCompEqEquilibre(NbEqEquilibreMax),    & ! identifiant des composant présents dans au moins 2 phases (donc concernés par égalité fugacités)
       NumIncPTC2NumIncComp_comp(NbIncPTCMax),  & ! Etant donné une ligne du "vecteur inconnu" quel est le composant
       NumIncPTC2NumIncComp_phase(NbIncPTCMax), &
                                ! Etant donné une ligne du "vecteur inconnu" quel est la phase ????
       Num2PhasesEqEquilibre(2, NbEqEquilibreMax), & ! phases impliquées dans l'équilibre (cf. NumCompEqEquilibre)
       NumIncComp2NumIncPTC(NbComp, NbPhase) ! matrice donnant pour chaque phase et chaque composant la ligne du "vecteur inconnu" 

  public :: &
       IncPrimSecd_allocate,   &
       IncPrimSecd_free,       &
       IncPrimSecd_compute, &
       IncPrimSecd_compPrim_nodes


  private :: &
       IncPrimSecd_compute_cv,             & ! all operations for one cv
       IncPrimSecd_init_cv,                & ! init infos according to ic (context) for each control volume (cv)
       IncPrimSecd_dFsurdX_cv,             & ! compute dF/dX for each cv
       IncPrimSecd_ps_cv,                  & ! choose first and second variables for each cv
       IncPrimSecd_dXssurdXp_cv              ! compute dFs/dXp
                  
contains

  ! main subroutine of this module
  ! compute dXssurdXp, SmdXs, NumIncPTCSPrim, NumIncPTCSecond
  ! loop of cell/frac/node
  subroutine IncPrimSecd_compute() &
        bind(C, name="IncPrimSecd_compute")

    integer :: k

    ! cell
    call IncPrimSecd_compute_cv( &
         NbCellLocal_Ncpus(commRank+1), &
         IncCell, CellRocktypeLocal, &
         dXssurdXpCell, SmdXsCell, &
         SmFCell,   &
                                !
         NumIncPTCSPrimCell, NumIncPTCSecondCell)

    ! frac
    call IncPrimSecd_compute_cv( &
         NbFracLocal_Ncpus(commRank+1), &
         IncFrac, FracRocktypeLocal, &
         dXssurdXpFrac, SmdXsFrac, &
         SmFFrac,   &
                                !
         NumIncPTCSPrimFrac, NumIncPTCSecondFrac)

    ! node
    call IncPrimSecd_compute_cv( &
         NbNodeLocal_Ncpus(commRank+1), &
         IncNode, NodeRocktypeLocal, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode,   &
                                !
         NumIncPTCSPrimNode, NumIncPTCSecondNode)


  end subroutine IncPrimSecd_compute


  ! all operations for a set of cv
  subroutine IncPrimSecd_compute_cv( &
       NbIncLocal, &
       inc, rt, &
       dXssurdXp, SmdXs, SmF, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! input
    integer, intent(in) :: NbIncLocal

    type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

    integer, intent(in) :: rt(IndThermique+1, NbIncLocal)

    ! output
    integer, intent(out) ::  &
         NumIncPTCSPrimCV (NbIncPTCSPrimMax, NbIncLocal),  &
         NumIncPTCSecondCV (NbEqFermetureMax, NbIncLocal)

    double precision, intent(out) :: &
         dXssurdXp (NbIncPTCSPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs (NbEqFermetureMax, NbIncLocal), &
         SmF (NbEqFermetureMax, NbIncLocal)
                      
    ! tmp
    integer :: k
    double precision :: &
         dFsurdX (NbIncPTCSMax, NbEqFermetureMax) ! (col,row) index order

    do k=1,NbIncLocal

         ! init tmp values for each cv
         call IncPrimSecd_init_cv(inc(k))

         ! compute dF/dX
         ! dFsurdX: (col, row) index order
         call IncPrimSecd_dFsurdX_cv(inc(k), rt(:,k), dFsurdX, SmF(:,k))

         ! choose inconnues prim and secd
         call IncPrimSecd_ps_cv(inc(k), dFsurdX, pschoice, &
            NumIncPTCSPrimCV(:,k), NumIncPTCSecondCV(:,k))

         ! compute dXssurdxp
         call IncPrimSecd_dXssurdXp_cv(inc(k), dFsurdX, SmF(:,k), &
            NumIncPTCSPrimCV(:,k), NumIncPTCSecondCV(:,k), &
            dXssurdXp(:,:,k), SmdXs(:,k))

    enddo

  end subroutine IncPrimSecd_compute_cv

  !> Update prim/secd arrays of nodes
  subroutine IncPrimSecd_compPrim_nodes

    integer :: k

    ! node
    call IncPrimSecd_compute_cv( &
         NbNodeLocal_Ncpus(commRank+1), &
         IncNode, NodeRocktypeLocal, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode,   &
                                !
         NumIncPTCSPrimNode, NumIncPTCSecondNode)

  end subroutine IncPrimSecd_compPrim_nodes


  subroutine IncPrimSecd_init_cv(inc)

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

  end subroutine IncPrimSecd_init_cv

  ! FIXME: reprendre notations de Xing et al. 2017
  ! F = closure equations
  !    *  molar fractions sum to 1
  !    * thermodynamic equilibrium between phases if any (fugacities equality)
  ! compute dFsurdX for each control volume
  subroutine IncPrimSecd_dFsurdX_cv(inc, rt, dFsurdX, SmF)

    type(TYPE_IncCVReservoir), intent(in) :: inc
    integer, intent(in) :: rt(IndThermique+1)
    double precision, intent(out) :: &  ! (col, row) index order
         dFsurdX(NbIncPTCSMax, NbEqFermetureMax)

    double precision, intent(out) :: &
         SmF(NbEqFermetureMax)

    integer :: i, mi, iph, iph1, iph2, icp, j, jph, numj, numc1, numc2
    double precision :: &
         f1, dPf1, dTf1, dCf1(NbComp), dSf1(NbPhase), &
         f2, dPf2, dTf2, dCf2(NbComp), dSf2(NbPhase)

    dFsurdX(:,:) = 0.d0
    SmF(:) = 0.d0

    ! --------------------------------------------------------------------------
    ! molar fractions sum to 1
    ! 1. F = sum_icp C_icp^iph(i) - 1, for i

    ! loop for rows associate with C_icp^iph(i) in dFsurdX

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
       SmF(i) = SmF(i) - 1.d0
    enddo

    mi = NbPhasePresente ! row is mi+i

    ! --------------------------------------------------------------------------
    ! thermodynamic equilibrium - fugacities equality
    ! 2. F = f_i^alpha * C_i^alpha - f_i^beta * C_i^beta
    do i=1, NbEqEquilibre !

       icp = NumCompEqEquilibre(i) ! component

       iph1 = Num2PhasesEqEquilibre(1,i) ! phase alpha
       iph2 = Num2PhasesEqEquilibre(2,i) ! phase beta

       numc1 = NumIncComp2NumIncPTC(icp,iph1) ! num of C_i^alpha in IncPTC
       numc2 = NumIncComp2NumIncPTC(icp,iph2) ! num of C_i^beta in IncPTC

       ! fugacity and div
       call f_Fugacity(rt, iph1, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph1), inc%Saturation, &
            f1, dPf1, dTf1, dCf1, dSf1)
       call f_Fugacity(rt, iph2, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph2), inc%Saturation, &
            f2, dPf2, dTf2, dCf2, dSf2)

       ! div Pression
       dFsurdX(1,i+mi) = dPf1*inc%Comp(icp,iph1) - dPf2*inc%Comp(icp,iph2)

#ifdef _THERMIQUE_
       ! div Temperature
       dFsurdX(2,i+mi) = dTf1*inc%Comp(icp,iph1) - dTf2*inc%Comp(icp,iph2)
#endif

       ! div Components
       ! d (f(P,T,C)*C_i)/dC_i = f + df/dC_i*C_i
       ! d (f(P,T,C)*C_i)/dC_j =     df/dC_j*C_i, j!=i
       dFsurdX(numc1, i+mi) = f1   ! first part of   d (f1(P,T,C)*C_i)/dC_i

       do j=1, NbComp     ! df1/dC_j*C_i    for every j which is in iph1
          if(MCP(j,iph1)==1) then ! phase iph1 contains component j

             numj = NumIncComp2NumIncPTC(j,iph1) ! num of C_j^iph1 in Inc
             dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
                  + inc%Comp(icp,iph1)*dCf1(j)
          end if
       end do

       dFsurdX(numc2,i+mi) = -f2   ! first part of   - d (f2(P,T,C)*C_i)/dC_i

       do j=1, NbComp     ! - df2/dC_j*C_i    for every j which is in iph2
          if(MCP(j,iph2)==1) then ! phase iph2 contains component j

             numj = NumIncComp2NumIncPTC(j,iph2) ! num of C_j^iph2 in Inc
             dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
                  - inc%Comp(icp,iph2)*dCf2(j)
          end if
       end do


       ! div Saturation
       do j=1, NbPhasePresente-1 ! row is i, col is j
         numj = j + NbIncPTC
         jph = NumPhasePresente(j)

         dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
           + dSf1(jph)*inc%Comp(icp,iph1) - dSf2(jph)*inc%Comp(icp,iph2)
       end do

       ! div Saturation secondaire
       jph = NumPhasePresente(NbPhasePresente)
       do j=1, NbPhasePresente-1 ! row is i, col is j
         numj = j + NbIncPTC

         dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
           - dSf1(jph)*inc%Comp(icp,iph1) + dSf2(jph)*inc%Comp(icp,iph2)
       end do

       ! SmF
       SmF(i+mi) = f1*inc%Comp(icp,iph1) - f2*inc%Comp(icp,iph2)
    end do

  end subroutine IncPrimSecd_dFsurdX_cv


  ! choose primary and secondary unknowns for each CV
  ! fill inc%Nb/NumIncPrim/Secd
  subroutine IncPrimSecd_ps_cv(inc, dFsurdX, pschoicecv, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    type(TYPE_IncCVReservoir), intent(in) :: inc
    ! dFsurdX may be used depending on choice method
    ! (e.g. projection on closure equations)
    ! ??? pourrait-on faire mieux qu'une extraction des composantes ???
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
       ! call IncPrimSecd_IncSecondGluton(inc, dFsurdX, &
       ! NumIncPTCSPrimCV, NumIncPTCSecondCV)
       call CommonMPI_abort("Glouton method is not implemented in IncPrimSecd")
    else if (pschoicecv==3) then ! Gauss method
       call IncPrimSecd_IncSecondGauss(inc, dFsurdX, &
            NumIncPTCSPrimCV, NumIncPTCSecondCV)
    end if

  end subroutine IncPrimSecd_ps_cv


  ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
  subroutine IncPrimSecd_dXssurdXp_cv(inc, dFsurdX, SmF, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, dXssurdXp, SmdXs)

    ! inputs
    type(TYPE_IncCVReservoir), intent(in) :: inc
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
    integer :: i, j, iph

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

       write(*,*)' inc ic ',inc%ic
       write(*,*)' inc P ',inc%Pression
       write(*,*)' inc T ',inc%Temperature
       write(*,*)' Sat ',inc%Saturation
       DO i=1,NbPhasePresente
         iph = NumPhasePresente(i)
         write(*,*)' phase  ', iph
         write(*,*)' C_i  ',inc%Comp(:,iph)
       ENDDO
       write(*,*)

    do i=1, NbEqFermeture ! = NbIncPTCSecond
       do j=1, NbEqFermeture
          write(*,*)' dFsurdXs ',i,j,dFsurdX(NumIncPTCSecondCV(j),i)
       enddo
       write(*,*)
    enddo

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

  end subroutine IncPrimSecd_dXssurdXp_cv

  ! Choix des inconnues primaires et secondaires
  !  a partir de la matrice dFsurdX par algorithme glouton
  !  minimisant les angles successifs
  subroutine IncPrimSecd_IncSecondGlouton(inc, dFsurdX, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc
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

  end subroutine IncPrimSecd_IncSecondGlouton


  !
  subroutine IncPrimSecd_IncSecondGauss(inc, dFsurdX, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc
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

  end subroutine IncPrimSecd_IncSecondGauss


  subroutine IncPrimSecd_PrimToSecd_C(increments_pointers) &
      bind(C, name="IncPrimSecd_PrimToSecd")

    type(Newton_increments_pointers), intent(in), value :: increments_pointers
    type(Newton_increments) :: increments

    call Newton_pointers_to_values(increments_pointers, increments)
    call IncPrimSecd_PrimToSecd( &
       increments%nodes, increments%fractures, increments%cells &
    )

  end subroutine IncPrimSecd_PrimToSecd_C

  ! compute prim values using secd values FIXME: l'inverse ?
  ! secd = SmdX - dXssurdXp * prim !! en fait delta_prim / élimination des secondaires
  ! v pour variation ?
  subroutine IncPrimSecd_PrimToSecd( &
       vnode, vfrac, vcell)

    real(c_double), dimension(:,:), intent(inout) :: &
         vnode, vfrac, vcell

    double precision :: &
         xp(NbCompThermique), &
         xs(NbEqFermetureMax)

    integer :: k, ic, i, iph
    integer :: &
         NbPhasePresente, NbEqFermeture, &
         NbIncPTC, NbIncPTCS, NbIncPTCPrim, NbIncPTCSPrim

    ! node
    call IncPrimSecd_PrimToSecd_cv( &
         IncNode, NbNodeLocal_Ncpus(commRank+1), &
         dXssurdXpNode, SmdXsNode, &
         NumIncPTCSPrimNode, NumIncPTCSecondNode, &
         vnode)

    ! frac
    call IncPrimSecd_PrimToSecd_cv( &
         IncFrac, NbFracLocal_Ncpus(commRank+1), &
         dXssurdXpFrac, SmdXsFrac, &
         NumIncPTCSPrimFrac, NumIncPTCSecondFrac, &
         vfrac)

    ! cell
    call IncPrimSecd_PrimToSecd_cv( &
         IncCell, NbCellLocal_Ncpus(commRank+1), &
         dXssurdXpCell, SmdXsCell, &
         NumIncPTCSPrimCell, NumIncPTCSecondCell, &
         vcell)

  end subroutine IncPrimSecd_PrimToSecd

  ! compute secd values using prim values for a set of cv 
  subroutine IncPrimSecd_PrimToSecd_cv( &
       inc, NbIncLocal, &
       dXssurdXp, SmdXs, &
       NumIncPTCSPrimCV, NumIncPTCSecondCV, &
       var_inc)

    ! input
    type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

    integer, intent(in) :: NbIncLocal

    double precision, intent(in) :: &
         dXssurdXp (NbIncPTCSPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs (NbEqFermetureMax, NbIncLocal)

    integer, intent(in) ::  &
         NumIncPTCSPrimCV (NbIncPTCSPrimMax, NbIncLocal),  &
         NumIncPTCSecondCV (NbEqFermetureMax, NbIncLocal)

    ! inout: need dvar_prim to fill dvar (both prim and secd)
    real(c_double), dimension(:,:), intent(inout) :: var_inc

    integer :: k, ic, i, iph
    integer :: &
         NbPhasePresente, NbEqFermeture, &
         NbIncPTC, NbIncPTCS, NbIncPTCPrim, NbIncPTCSPrim

    double precision :: &
         xp(NbCompThermique), &
         xs(NbEqFermetureMax)

    do k=1,NbIncLocal

         ic = inc(k)%ic
         NbPhasePresente = NbPhasePresente_ctx(ic)
         NbEqFermeture = NbEqFermeture_ctx(ic)
         NbIncPTC  = NbIncPTC_ctx(ic)
         NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente
         NbIncPTCSPrim = NbIncPTCSPrim_ctx(ic)
         NbIncPTCPrim = NbIncPTC - NbEqFermeture

         xp(1:NbCompThermique) = var_inc(1:NbCompThermique,k)
         xs(1:NbEqFermeture) = SmdXs(1:NbEqFermeture,k)

         ! http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
         ! y := alpha*A*x + beta*y
         ! Compute xs = xs - dXssurdXp * xp, where xs=SmdXs and xp=var_prim
         call dgemv('T', NbIncPTCSPrim, NbEqFermeture, &
         -1.d0, dXssurdXp(:,:,k), NbIncPTCSPrimMax, &
         xp(:), 1, -1.d0, xs(:), 1)

         !-----------------------------------------------------
         ! update var_inc with the primary/secondary increments
         !-----------------------------------------------------

         var_inc(:,k) = 0.d0

         ! copy prim P,T,C,S
         do i=1, NbIncPTCSPrim
            var_inc(NumIncPTCSPrimCV(i,k),k) = xp(i)
         end do

         ! copy secd P,T,C
         do i=1, NbEqFermeture
            var_inc(NumIncPTCSecondCV(i,k),k) = xs(i)
         end do

         ! fill secd S
         ! if NbPhasePresente=1,
         !    then this phase must be secd
         ! else last saturation is secd, the others are prim
         !    secd S = - sum_{S^alpha is prim} S^alpha
         if( NbPhasePresente==1) then

            ! iph is this phase in (P,T,C,S,n_i)
            iph = NbIncPTC + NumPhasePresente_ctx(1,ic)
            var_inc(iph,k) = 0.d0
         else

            ! iph is last present phase in vector (P,T,C,S,n_i)
            iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente,ic)
            var_inc(iph,k) = 0.d0
            do i=1, NbPhasePresente-1
               var_inc(iph,k) = var_inc(iph,k) - xp(NbIncPTCPrim+i)
            end do
         end if

         ! term prim n_k(X_j^n) ????
         ! copy prim n_i
         do i=1, NbCompCtilde_ctx(ic) ! =NbCompThermique-NbIncPTCSPrim
            var_inc(NbIncPTCS+i,k) = xp(NbIncPTCSPrim+i)
         end do

    end do ! loop over local inc

  end subroutine IncPrimSecd_PrimToSecd_cv 



  ! allocate
  subroutine IncPrimSecd_allocate

    integer :: nbCell, nbFrac, nbNode, nbNodeInj

    nbCell = NbCellLocal_Ncpus(commRank+1)
    nbFrac = NbFracLocal_Ncpus(commRank+1)
    nbNode = NbNodeLocal_Ncpus(commRank+1)
    nbNodeInj = NodeByWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)
    ! print*, 'IncPrimSecd_allocate', nbCell, nbFrac, nbNode, nbNodeInj

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

  end subroutine IncPrimSecd_allocate


  ! free
  subroutine IncPrimSecd_free

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

  end subroutine IncPrimSecd_free

end module IncPrimSecd
