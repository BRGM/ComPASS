!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecdFreeFlow

  ! FIXME: add only: check that all members are necessary
  use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort
  use DefModel, only: &
     NbComp, NbPhase, MCP, GAS_PHASE, &
     NbIncTotalPrimMax, NbEqFermetureMax, NbIncTotalMax, &
     pschoice, psprim, pssecd, &
     NumPhasePresente_ctx, NbPhasePresente_ctx, &
     NbIncPTCMax, IndThermique, NbCompThermique, NbEqEquilibreMax
  use Thermodynamics, only: f_Fugacity, f_PressionCapillaire
  use NumbyContext, only: &
     NbEqFermeture_ctx, NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, &
     NumIncComp2NumIncPTC_ctx, NbIncTotalPrim_ctx, NumCompCtilde_ctx, &
     NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
     NbEqEquilibre_ctx, NbIncPTC_ctx, NbIncTotal_ctx, NbCompCtilde_ctx
  use IncCVReservoir, only: &
     TYPE_IncCVReservoir, &
     IncCell, IncFrac, IncNode
  use IncPrimSecd, only: IncPrimSecd_ps_cv, &
      dXssurdXpCell, dXssurdXpFrac, dXssurdXpNode, &
      SmdXsCell, SmdXsFrac, SmdXsNode, &
      SmFCell, SmFFrac, SmFNode, &
      NumIncTotalPrimCell, NumIncTotalPrimFrac, NumIncTotalPrimNode, &
      NumIncTotalSecondCell, NumIncTotalSecondFrac, NumIncTotalSecondNode
  use Physics, only: atm_pressure
  use MeshSchema, only: &
     NbCellLocal_Ncpus, NbFracLocal_Ncpus, NbNodeLocal_Ncpus, &
     NodeRocktypeLocal, CellRocktypeLocal, FracRocktypeLocal, &
     NodeByWellInjLocal


  implicit none


   
  ! tmp values to simpfy notations of numerotation
  ! ex. NbPhasePresente = NbPhasePresente_ctx(inc%ic)
  integer, private :: &
       NbPhasePresente, NbCompCtilde, &
       NbEqFermeture, NbEqEquilibre,  & ! NbEqEquilibre -> Nombre d'Equation d'Equilibre thermodynamique fct du contexte (i.e. égalité des fugacités)
       NbIncPTC, NbIncTotal, &
       NbIncTotalPrim, &
                                !
       NumPhasePresente(NbPhase),               & ! Num -> identifiants de la phase présente
       NumCompCtilde(NbComp),                   & ! Num -> identifiants des composants absents
       NumCompEqEquilibre(NbEqEquilibreMax),    & ! identifiant des composant présents dans au moins 2 phases (donc concernés par égalité fugacités)
       Num2PhasesEqEquilibre(2, NbEqEquilibreMax), & ! phases impliquées dans l'équilibre (cf. NumCompEqEquilibre)
       NumIncComp2NumIncPTC(NbComp, NbPhase) ! matrice donnant pour chaque phase et chaque composant la ligne du "vecteur inconnu" 

  public :: &
       IncPrimSecdFreeFlow_compute, &
       IncPrimSecdFreeFlow_compPrim_nodes


  private :: &
       IncPrimSecdFreeFlow_compute_cv,             & ! all operations for one cv
       IncPrimSecdFreeFlow_init_cv,                & ! init infos according to ic (context) for each control volume (cv)
       IncPrimSecdFreeFlow_dFsurdX_cv,             & ! compute dF/dX for each cv
       IncPrimSecdFreeFlow_dXssurdXp_cv              ! compute dFs/dXp
                  
contains

  ! main subroutine of this module
  ! NumIncTotalPrim, NumIncTotalSecond filled in IncPrimSecd
  ! compute dXssurdXp, SmdXs
  ! of the FreeFlow d.o.f. (nodes only)
  subroutine IncPrimSecdFreeFlow_compute() &
        bind(C, name="IncPrimSecdFreeFlow_compute")

    ! FreeFlow BC nodes
    call IncPrimSecdFreeFlow_compute_cv( &
         NbNodeLocal_Ncpus(commRank+1), &
         IncNode, NodeRocktypeLocal, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode,   &
                                !
         NumIncTotalPrimNode, NumIncTotalSecondNode)


  end subroutine IncPrimSecdFreeFlow_compute


  ! all operations for a set of cv
  subroutine IncPrimSecdFreeFlow_compute_cv( &
       NbIncLocal, &
       inc, rt, &
       dXssurdXp, SmdXs, SmF, &
       NumIncTotalPrimCV, NumIncTotalSecondCV)

    ! input
    integer, intent(in) :: NbIncLocal

    type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

    integer, intent(in) :: rt(IndThermique+1, NbIncLocal)

    ! output
    integer, intent(out) ::  &
         NumIncTotalPrimCV (NbIncTotalPrimMax, NbIncLocal),  &
         NumIncTotalSecondCV (NbEqFermetureMax, NbIncLocal)

    double precision, intent(out) :: &
         dXssurdXp (NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs (NbEqFermetureMax, NbIncLocal), &
         SmF (NbEqFermetureMax, NbIncLocal)
                      
    ! tmp
    integer :: k
    double precision :: &
         dFsurdX (NbIncTotalMax, NbEqFermetureMax) ! (col,row) index order

    do k=1,NbIncLocal

         ! done only for ff dof
         if(inc(k)%ic<2**NbPhase) cycle ! FIXME: loop over freeflow dof only, avoid reservoir node

         ! init tmp values for each cv
         call IncPrimSecdFreeFlow_init_cv(inc(k))

         ! compute dF/dX
         ! dFsurdX: (col, row) index order
         call IncPrimSecdFreeFlow_dFsurdX_cv(inc(k), rt(:,k), dFsurdX, SmF(:,k))

         ! FIXME: si je peux faire IncPrimSecd_compute sur les dof FF, enlever ce call et remettre protected à NumIncTotal... 
         ! (mais cela m'étonnerait car le numb d'eq de fermeture ne correspond pas au nb d'inconnus secds dans ce cas)
         ! choose prim and sced unknowns
         call IncPrimSecd_ps_cv(inc(k), dFsurdX, pschoice, &
            NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k))

         ! compute dXssurdxp
         call IncPrimSecdFreeFlow_dXssurdXp_cv(inc(k), dFsurdX, SmF(:,k), &
            NumIncTotalPrimCV(:,k), NumIncTotalSecondCV(:,k), &
            dXssurdXp(:,:,k), SmdXs(:,k))

    enddo

  end subroutine IncPrimSecdFreeFlow_compute_cv

  !> Update prim/secd arrays of nodes (same as IncPrimSecdFreeFlow_compute)
  ! FIXME NumIncTotalPrimNode, NumIncTotalSecondNode are updated in IncPrimSecd.F90, is it necessary to update dXssurdXpNode, SmdXsNode ?
  subroutine IncPrimSecdFreeFlow_compPrim_nodes

    ! node
    call IncPrimSecdFreeFlow_compute_cv( &
         NbNodeLocal_Ncpus(commRank+1), &
         IncNode, NodeRocktypeLocal, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode,   &
                                !
         NumIncTotalPrimNode, NumIncTotalSecondNode)

  end subroutine IncPrimSecdFreeFlow_compPrim_nodes


  subroutine IncPrimSecdFreeFlow_init_cv(inc)

    type(TYPE_IncCVReservoir), intent(in) :: inc

    NbPhasePresente = NbPhasePresente_ctx(inc%ic)
    NbCompCtilde = NbCompCtilde_ctx(inc%ic)

    NbEqFermeture = NbEqFermeture_ctx(inc%ic)
    NbEqEquilibre = NbEqEquilibre_ctx(inc%ic)

    NbIncPTC = NbIncPTC_ctx(inc%ic)
    NbIncTotal = NbIncTotal_ctx(inc%ic)

    NbIncTotalPrim = NbIncTotalPrim_ctx(inc%ic)

    NumPhasePresente(:) = NumPhasePresente_ctx(:,inc%ic)
    NumCompCtilde(:) = NumCompCtilde_ctx(:,inc%ic)
    NumCompEqEquilibre(:) = NumCompEqEquilibre_ctx(:,inc%ic)

    Num2PhasesEqEquilibre(:,:) = Num2PhasesEqEquilibre_ctx(:,:,inc%ic)
    NumIncComp2NumIncPTC(:,:) = NumIncComp2NumIncPTC_ctx(:,:,inc%ic)

  end subroutine IncPrimSecdFreeFlow_init_cv

  ! FIXME: reprendre notations de Xing et al. 2017
  ! F = closure equations
  !    * molar fractions sum to 1 for present phases
  !    * thermodynamic equilibrium between phases if any (fugacities equality)
  !    * P^g = P^atm
  ! Remark: Sg+Sl=1 is not a closure law, as well as Pphase=Pref+Pc(phase) and q^l_atm * S^g = 0, it is forced in the implementation
  ! compute dFsurdX for each control volume (Careful: the lines of the derivatives must coincide with the index of unknowns in DefModel.F90)
  !      dFsurdX(1,:)                                            derivative Pressure
  !      #ifdef _THERMIQUE_ dFsurdX(2,:)                         derivative Temperature
  !      dFsurdX(2+IndThermique:NbEquilibre+IndThermique+1,:)    derivative Components
  !      dFsurdX(NbIncPTC+1:NbIncPTC+NbPhasePresente+1, :)       derivative principal Saturations
  !      dFsurdX(, :)                                            derivative freeflow flowrate
  subroutine IncPrimSecdFreeFlow_dFsurdX_cv(inc, rt, dFsurdX, SmF)

    type(TYPE_IncCVReservoir), intent(in) :: inc
    integer, intent(in) :: rt(IndThermique+1)
    double precision, intent(out) :: &  ! (col, row) index order
         dFsurdX(NbIncTotalMax, NbEqFermetureMax)

    double precision, intent(out) :: &
         SmF(NbEqFermetureMax)

    integer :: i, mi, iph, iph1, iph2, icp, j, jph, jph_scd, numj, numc1, numc2
    double precision :: &
         f1, dPf1, dTf1, dCf1(NbComp), dSf1(NbPhase), &
         f2, dPf2, dTf2, dCf2(NbComp), dSf2(NbPhase), &
         Pc, dSPc(NbPhase)

    dFsurdX(:,:) = 0.d0  ! local to this file, cannot rely on previous computation in IncPrimSecd
    SmF(:) = 0.d0

    ! --------------------------------------------------------------------------
    ! molar fractions sum to 1
    ! 1. F = sum_icp C_icp^iph(i) - 1, for i

    ! loop for rows associate with C_icp^iph(i) in dFsurdX
    do i=1, NbPhasePresente  ! row is i, col is j
       iph = NumPhasePresente(i)

       do icp=1, NbComp ! loop for cols
          if(MCP(icp,iph)==1) then
             j = NumIncComp2NumIncPTC(icp,iph)  ! num of C_icp^iph in IncPTC
             dFsurdX(j,i) = 1.d0    ! derivative wrt C

             SmF(i) = SmF(i) + inc%Comp(icp,iph)
          end if
       end do
       SmF(i) = SmF(i) - 1.d0
    enddo

    mi = NbPhasePresente ! nb of closure eq already stored, mi+1 futur row of dFsurdX and SmF

    ! --------------------------------------------------------------------------
    ! thermodynamic equilibrium - fugacities equality
    ! 2. F = f_i^alpha * C_i^alpha - f_i^beta * C_i^beta
    do i=1, NbEqEquilibre !

       icp = NumCompEqEquilibre(i) ! component i

       iph1 = Num2PhasesEqEquilibre(1,i) ! phase alpha
       iph2 = Num2PhasesEqEquilibre(2,i) ! phase beta

       numc1 = NumIncComp2NumIncPTC(icp,iph1) ! num of C_i^alpha in IncPTC
       numc2 = NumIncComp2NumIncPTC(icp,iph2) ! num of C_i^beta in IncPTC

       ! fugacity
       call f_Fugacity(rt, iph1, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph1), inc%Saturation, &
            f1, dPf1, dTf1, dCf1, dSf1)
       call f_Fugacity(rt, iph2, icp, inc%Pression, inc%Temperature, &
            inc%Comp(:,iph2), inc%Saturation, &
            f2, dPf2, dTf2, dCf2, dSf2)

       ! derivative Pressure
       dFsurdX(1,i+mi) = dPf1*inc%Comp(icp,iph1) - dPf2*inc%Comp(icp,iph2)

#ifdef _THERMIQUE_
       ! derivative Temperature
       dFsurdX(2,i+mi) = dTf1*inc%Comp(icp,iph1) - dTf2*inc%Comp(icp,iph2)
#endif

       ! derivative Components
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

       ! derivative primary Saturations
       ! with contribution of secondary Saturation
       ! because sum(saturations)=1 is eliminated
       jph_scd = NumPhasePresente(NbPhasePresente)
       do j=1, NbPhasePresente-1
         numj = j + NbIncPTC
         jph = NumPhasePresente(j)

         dFsurdX(numj,i+mi) = dFsurdX(numj,i+mi) &
           + dSf1(jph)*inc%Comp(icp,iph1) - dSf2(jph)*inc%Comp(icp,iph2) &
           - dSf1(jph_scd)*inc%Comp(icp,iph1) + dSf2(jph_scd)*inc%Comp(icp,iph2)
       end do

       ! SmF
       SmF(i+mi) = f1*inc%Comp(icp,iph1) - f2*inc%Comp(icp,iph2)
    end do

    mi = mi + NbEqEquilibre ! mi+1 futur row of dFsurdX and SmF

    ! --------------------------------------------------------------------------
    ! 3. P^g - P^atm = 0     ie      Pref + Pc(GAS_PHASE) - P^atm = 0 
    call f_PressionCapillaire(rt,GAS_PHASE,inc%Saturation,Pc,DSPc)

    ! derivative Pressure
    dFsurdX(1,mi+1) = 1.d0

    ! derivative primary Saturations
    ! with contribution of secondary Saturation
    ! because sum(saturations)=1 is eliminated
    jph_scd = NumPhasePresente(NbPhasePresente)
    do j=1, NbPhasePresente-1
      numj = j + NbIncPTC
      jph = NumPhasePresente(j)

      dFsurdX(numj,mi+1) = DSPc(jph) - DSPc(jph_scd)
    end do

    ! SmF
    SmF(mi+1) = inc%Pression + Pc - atm_pressure

    
   end subroutine IncPrimSecdFreeFlow_dFsurdX_cv

  ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
  subroutine IncPrimSecdFreeFlow_dXssurdXp_cv(inc, dFsurdX, SmF, &
       NumIncTotalPrimCV, NumIncTotalSecondCV, dXssurdXp, SmdXs)

    ! inputs
    type(TYPE_IncCVReservoir), intent(in) :: inc
    double precision, intent(in) ::  & ! (col, row) index order
         dFsurdX(NbIncTotalMax, NbEqFermetureMax), &
         SmF(NbEqFermetureMax)

    integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

    ! output
    double precision, intent(out) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

    ! tmp
    double precision :: & ! (row,col) index order, lapack
         dFsurdX_prim(NbEqFermetureMax, NbIncTotalPrimMax), &
         dFsurdX_secd(NbEqFermetureMax, NbEqFermetureMax)

    ! parameters for lapack
    integer :: ipiv(NbEqFermetureMax), info, Ierr, errcode
    integer :: i, j, iph

    ! from dFsurdX, take out the cols of prim and secd variable
    do j=1, NbIncTotalPrim
       do i=1, NbEqFermeture
          dFsurdX_prim(i,j) = dFsurdX(NumIncTotalPrimCV(j),i)
       enddo
    enddo

    do j=1, NbEqFermeture ! = Nb of secd unknowns
       do i=1, NbEqFermeture
          dFsurdX_secd(i,j) = dFsurdX(NumIncTotalSecondCV(j),i)
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

    do i=1, NbEqFermeture ! = Nb of secd unknowns
       do j=1, NbEqFermeture
          write(*,*)' dFsurdXs ',i,j,dFsurdX(NumIncTotalSecondCV(j),i)
       enddo
       write(*,*)
    enddo

       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if ! info

    call dgetrs('N', NbEqFermeture, NbIncTotalPrim, &
         dFsurdX_secd, NbEqFermetureMax, &
         ipiv, dFsurdX_prim, NbEqFermetureMax, info)
    if(info /=0) then
       write(0,'(A,I0)') "dgetrs error in dXssurdxp, info = ", info
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    do j=1, NbIncTotalPrim
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

  end subroutine IncPrimSecdFreeFlow_dXssurdXp_cv

end module IncPrimSecdFreeFlow
