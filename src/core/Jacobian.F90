!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Jacobian

  ! workflow:
  !   1. sturctures of Jacobian before and after Schur
  !   2. Jacobian -> Regularization
  !      -> Alignment -> Schur

  use iso_c_binding, only: c_double
  use mpi, only: MPI_Abort
  use CommonType, only: CSRArray2dble
  use CommonMPI, only: &
     commRank, ComPASS_COMM_WORLD, CommonMPI_abort

  use DefModel, only: &
     NumPhasePresente_ctx, NbPhasePresente_ctx, &
     NbComp, NbPhase, NbCompThermique, MCP, aligmat, aligmethod

  use LoisThermoHydro, only: &
     DensiteMolaireKrViscoCompWellInj, DensiteMolaireKrViscoEnthalpieWellInj, &
     DensiteMolaireKrViscoCompNode, DensiteMolaireKrViscoCompCell, DensiteMolaireKrViscoCompFrac, &
     DensiteMolaireKrViscoEnthalpieNode, DensiteMolaireKrViscoEnthalpieCell, DensiteMolaireKrViscoEnthalpieFrac, &
     DensiteMolaireSatComp, &
     DensiteMolaireEnergieInterneSat, &
     divDensiteMolaireKrViscoCompNode, divDensiteMolaireKrViscoCompCell, divDensiteMolaireKrViscoCompFrac, &
     divDensiteMolaireKrViscoEnthalpieNode, divDensiteMolaireKrViscoEnthalpieCell, divDensiteMolaireKrViscoEnthalpieFrac, &
     divTemperatureNode, divTemperatureCell, divTemperatureFrac, &
     SmTemperatureNode, SmTemperatureCell, SmTemperatureFrac, &
     SmDensiteMassiqueNode, SmDensiteMassiqueCell, SmDensiteMassiqueFrac, &
     SmPressionNode, SmPressionCell, SmPressionFrac, &
#ifdef _WIP_FREEFLOW_STRUCTURES_
     divFreeFlowMolarFlowrateNode, SmFreeFlowMolarFlowrateNode, &
     FreeFlowMolarFlowrateCompNode, divFreeFlowMolarFlowrateCompNode, SmFreeFlowMolarFlowrateCompNode, &
     FreeFlowHmCompNode, divFreeFlowHmCompNode, SmFreeFlowHmCompNode, &
     divFreeFlowHTTemperatureNetRadiationNode, SmFreeFlowHTTemperatureNetRadiationNode, &
     FreeFlowMolarFlowrateEnthalpieNode, divFreeFlowMolarFlowrateEnthalpieNode, SmFreeFlowMolarFlowrateEnthalpieNode, &
     AtmEnthalpieNode, &
#endif
     SmDensiteMolaireKrViscoCompNode, SmDensiteMolaireKrViscoCompCell, SmDensiteMolaireKrViscoCompFrac, &
     SmDensiteMolaireKrViscoEnthalpieNode, SmDensiteMolaireKrViscoEnthalpieCell, SmDensiteMolaireKrViscoEnthalpieFrac, &
     SmDensiteMolaireEnergieInterneSatNode, SmDensiteMolaireEnergieInterneSatCell, SmDensiteMolaireEnergieInterneSatFrac, &
     divDensiteMassiqueNode, divDensiteMassiqueCell, divDensiteMassiqueFrac, &
     divPressionCapNode, divPressionCapCell, divPressionCapFrac, &
     divDensiteMolaireKrViscoEnthalpieWellInj, divDensiteMolaireKrViscoCompWellInj, &
     divDensiteMolaireSatCompNode, divDensiteMolaireSatCompCell, divDensiteMolaireSatCompFrac, &
     divDensiteMolaireEnergieInterneSatNode, divDensiteMolaireEnergieInterneSatCell, divDensiteMolaireEnergieInterneSatFrac, &
     divPressionNode, divPressionCell, divPressionFrac, &
     divPressionCapNode, divPressionCapCell, divPressionCapFrac, &
     SmDensiteMolaireSatComp

   use NumbyContext, only: &
      NbCompCtilde_ctx, NumCompCtilde_ctx, NbIncPTC_ctx

  use Physics, only: gravity, CpRoche, atm_comp

  use Newton, only: Newton_increments_pointers, Newton_increments, Newton_pointers_to_values
  use SchemeParameters, only: eps

  use IncCVReservoir, only: &
     IncNode, IncCell, IncFrac
  use IncCVWells, only: &
     PerfoWellInj, DataWellInjLocal, &
     PerfoWellProd, PerfoWellInj
  use VAGFrac, only: &
     TkLocal_Darcy, TkLocal_Fourier, TkFracLocal_Darcy, TkFracLocal_Fourier, &
     VolDarcy, &
     PoroVolDarcy, &
     PoroVolFourier, &
     Poro_1VolFourier

  use IncPrimSecd, only: &
     NbIncTotalPrim_ctx, &
     NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
     NumIncTotalPrimNode, NumIncTotalPrimCell, NumIncTotalPrimFrac, &
     NbIncTotalPrimMax

  use MeshSchema, only: &
     IdNodeLocal, &
#ifdef _WIP_FREEFLOW_STRUCTURES_
     IdFFNodeLocal, &
#endif
     NodebyCellLocal, FracbyCellLocal, NodebyFaceLocal, FaceToFracLocal, &
     NodebyWellProdLocal, NodeDatabyWellProdLocal, &
     NodebyWellInjLocal, NodeDatabyWellInjLocal, &
     NbCellLocal_Ncpus, NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbNodeOwn_Ncpus, &
     NbFracOwn_Ncpus, NbWellProdLocal_Ncpus, NbFaceOwn_Ncpus, &
     NbWellInjLocal_Ncpus, NbWellInjOwn_Ncpus, NbWellProdOwn_Ncpus, DataWellProdLocal, &
     CellbyNodeOwn, NodebyFracOwn, NodebyNodeOwn, FracbyNodeOwn, FracbyFracOwn, &
     WellInjbyNodeOwn, WellProdbyNodeOwn, CellbyFracOwn, &
     NbFracCellMax, NbNodeCellMax, NbNodeFaceMax, &
     FracToFaceLocal, XNodeLocal, XCellLocal, XFaceLocal, &
     SurfFreeFlowLocal

  use Flux, only: &
     FluxDarcyKI, FluxDarcyFI
  use Residu, only: &
     ResiduNode, ResiduCell, ResiduFrac, &
     ResiduWellInj, ResiduWellProd

  implicit none

  !            node, frac, cell, wellinj, wellprod
  !           | a11, a12, a13, a14, a15 | node own
  ! JacBigA = | a21, a22, a23, 0,   0   | frac own
  !           | a31, a32, a33, 0,   0   | cell (own and ghost)
  !           ! a41, 0,   0,   a44, 0   | wellinj own
  !           ! a51, 0,   0,   0,   a55 | wellprod own

  !            node, frac, wellinj, wellprod
  ! JacA =    | A11, A12, A13, A14 | node own
  !           | A21, A22, 0,   0   | frac own
  !           | A31, 0,   A33, 0   | wellinj own,
  !           | A41, 0,   0,   A44 | wellprod own
  type(CSRArray2dble), public :: JacBigA
  type(CSRArray2dble), public :: JacA

  ! second membre before and after Schur complement
  double precision, allocatable, dimension(:,:), public :: &
       bigSm, Sm

  ! ipiv(:,k): pivot indices for cell k
  integer, allocatable, dimension(:,:), private :: ipiv

  ! JacBigA is a sparse matrix.
  ! for an element in JacBigA, we need to known its num in JacBigA%Val.
  ! the following vectors are used for this purpose
  integer, allocatable, dimension(:), private :: csrK , csrSR

  public :: &
       Jacobian_StrucJacBigA,  & !< non-zero structure of Jacobian before Schur
       Jacobian_StrucJacA,     & !< non-zero structure of Jacobian after Schur
       Jacobian_JacBigA_BigSm, & !< Jacobian and second member
       Jacobian_Schur,       & !< Schur complement
       Jacobian_free

  private :: &
       ! the fill of Jacobian/Second membre is decomposed into three subroutines
       !   (1) Jacobian_JacBigA_BigSm_accmolaire: term n_k(X_j^n)
       !   (2) Jacobian_JacBigA_BigSm_cell: loop of cell
       !      (2.1) nodes in cell
       !      (2.2) fracs in cell
       !   (3) Jacobian_JacBigA_BigSm_frac: loop of frac
       Jacobian_JacBigA_BigSm_accmolaire, &  ! (1)
       Jacobian_JacBigA_BigSm_cell,       &  ! (2)
       Jacobian_JacBigA_BigSm_frac,       &  ! (3)
       Jacobian_JacBigA_BigSm_wellinj,    &  !
       Jacobian_JacBigA_BigSm_wellprod,   &  !
       !
       ! In (2), we compute div( DensiteMolaire*Kr*Visco*Comp)*V_{k,s}
       !                  + DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
       !                  + div( DensiteMolaire*Kr*Visco*Enthalpie*FluxFourier ) if thermique
       !         V_{k,s} is Darcy Flux: k is cell, s is nodes/fracs; (2.1)
       !                                k is frac, s is nodes (2.2)
       !
       ! To compute DensiteMolaire*Kr*Visco*Comp*div(V_{k,s}),
       !    first div(pho) in the loop of k,s
       !    then  div(V_{k,s}) in the loop of k,s
       !    last  DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
       !
       ! div(DensiteMolaire*Kr*Visco*Comp)*V_{k,s}
       Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellnode, &  ! k is cell, s is node
       Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellfrac, &  ! k is cell, s is frac
       Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_fracnode, &  ! k is frac, s is node

       ! DensiteMolaire*Kr*Visco*Comp*div(V_{k,s})
       Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellnode,     &  ! k is cell, s is node
       Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellfrac,     &  ! k is cell, s is frac
       Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_fracnode,     &  ! k is frac, s is node

       ! div(V_{k,s})
       Jacobian_divDarcyFlux_cellnode,   & ! k is cell, s is node
       Jacobian_divDarcyFlux_cellfrac,   & ! k is cell, s is frac
       Jacobian_divDarcyFlux_fracnode,   & ! k is frac, s is node

       ! div(rho)
       Jacobian_divrho_cellnode, &
       Jacobian_divrho_cellfrac, &
       Jacobian_divrho_fracnode, &

       ! For thermique, we compute div(FluxFourier)
       ! then div( DensiteMolaire*Kr*Visco*Enthalpie*FluxFourier )
       Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellnode, &
       Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellfrac, &
       Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_fracnode, &

       ! div(FourierFlux)
       Jacobian_divFourierFlux_cellnode, &
       Jacobian_divFourierFlux_cellfrac, &
       Jacobian_divFourierFlux_fracnode, &

#ifdef _WIP_FREEFLOW_STRUCTURES_
       ! div(FreeFlow)
       Jacobian_JacBigA_BigSm_FF_node, &
       Jacobian_divMolarFreeFlow_node, &  ! k is cell, s is node
       Jacobian_divThermalFreeFlow_node, &
#endif

       Jacobian_RowCol_KSR, &
       Jacobian_RowCol_FR,  &

       ! regularization/alignment
       Jacobian_Regularization,     &
       Jacobian_Regularization_row, &
       Jacobian_Alignment_diag,     &
       Jacobian_Alignment_diag_row, &
       Jacobian_Alignment_man,      &
       Jacobian_Alignment_man_row,  &
      
       Jacobian_JacBigA_locate_frac_row

contains

    subroutine dump_jacobian(specific_row, specific_col)

    integer, optional, intent(in) :: specific_row, specific_col
    integer :: i, j, s, n
    double precision :: a

    do s=1, NbNodeOwn_Ncpus(commRank+1)
        if( .not.present(specific_row) .or. s==specific_row) then
            do n=JacBigA%Pt(s)+1, JacBigA%Pt(s+1)
                if( .not.present(specific_col) .or. JacBigA%Num(n)==specific_col) then
                    write(*,*) 'nodes', s, JacBigA%Num(n), 'nonzero', n, 'on proc', commRank
                    do j=1, NbCompThermique
                        do i=1, NbCompThermique
                            a = JacBigA%Val(i,j,n)
                            if(a==0. .or. a==1.) then
                                write(*,"(I15)",advance="no") int(a)
                            else
                                write(*,"(E15.7)",advance="no") a
                            end if
                        end do
                        write(*,*)
                    end do
                end if
            end do
        end if
    end do

    end subroutine dump_jacobian

  !> \brief Compute Jacobian and rhs
  subroutine Jacobian_ComputeJacSm(Delta_t)  &
        bind(C, name="Jacobian_ComputeJacSm")

    real(c_double), intent(in), value :: Delta_t
    integer :: errcode, Ierr

    ! Jacobian and second member
    call Jacobian_JacBigA_BigSm(Delta_t)

    ! Regularization of Jacobian
    call Jacobian_Regularization

    ! Schur complement of Jacobian and second member
    call Jacobian_Schur

    ! Alignment of Jacobian after Schur
    if (aligmethod==1) then ! manually

       call Jacobian_Alignment_man

    else if (aligmethod==2) then ! diag inverse

       call Jacobian_Alignment_diag

    else
       write(0,'(A)') "Error: Alignment method is not defined"
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

  end subroutine Jacobian_ComputeJacSm

  !> \brief Fill BigSm with the residual values     
  !!
  !!   If the node is dir, we suppose that the residu is 0,
  !!     in other words, we set its second member as 0.
  !!
  !!   The Jacobian corresponding to dir nodes are Id matrix,
  !!     setting bigSm(dirichlet)=0 or not has no influence mathematically,
  !!     however, it could have influence for linear solver
  !!     since the values could be very large 10**7.
  !!     It is observed when debugging!!!
  subroutine Jacobian_JacBigA_BigSm_init_from_residual()

    integer :: j, start


    do j=1, NbNodeOwn_Ncpus(commRank+1)
       if(IdNodeLocal(j)%P=="d") then
          bigSm(:,j) = 0.d0
       else
          bigSm(:,j) = - ResiduNode(:,j)
       end if
    end do

    start = NbNodeOwn_Ncpus(commRank+1)
    do j=1, NbFracOwn_Ncpus(commRank+1)
       bigSm(:,j+start) = - ResiduFrac(:,j)
    end do

    start = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)
    do j=1, NbCellLocal_Ncpus(commRank+1)
       bigSm(:,j+start) = - ResiduCell(:,j)
    end do

    start = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1) &
         + NbCellLocal_Ncpus(commRank+1)
    do j=1, NbWellInjOwn_Ncpus(commRank+1)
       bigSm(1,j+start) = - ResiduWellInj(j)   ! Pressure
       bigSm(2:NbCompThermique,j+start) = 0.d0 ! not used
    end do

    start = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1) &
         + NbCellLocal_Ncpus(commRank+1) + NbWellInjOwn_Ncpus(commRank+1)
    do j=1, NbWellProdOwn_Ncpus(commRank+1)
       bigSm(1,j+start) = - ResiduWellProd(j)  ! Pressure
       bigSm(2:NbCompThermique,j+start) = 0.d0 ! not used
    end do

   end subroutine Jacobian_JacBigA_BigSm_init_from_residual

  !> \brief fill Jacobian and second member before Schur: main subroutine
  subroutine Jacobian_JacBigA_BigSm(Delta_t)

    double precision, intent(in) :: Delta_t

    integer :: j, start, s, i, nz

    !> 1. init second membre
    call Jacobian_JacBigA_BigSm_init_from_residual
    
    JacBigA%Val(:,:,:) = 0.d0

    !> 2. div prim and Sm for term n_k(X_j^n)
    call Jacobian_JacBigA_BigSm_accmolaire(Delta_t)

    !> 3.1 loop of cell
    call Jacobian_JacBigA_BigSm_cell

    !> 3.2 loop of frac
    call Jacobian_JacBigA_BigSm_frac

    !> 3.3 loop of well inj
    call Jacobian_JacBigA_BigSm_wellinj

    !> 3.4 loop of well prod
    call Jacobian_JacBigA_BigSm_wellprod

#ifdef _WIP_FREEFLOW_STRUCTURES_
    !> 3.5 loop of FreeFlow Nodes
    call Jacobian_JacBigA_BigSm_FF_node
#endif

    !> 4. Dirichlet in Jacobian
    do s=1, NbNodeOwn_Ncpus(commRank+1)

       ! Darcy dir
       if(IdNodeLocal(s)%P=="d") then

          ! the diagonal element (s,s) is JacBigA(nz)
          do i=JacBigA%Pt(s)+1, JacBigA%Pt(s+1)
             if( JacBigA%Num(i)==s) then
                nz = i
                exit
             end if
          end do

          ! Identity matrix for Darcy
          ! JacBigA%Val(:,:,nz) = 0.d0
          do j=1, NbComp
             JacBigA%Val(j,j,nz) = 1.d0
          end do
       end if

#ifdef _THERMIQUE_

       ! Fourier dir
       if(IdNodeLocal(s)%T=="d") then

          ! the diagonal element (s,s) is JacBigA(nz)
          do i=JacBigA%Pt(s)+1, JacBigA%Pt(s+1)
             if( JacBigA%Num(i)==s) then
                nz = i
                exit
             end if
          end do

          ! Identity for Fourier
          JacBigA%Val(NbComp+1,NbComp+1,nz) = 1.d0
       end if
#endif

    end do

    ! if(commRank==1) then

    !    print*, ""
    !    open(unit=11, file='res.txt', status='unknown')

    !    print*, 5888.75d0-XCellLocal(1,15813), &
    !         XCellLocal(2,15813), XCellLocal(3,15813)

    !    do i=15813+NbNodeOwn_Ncpus(commRank+1), &
    !         15813+NbNodeOwn_Ncpus(commRank+1)

    !       write(*,'(ES22.13)') bigSm(:,i)
    !       write(11,'(ES22.13)') bigSm(:,i)
    !       print*, ""
    !    end do

    !    close(11)
    ! end if

    ! if(commRank==0) then

    !    print*, ""

    !    open(unit=11, file='res.txt', status='unknown')

    !    do s=1, JacBigA%Nb
    !       do nz=JacBigA%Pt(s)+1, JacBigA%Pt(s+1)

    !          if(s==8 .and. JacBigA%Num(nz)==19) then
    !             do i=1, NbCompThermique
    !                do j=1, NbIncTotalPrim_ctx(1)
    !                   write(*,'(ES22.13)') JacBigA%Val(j,i,nz)
    !                   write(11,'(ES22.13)') JacBigA%Val(j,i,nz)
    !                end do
    !                print*, ""
    !                write(11,*) ""
    !             end do
    !          end if

    !       end do
    !    end do

    !    close(11)
    ! end if

  end subroutine Jacobian_JacBigA_BigSm


  !> \brief Sub-subroutine of Jacobian_JacBigA_BigSm for term n_k(X_j^n)
  subroutine Jacobian_JacBigA_BigSm_accmolaire(Delta_t)

    double precision, intent(in) :: Delta_t

    integer :: k, rowk, i, icp, nz, j, m, mph

    ! 2.1 div prim: n_k(X_j^n), k is node
    do k=1, NbNodeOwn_Ncpus(commRank+1) ! node

       ! look for diagonal
#ifdef _THERMIQUE_
       if ((IdNodeLocal(k)%P /= "d") .or. (IdNodeLocal(k)%T /= "d")) then
#else
       if ((IdNodeLocal(k)%P /= "d")) then
#endif
          rowk = k ! row of node k in vector (node own, frac own, cell)

          ! the diagonal element (k,k) is JacBigA(nz)
          do i=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
             if( JacBigA%Num(i)==k) then
                nz = i
                exit
             end if
          end do
       end if

       if (IdNodeLocal(k)%P /= "d") then

          ! loop of components in Ctilde, n_k is an unknown independent
          do i=1, NbCompCtilde_ctx(IncNode(k)%ic)
             icp = NumCompCtilde_ctx(i, IncNode(k)%ic)

             j = NbIncTotalPrim_ctx(IncNode(k)%ic) + i
             JacBigA%val(j,icp,nz) = VolDarcy%nodes(k) / Delta_t
          end do

          ! loop of component, n_k is not an unknown independent
          do m=1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k, k is node
             mph = NumPhasePresente_ctx(m,IncNode(k)%ic)

             do icp=1, NbComp
                if( MCP(icp,mph)==1) then ! Q_k \cap P_i

                   do j=1, NbIncTotalPrim_ctx(IncNode(k)%ic)
                      JacBigA%Val(j,icp,nz) = JacBigA%Val(j,icp,nz) &
                           + divDensiteMolaireSatCompNode(j,icp,m,k) * PoroVolDarcy%nodes(k) / Delta_t
                   end do

                   bigSm(icp,rowk) = bigSm(icp,rowk) &
                        - SmDensiteMolaireSatComp%nodes(icp,m,k) * PoroVolDarcy%nodes(k) / Delta_t

                end if
             end do

          end do ! end of loop n_k
       end if

#ifdef _THERMIQUE_

       if (IdNodeLocal(k)%T /= "d") then

          do m=1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k, k is node
             mph = NumPhasePresente_ctx(m,IncNode(k)%ic)

             do j=1, NbIncTotalPrim_ctx(IncNode(k)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + PoroVolFourier%nodes(k) * divDensiteMolaireEnergieInterneSatNode(j,m,k) / Delta_t
             end do

             bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
                  - PoroVolFourier%nodes(k) * SmDensiteMolaireEnergieInterneSatNode(m,k) / Delta_t
          end do

          do j=1, NbIncTotalPrim_ctx(IncNode(k)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                  + Poro_1VolFourier%nodes(k) * CpRoche * divTemperatureNode(j,k) / Delta_t
          end do

          bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
               - Poro_1VolFourier%nodes(k) * CpRoche * SmTemperatureNode(k) / Delta_t
       end if

#endif
    end do ! end of div prim n_k, node


    ! 2.2 div prim: n_k(X_j^n), k is frac
    do k=1, NbFracOwn_Ncpus(commRank+1) ! node

       rowk = k + NbNodeOwn_Ncpus(commRank+1) ! row of frac k in vector (node own, frac own, cell)

       ! the diagonal element (k,k) is KacBigA(nz)
       do i=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
          if( (JacBigA%Num(i)-NbNodeLocal_Ncpus(commRank+1))==k) then
             nz = i
             exit
          end if
       end do

       ! loop of composants in Ctilde, n_k is an unknown independent
       do i=1, NbCompCtilde_ctx(IncFrac(k)%ic)
          icp = NumCompCtilde_ctx(i, IncFrac(k)%ic)

          j = NbIncTotalPrim_ctx(IncFrac(k)%ic) + i
          JacBigA%Val(j,icp,nz) = VolDarcy%fractures(k) / Delta_t
       end do

       ! loop of composants, n_k is not an unknown independent
       do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is node
          mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

          do icp=1, NbComp
             if( MCP(icp,mph)==1) then ! Q_k \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                   JacBigA%Val(j,icp,nz) = JacBigA%Val(j,icp,nz) &
                        + divDensiteMolaireSatCompFrac(j,icp,m,k) * PoroVolDarcy%fractures(k) / Delta_t
                end do

                bigSm(icp,rowk) = bigSm(icp,rowk) &
                     - SmDensiteMolaireSatComp%fractures(icp,m,k) * PoroVolDarcy%fractures(k) / Delta_t
             end if
          end do

       end do ! end of loop phase

#ifdef _THERMIQUE_

       do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is frac
          mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

          do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                  + PoroVolFourier%fractures(k) * divDensiteMolaireEnergieInterneSatFrac(j,m,k) / Delta_t
          end do

          bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
               - PoroVolFourier%fractures(k) * SmDensiteMolaireEnergieInterneSatFrac(m,k) / Delta_t
       end do

       do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
          JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
               + Poro_1VolFourier%fractures(k) * CpRoche * divTemperatureFrac(j,k) / Delta_t
       end do

       bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
            - Poro_1VolFourier%fractures(k) * CpRoche * SmTemperatureFrac(k) / Delta_t
#endif

    end do ! end of div prim n_k, frac

    ! 2.3 div prim: n_k(X_j^n), k is cell
    do k=1, NbCellLocal_Ncpus(commRank+1) ! cell

       rowk = k + NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1) ! row of cell k in vector (node own, frac own, cell)

       ! the diagonal element (k,k) is KacBigA(nz)
       do i=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
          if( (JacBigA%Num(i)-NbNodeLocal_Ncpus(commRank+1)-NbFracLocal_Ncpus(commRank+1))==k) then
             nz = i
             exit
          end if
       end do

       ! loop of composants in Ctilde, n_k is an unknown independent
       do i=1, NbCompCtilde_ctx(IncCell(k)%ic)
          icp = NumCompCtilde_ctx(i, IncCell(k)%ic)

          j = NbIncTotalPrim_ctx(IncCell(k)%ic) + i
          JacBigA%Val(j,icp,nz) = VolDarcy%cells(k) / Delta_t
       end do

       ! loop of composants, n_k is not an unknown independent
       do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k, k is cell
          mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

          do icp=1, NbComp
             if( MCP(icp,mph)==1) then ! Q_k \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                   JacBigA%Val(j,icp,nz) = JacBigA%Val(j,icp,nz) &
                        + divDensiteMolaireSatCompCell(j,icp,m,k) * PoroVolDarcy%cells(k) / Delta_t
                end do

                bigSm(icp,rowk) = bigSm(icp,rowk) &
                     - SmDensiteMolaireSatComp%cells(icp,m,k) * PoroVolDarcy%cells(k) / Delta_t

                ! if(commRank==0 .and. k==1 .and. s==1 .and. m==1) then
                !    print*, SmDensiteMolaireSatComp%cells(icp,m,k) * PoroVolDarcy%cells(k) / Delta_t
                ! end if

             end if
          end do

       end do ! end of loop Q_k

#ifdef _THERMIQUE_

       do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k, k is cell
          mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                  + PoroVolFourier%cells(k) * divDensiteMolaireEnergieInterneSatCell(j,m,k) / Delta_t
          end do

          bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
               - PoroVolFourier%cells(k) * SmDensiteMolaireEnergieInterneSatCell(m,k) / Delta_t

          ! if(commRank==0 .and. k==1 .and. s==1) then
          !    print*, PoroVolFourier%cells(k) * SmDensiteMolaireEnergieInterneSatCell(m,k) / Delta_t
          ! end if
       end do

       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
          JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
               + Poro_1VolFourier%cells(k) * CpRoche * divTemperatureCell(j,k) / Delta_t
       end do

       bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
            - Poro_1VolFourier%cells(k) * CpRoche * SmTemperatureCell(k) / Delta_t
#endif

    end do ! end of div prim n_k, cell

  end subroutine Jacobian_JacBigA_BigSm_accmolaire


  !> \brief Loop over the cells and node by cell to compute the Jacobian                         <br>
  !! loop of cell, index is k                                                                    <br>
  !! {                                                                                           <br>
  !!   1. loop of nodes s in cell k                                                             <br>
  !!      {                                                                                         <br>
  !!        1.1 div( B * DarcyFlux_{k,s}), k is cell, s is node;                                  <br>
  !!               where B is DensiteMolaire*Kr/Viso or DensiteMolaire/Viso*Enthalpie             <br>
  !!            this term has three contributions to Jacobian:                                      <br>
  !!               A_kr (k is row cell, r is col cell k)                                              <br>
  !!               A_kr (k is row cell, r is col node s)                                              <br>
  !!               A_kr (k is row cell, r is col nodes/fracs in cell k)                                <br>
  !!        1.2 A_sk, k is cell, s is node own;                                                            <br>
  !!      }                                                                                            <br>
  !!
  !!   2. loop of fracs in cell k                                                                        <br>
  !!      {                                                                                               <br>
  !!        A_ks, k is cell, s is frac;                                                                  <br>
  !!        A_sk, k is cell, s is frac own;                                                                  <br>
  !!      }                                                                                                  <br>
  !! }
  subroutine Jacobian_JacBigA_BigSm_cell

    ! div prims and Sm from term div(DensiteMolaire*Kr/Visco*Comp)*FluxDarcy
    double precision :: &
         divK1( NbIncTotalPrimMax, NbComp), & ! k for cell, represent k in paper
         divS1( NbIncTotalPrimMax, NbComp), & ! s for node/frac, represent s in paper
         Sm1( NbComp)

    ! div prims and Sm from term DensiteMolaire*Kr/Visco*Comp*div(FluxDarcy)
    ! three contributions from this item
    double precision :: &
         divK2( NbIncTotalPrimMax, NbComp), & ! A_kr, k is row cell, r is col cell k
         divS2( NbIncTotalPrimMax, NbComp), & ! A_kr, k is row cell, r is col col node s
         divR2( NbIncTotalPrimMax, NbComp, NbNodeCellMax+NbFracCellMax), & ! A_kr, k is row cell, r is nodes/fracs in cell k
         Sm2( NbComp)

    ! div prim of Darcy flux
    double precision :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), & !
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), &
         SmDarcyFlux( NbPhase)

#ifdef _THERMIQUE_
    ! div prim of Fourier flux
    double precision :: &
         divFourierFlux_k( NbIncTotalPrimMax), &
         divFourierFlux_s( NbIncTotalPrimMax), &
         divFourierFlux_r( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), & ! r represent s' in paper
         SmFourierFlux

    ! div prims and Sm from term div(DensiteMolaire*Kr*Enthalpie/Visco*FluxDarcy)
    double precision :: &
         divEgK( NbIncTotalPrimMax), & ! K for cell, represent k in paper
         divEgS( NbIncTotalPrimMax), & ! S for node/frac, represent s in paper
         divEgR( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), & ! R for node/frac, represent s' in paper
         SmEg
#endif

    integer :: k, s, nums, sf, r, numr, rf, i, j
    integer :: rows, cols, colr, nz
    integer :: m
    integer :: NbNodeCell, NbFracCell

    integer :: rowK, colk, &                    ! row/col of cell k in JacBigA
         rowSR( NbNodeCellMax+NbFracCellMax), & ! rows (in JacBigA) of nodes/frac in cell k
         colSR( NbNodeCellMax+NbFracCellMax)    ! cols (in JacBigA) of nodes/frac in cell k

    divK1(:,:) = 0.d0
    divS1(:,:) = 0.d0
    divK2(:,:) = 0.d0
    divS2(:,:) = 0.d0
    divR2(:,:,:) = 0.d0

    csrK(:) = 0
    csrSR(:) = 0

    ! main loop of cell
    do k=1, NbCellLocal_Ncpus(commRank+1)

       nbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
       nbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

       ! (rowk, colk) is row/col of cell k in JacBigA
       ! (rowSR, colSR) are rows/cols (in JacBigA) of nodes/frac connected to cell k

       call Jacobian_RowCol_KSR(k, nbNodeCell, nbFracCell, &
            rowk, colk, rowSR, colSR)

       do m=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
          csrK( JacBigA%Num(m) ) = m - JacBigA%Pt(rowk)
       end do

       ! two parts: 1. s is node; 2. s is frac

       ! 1. s is node
       do s=1, NbNodeCell
          nums = NodebyCellLocal%Num( NodebyCellLocal%Pt(k)+s) ! num node of s

          ! compute \sum{P_i \cap Q_{k/s} } div(DensiteMolaire*Kr/Viso) * DarcyFlux
          ! this term has two contributions: A_kr, r is k -> divK1
          !                                  A_kr, r is s -> divS1
          call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellnode(k,s,nums, &
               divK1, divS1, Sm1)

          ! compute div Darcy flux
          call Jacobian_divDarcyFlux_cellnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

          ! compute DensiteMolaire*Kr/Viso * div(Flux) using div(DarcyFlux)
          call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

          ! compute div (DensiteMolaire*Enthalpie/Viso * DarcyFlux) using div(DarcyFlux)
          call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

          ! compute div FluxFourier
          call Jacobian_divFourierFlux_cellnode(k,s,nums, &
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

          ! divK, divS, divR to JacBigA
          ! row of cell k then row of node s

          ! A_kk
          nz = JacBigA%Pt(rowk) + csrK(colk)
          do i=1, NbComp
             do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divK1(j,i) + divK2(j,i)
             end do
          end do

#ifdef _THERMIQUE_

          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                  + divEgK(j) + divFourierFlux_k(j)
          end do
#endif

          ! A_ks
          cols = colSR(s)
          nz = JacBigA%Pt(rowk) + csrK(cols)
          do i=1, NbComp
             do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divS1(j,i) + divS2(j,i)
             end do
          end do

#ifdef _THERMIQUE_

          ! ps. divFourierFlux_s = 0
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) + divEgS(j)
          end do
#endif

          ! A_kr, r is node
          do r=1, NbNodeCell ! r represent s' in paper, r is node
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             colr = colSR(r) ! col
             nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divR2(j,i,r)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + divEgR(j,r) + divFourierFlux_r(j,r)
             end do
#endif
          end do

          ! A_kr, r is frac
          do r=1, NbFracCell ! r represent s' in paper, r is frac
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             colr = colSR(rf)
             nz = JacBigA%Pt(rowk) + csrK(colr)
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divR2(j,i,rf)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + divEgR(j,rf) + divFourierFlux_r(j,rf)
             end do
#endif
          end do

          ! Sm
          bigSm(1:NbComp,rowk) = bigSm(1:NbComp,rowk) &
               - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_

          bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
               - SmEg - SmFourierFlux
#endif

          ! we consider the jacobian rows corresponding to node s when node s is own
          if( IdNodeLocal(nums)%Proc=="o") then

             rows = rowSR(s)
             do m=JacBigA%Pt(rows)+1, JacBigA%Pt(rows+1)
                csrSR( JacBigA%Num(m) ) = m - JacBigA%Pt(rows)
             end do

             ! A_sk, s is node, k is cell
             nz = JacBigA%Pt(rows) + csrSR(colk)

             if( IdNodeLocal(nums)%P /= "d" ) then
                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divK1(j,i) - divK2(j,i)
                   end do
                end do
             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgK(j) - divFourierFlux_k(j)
                end do
             end if
#endif

             ! A_ss, s is node
             cols = colSR(s)
             nz = JacBigA%Pt(rows) + csrSR(cols)

             if( IdNodeLocal(nums)%P /= "d" ) then
                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) &
                                             - divS1(j,i) - divS2(j,i)
                   end do
                end do
             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then

                ! ps. divFourierFlux_s=0
                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgS(j)
                end do
             end if
#endif

             ! A_sr, s is node, r is node
             do r=1, NbNodeCell ! r represent s' in paper
                numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                colr = colSR(r)
                nz = JacBigA%Pt(rows) + csrSR(colr)

                if( IdNodeLocal(nums)%P /= "d" ) then

                   do i=1, NbComp
                      do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                         JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divR2(j,i,r)
                      end do
                   end do
                end if

#ifdef _THERMIQUE_

                if( IdNodeLocal(nums)%T /= "d" ) then

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                           - divEgR(j,r) - divFourierFlux_r(j,r)
                   end do
                end if
#endif
             end do

             ! A_sr, s is node, r is frac
             do r=1, NbFracCell ! r represent s' in paper
                numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))
                rf = r + NbNodeCell

                colr = colSR(rf)
                nz = JacBigA%Pt(rows) + csrSR(colr)

                if( IdNodeLocal(nums)%P /= "d" ) then

                   do i=1, NbComp
                      do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                         JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divR2(j,i,rf)
                      end do
                   end do
                end if

#ifdef _THERMIQUE_

                if( IdNodeLocal(nums)%T /= "d" ) then

                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                           - divEgR(j,rf) - divFourierFlux_r(j,rf)
                   end do
                end if
#endif
             end do

             ! Sm
             if( IdNodeLocal(nums)%P /= "d" ) then
                bigSm(1:NbComp,rows) = bigSm(1:NbComp,rows) &
                     + Sm1(1:NbComp) + Sm2(1:NbComp)
             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then
                bigSm(NbComp+1,rows) = bigSm(NbComp+1,rows) &
                     + SmEg + SmFourierFlux
             end if
#endif

          end if

       end do ! end of row s

       ! if(k==1 .and. commRank==0) then
       !    do icp=1, NbComp
       !       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
       !          write(*,'(ES22.14)') JacBigA%Val(j,icp,nz)
       !       end do
       !       print*, ""
       !    end do
       ! end if


       ! 2. s is frac
       do s=1, NbFracCell
          nums = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+s)) ! nums is frac num

          sf = s + NbNodeCell

          ! compute div(DensiteMolaire*Kr/Viso) * DarcyFlux
          call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellfrac(k,s,nums, &
               divK1, divS1, Sm1)

          ! compute div Darcy flux
          call Jacobian_divDarcyFlux_cellfrac(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

          ! compute DensiteMolaire*Kr/Viso * div(Flux) using div(DarcyFlux)
          call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellfrac(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

          ! compute div (DensiteMolaire*Enthalpie/Viso * DarcyFlux) using div(DarcyFlux)
          call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellfrac(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

          ! compute div FluxFourier
          call Jacobian_divFourierFlux_cellfrac(k,s,nums, &
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

          ! divK, divS, divR to JacBigA
          ! row of cell k then row of frac s (nums)

          ! A_kk
          nz = JacBigA%Pt(rowk) + csrK(colk)
          do i=1, NbComp
             do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divK1(j,i) + divK2(j,i)
             end do
          end do

#ifdef _THERMIQUE_

          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                  + divEgK(j) + divFourierFlux_k(j)
          end do
#endif

          ! A_ks, s is frac
          cols = colSR(sf)
          nz = JacBigA%Pt(rowk) + csrK(cols)
          do i=1, NbComp
             do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divS1(j,i) + divS2(j,i)
             end do
          end do

#ifdef _THERMIQUE_

          ! ps. divFourierFlux_s = 0
          do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
             JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) + divEgS(j)
          end do
#endif

          ! A_kr, r is node
          do r=1, NbNodeCell ! r represent s' in paper, r is node
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             colr = colSR(r) ! col
             nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divR2(j,i,r)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + divEgR(j,r) + divFourierFlux_r(j,r)
             end do
#endif
          end do

          ! A_kr, r is frac
          do r=1, NbFracCell ! r represent s' in paper, r is frac
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             colr = colSR(rf)
             nz = JacBigA%Pt(rowk) + csrK(colr)
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divR2(j,i,rf)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + divEgR(j,rf) + divFourierFlux_r(j,rf)
             end do
#endif
          end do

          ! Sm
          bigSm(1:NbComp,rowk) = bigSm(1:NbComp,rowk) &
               - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_
          bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
               - SmEg - SmFourierFlux
#endif

          ! frac own
          ! its corresponding face num <= NbFaceOwn
          if( FracToFaceLocal(nums)<=NbFaceOwn_Ncpus(commRank+1)) then

             rows = rowSR(sf)

             do m=JacBigA%Pt(rows)+1, JacBigA%Pt(rows+1)
                csrSR( JacBigA%Num(m) ) = m - JacBigA%Pt(rows)
             end do

             ! A_sk, s is frac (own), k is cell
             nz = JacBigA%Pt(rows) + csrSR(colk)
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divK1(j,i) - divK2(j,i)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     - divEgK(j) - divFourierFlux_k(j)
             end do
#endif

             ! A_ss, s is frac (own)
             cols = colSR(sf)
             nz = JacBigA%Pt(rows) + csrSR(cols)

             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divS1(j,i) - divS2(j,i)
                end do
             end do

#ifdef _THERMIQUE_

             ! ps. divFourierFlux_s=0
             do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     - divEgS(j)
             end do
#endif

             ! A_sr, s is frac (own), r is node
             do r=1, NbNodeCell ! r represent s' in paper
                numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                colr = colSR(r)
                nz = JacBigA%Pt(rows) + csrSR(colr)

                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divR2(j,i,r)
                   end do
                end do

#ifdef _THERMIQUE_
                do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgR(j,r) - divFourierFlux_r(j,r)
                end do
#endif
             end do

             ! A_sr, s is frac (own), r is frac
             do r=1, NbFracCell ! r represent s' in paper
                numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))
                rf = r + NbNodeCell

                colr = colSR(rf)
                nz = JacBigA%Pt(rows) + csrSR(colr)

                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divR2(j,i,rf)
                   end do
                end do

#ifdef _THERMIQUE_

                do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgR(j,rf) - divFourierFlux_r(j,rf)
                end do
#endif
             end do

             ! Sm
             bigSm(1:NbComp,rows) = bigSm(1:NbComp,rows) &
                  + Sm1(1:NbComp) + Sm2(1:NbComp)

#ifdef _THERMIQUE_

             bigSm(NbComp+1,rows) = bigSm(NbComp+1,rows) &
                  + SmEg + SmFourierFlux
#endif

          end if ! end of frac won

       end do ! end of s frac

    end do ! end of cell k

  end subroutine Jacobian_JacBigA_BigSm_cell


  subroutine Jacobian_JacBigA_locate_frac_row(row)
  
    integer, intent(in) :: row

    integer :: rowk, colk, &      ! row/col of frac k in JacBigA
         rowSR( NbNodeFaceMax), & ! rows (in JacBigA) of nodes in frac k
         colSR( NbNodeFaceMax)    ! cols (in JacBigA) of nodes in frac k

    integer :: k, fk, i
    integer :: nbNodeFrac

    !open(UNIT=17,FILE='fractures.dat',status='old')
    !do k=1, NbFracLocal_Ncpus(commRank+1)
    !    fk = FracToFaceLocal(k) ! fk is face num
    !    nbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)
    !    write (17,*) "face", k, fk, NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+1:NodebyFaceLocal%Pt(fk+1))
    !end do
    !close(unit=17)

 do k=1, NbFracLocal_Ncpus(commRank+1)

       fk = FracToFaceLocal(k) ! fk is face num
       nbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

       !print *, "face", k, fk, NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+1:NodebyFaceLocal%Pt(fk+1))
       
       call Jacobian_RowCol_FR(k, nbNodeFrac, &
            rowk, colk, rowSR, colSR)

       if(row>=rowk .and. row<rowk+NbCompThermique) then
            print *, "row:", row, "frac:", k, "face:", fk
       end if

       do i=1, NbNodeFaceMax
           if(row>=rowSR(i) .and. row<rowSR(i)+NbCompThermique) then
                print *, "row:", row, "frac:", k, "face:", fk, "node", i, rowSR
           end if
       end do
       
   end do
       
  end subroutine Jacobian_JacBigA_locate_frac_row

  ! loop of frac, index is k
  ! {
  !   loop of nodes in frac k
  !   {
  !     A_ks, k is frac, s is node; (rowk/colk is row/col of cell k)
  !     A_sk, k is frac, s is node own;
  !   }
  ! }
  subroutine Jacobian_JacBigA_BigSm_frac


    ! div prims and Sm from term div(DensiteMolaire*Kr/Visco*Comp)*FluxDarcy
    double precision :: &
         divK1( NbIncTotalPrimMax, NbComp), & ! K for frac, represent k in paper
         divS1( NbIncTotalPrimMax, NbComp), & ! S for node, represent s in paper
         Sm1( NbComp)

    ! div prims and Sm from term DensiteMolaire*Kr/Visco*Comp*div(FluxDarcy)
    double precision :: &
         divK2( NbIncTotalPrimMax, NbComp), &
         divS2( NbIncTotalPrimMax, NbComp), &
         divR2( NbIncTotalPrimMax, NbComp, NbNodeFaceMax), & ! R for node/frac, represent s' in paper
         Sm2( NbComp)

    double precision :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeFaceMax), & ! r represen
         SmDarcyFlux( NbPhase)

#ifdef _THERMIQUE_

    ! div FluxFourier
    double precision :: &
         divFourierFlux_k( NbIncTotalPrimMax), &
         divFourierFlux_s( NbIncTotalPrimMax), &
         divFourierFlux_r( NbIncTotalPrimMax, NbNodeFaceMax), & ! r represent s' in paper
         SmFourierFlux

    ! div prims and Sm from term div(DensiteMolaire*Kr*Enthalpie/Visco*FluxDarcy)
    double precision :: &
         divEgK( NbIncTotalPrimMax), & ! K for frac, represent k in paper
         divEgS( NbIncTotalPrimMax), & ! S for node, represent s in paper
         divEgR( NbIncTotalPrimMax, NbNodeFaceMax), & ! R for node/frac, represent s' in paper
         SmEg
#endif

    integer :: rowK, colk, &      ! row/col of frac k in JacBigA
         rowSR( NbNodeFaceMax), & ! rows (in JacBigA) of nodes in frac k
         colSR( NbNodeFaceMax)    ! cols (in JacBigA) of nodes in frac k

    integer :: k, s, nums, r, numr, fk, i, j
    integer :: rows, cols, colr, nz
    integer :: m
    integer :: nbNodeFrac

    divK1(:,:) = 0.d0
    divS1(:,:) = 0.d0
    divK2(:,:) = 0.d0
    divS2(:,:) = 0.d0
    divR2(:,:,:) = 0.d0

    csrK(:) = 0
    csrSR(:) = 0

    do k=1, NbFracLocal_Ncpus(commRank+1)

       fk = FracToFaceLocal(k) ! fk is face num
       nbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

       call Jacobian_RowCol_FR(k, nbNodeFrac, &
            rowk, colk, rowSR, colSR)

       do m=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
          csrK( JacBigA%Num(m) ) = m - JacBigA%Pt(rowk) ! CHECKME: m: 1 -> NbCompThermique ?
       end do

       do s=1, nbNodeFrac
          nums = NodebyFaceLocal%Num( NodebyFaceLocal%Pt(fk)+s)

          ! compute div(DensiteMolaire*Kr/Viso) * DarcyFlux
          call Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_fracnode(k,s,nums,&
               divK1, divS1, Sm1)

          ! compute div Darcy flux
          call Jacobian_divDarcyFlux_fracnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux)

          ! compute DensiteMolaire*Kr/Viso * div(Flux) using div(DarcyFlux)
          call Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_fracnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divK2, divS2, divR2, Sm2)

#ifdef _THERMIQUE_

          ! compute div (DensiteMolaire*Enthalpie/Viso * DarcyFlux) using div(DarcyFlux)
          call Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_fracnode(k,s,nums, &
               divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
               divEgK, divEgS, divEgR, SmEg)

          ! compute div FluxFourier
          call Jacobian_divFourierFlux_fracnode(k,s,nums, &
               divFourierFlux_k, divFourierFlux_r, SmFourierFlux)
#endif

          ! line with frac own
          if(k<=NbFracOwn_Ncpus(commRank+1)) then

             ! A_kk
             nz = JacBigA%Pt(rowk) + csrK(colk)
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divK1(j,i) + divK2(j,i)
                end do
             end do

#ifdef _THERMIQUE_

             do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                     + divEgK(j) + divFourierFlux_k(j)
             end do
#endif

             ! A_ks
             cols = colSR(s)
             nz = JacBigA%Pt(rowk) + csrK(cols)
             do i=1, NbComp
                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                   JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divS1(j,i) + divS2(j,i)
                end do
             end do

#ifdef _THERMIQUE_

             ! ps. divFourierFlux_s = 0
             do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) + divEgS(j)
             end do
#endif

             ! A_ks'
             do r=1, NbNodeFrac ! r represent s' in paper, r is node

                numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

                colr = colSR(r) ! col
                nz = JacBigA%Pt(rowk) + csrK(colr) ! nnz
                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divR2(j,i,r)
                   end do
                end do

#ifdef _THERMIQUE_

                do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        + divEgR(j,r) + divFourierFlux_r(j,r)
                end do
#endif
             end do

             ! Sm
             bigSm(1:NbComp,rowk) = bigSm(1:NbComp,rowk) &
                  - Sm1(1:NbComp) - Sm2(1:NbComp)

#ifdef _THERMIQUE_

             bigSm(NbComp+1,rowk) = bigSm(NbComp+1,rowk) &
                  - SmEg - SmFourierFlux
#endif

          end if

          ! node s is own
          if(IdNodeLocal(nums)%Proc=="o") then

             rows = rowSR(s)
             do m=JacBigA%Pt(rows)+1, JacBigA%Pt(rows+1)
                csrSR( JacBigA%Num(m) ) = m - JacBigA%Pt(rows)
             end do

             ! A_sk
             nz = JacBigA%Pt(rows) + csrSR(colk)

             if( IdNodeLocal(nums)%P /= "d" ) then

                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divK1(j,i) - divK2(j,i)
                   end do
                end do
             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then

                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgK(j) - divFourierFlux_k(j)
                end do
             end if
#endif

             ! A_ss
             cols = colSR(s)
             nz = JacBigA%Pt(rows) + csrSR(cols)

             if( IdNodeLocal(nums)%P /= "d" ) then

                do i=1, NbComp
                   do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                      JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divS1(j,i) - divS2(j,i)
                   end do
                end do

             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then

                ! ps. divFourierFlux_s=0
                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                   JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                        - divEgS(j)
                end do
             end if
#endif

             ! A_ss', s' is node
             do r=1, NbNodeFrac ! r represent s' in paper
                numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

                colr = colSR(r)
                nz = JacBigA%Pt(rows) + csrSR(colr)

                if( IdNodeLocal(nums)%P /= "d" ) then
                   do i=1, NbComp
                      do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                         JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) - divR2(j,i,r)
                      end do
                   end do
                end if

#ifdef _THERMIQUE_

                if( IdNodeLocal(nums)%T /= "d" ) then
                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) &
                           - divEgR(j,r) - divFourierFlux_r(j,r)
                   end do
                end if
#endif
             end do

             ! Sm
             if( IdNodeLocal(nums)%P /= "d" ) then
                bigSm(1:NbComp,rows) = bigSm(1:NbComp,rows) &
                     + Sm1(1:NbComp) + Sm2(1:NbComp)
             end if

#ifdef _THERMIQUE_

             if( IdNodeLocal(nums)%T /= "d" ) then
                bigSm(NbComp+1,rows) = bigSm(NbComp+1,rows) &
                     + SmEg + SmFourierFlux
             end if
#endif

          end if

       end do ! s in frac k
    end do ! frac k

  end subroutine Jacobian_JacBigA_BigSm_frac

  ! loop of injection well
  ! (qw - qmol) * (Pwmax - Pw) = - qmol * (Pwmax-Pw) + qw * (Pwmax - Pw)
  !   where qw = sum_{s} sum_{i} q_{w,s,i}
  subroutine Jacobian_JacBigA_BigSm_wellinj

    integer :: k, rowk, colk, s, nums, rows, cols, m, icp, nz
    double precision :: Tws, Ps_Pws, Ts, WIDws, WIFws
    double precision :: dP_w(NbComp), dP_s(NbComp), dP_ER_w, dP_ER_s
    logical :: something_is_injected

    nz = -1
    
    do k=1, NbWellInjLocal_Ncpus(commRank+1)

              something_is_injected = .false.

              ! A_kk, k is well
       rowk = k + NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1) &
            + NbCellLocal_Ncpus(commRank+1)
       colk = k + NbNodeLocal_Ncpus(commRank+1) + NbFracLocal_Ncpus(commRank+1) &
            + NbCellLocal_Ncpus(commRank+1)

       ! assembly equaitons of own wells
       if(k<=NbWellInjOwn_Ncpus(commRank+1)) then
          do m=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
             csrK( JacBigA%Num(m) ) = m - JacBigA%Pt(rowk)
          end do
       end if

       ! Pwmax_Pw = DataWellInjLocal(k)%PressionMax - IncPressionWellInj(k) ! P_w^{max}-P_w

       ! Two cases:
       !   if indwell=="p": well equation is Pwmax-Pw = 0
       !   if indwell=="f": well equation is qw-qmol = 0, where qw = sum_{s} sum_{i} q_{w,s,i}

       ! Two steps:
       !   Step 1. well equation if indwell=="p"
       !   Step 2. loop of nodes of well k
       !           - well equation if indwell=="f"
       !           - reservoir equation

       ! Step 1. well equation: Pwmax - Pw = 0
       if( DataWellInjLocal(k)%IndWell == 'p') then

          if(k<=NbWellInjOwn_Ncpus(commRank+1)) then ! own injection well
             nz = JacBigA%Pt(rowk) + csrK(colk)
             JacBigA%Val(1,1,nz) = -1.d0
          end if
          something_is_injected = .true.
       end if

       ! Step 2.

       ! nodes in well k
       do s = NodebyWellInjLocal%Pt(k)+1, NodebyWellInjLocal%Pt(k+1)
          nums = NodebyWellInjLocal%Num(s) ! num node

          Ps_Pws = IncNode(nums)%Pression - PerfoWellInj(s)%Pression ! P_s - P_{w,s}
          Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
          Ts = IncNode(nums)%Temperature    ! T_s

          WIDws = NodeDatabyWellInjLocal%Val(s)%WID
          WIFws = NodeDatabyWellInjLocal%Val(s)%WIF

          dP_w(:) = 0.d0
          dP_s(:) = 0.d0
#ifdef _THERMIQUE_
          dP_ER_w = 0.d0
          dP_ER_s = 0.d0
#endif

          if(IdNodeLocal(nums)%Proc=="o") then ! ps. this node can be in the boundary
             rows = nums
             do m=JacBigA%Pt(rows)+1, JacBigA%Pt(rows+1)
                csrSR( JacBigA%Num(m) ) = m - JacBigA%Pt(rows)
             end do
          end if

          cols = nums
          if(Ps_Pws < 0.d0) then ! if >0, this term is 0
  
              something_is_injected = .true.

             do icp=1, NbComp
                dP_w(icp) = divDensiteMolaireKrViscoCompWellInj(icp,s) * WIDws * Ps_Pws &
                     - DensiteMolaireKrViscoCompWellInj(icp,s) * WIDws
             end do

             do icp=1, NbComp
                dP_s(icp) = DensiteMolaireKrViscoCompWellInj(icp,s) * WIDws
             end do

#ifdef _THERMIQUE_
             dP_ER_w = divDensiteMolaireKrViscoEnthalpieWellInj(s) * WIDws * Ps_Pws &
                  - DensiteMolaireKrViscoEnthalpieWellInj(s) * WIDws
             dP_ER_s = DensiteMolaireKrViscoEnthalpieWellInj(s) * WIDws
#endif

             if(k<=NbWellInjOwn_Ncpus(commRank+1)) then ! own injection well

                if( DataWellInjLocal(k)%IndWell == 'f') then

                   ! A_kk, k is own injection well
                   nz = JacBigA%Pt(rowk) + csrK(colk)
                   do icp=1, NbComp
                      JacBigA%Val(1,1,nz) = JacBigA%Val(1,1,nz) + dP_w(icp)
                   end do

                   ! A_ks, k is own injection well, s is node
                   nz = JacBigA%Pt(rowk) + csrK(cols)
                   do icp=1, NbComp
                      JacBigA%Val(1,1,nz) = JacBigA%Val(1,1,nz) + dP_s(icp)
                   end do
                end if
             end if

             if(IdNodeLocal(nums)%Proc=="o") then ! ps. this node can be in the boundary

                ! Ask, s is node, k is injection well
                nz = JacBigA%Pt(rows) + csrSR(colk)
                JacBigA%Val(1,1:NbComp,nz) = JacBigA%Val(1,1:NbComp,nz) + dP_w(:) ! term q_{w,s,i}, derivative of P_w

#ifdef _THERMIQUE_
                JacBigA%Val(1,NbComp+1,nz) = JacBigA%Val(1,NbComp+1,nz) + dP_ER_w
#endif

                ! Ass, s is node
                nz = JacBigA%Pt(rows) + csrSR(cols)
                JacBigA%Val(1,1:NbComp,nz) = JacBigA%Val(1,1:NbComp,nz) + dP_s(:) ! term q_{w,s,i}, derivative of P_w

#ifdef _THERMIQUE_
                JacBigA%Val(1,NbComp+1,nz) = JacBigA%Val(1,NbComp+1,nz) + dP_ER_s
#endif
             end if

          end if

!           ! (3) WIF_{w,s} (T_s - T_{w,s})
! #ifdef _THERMIQUE_
!           if(IdNodeLocal(nums)%Proc=="o") then ! ps. this node can be in the boundary

!              ! Ass, s is node
!              nz = JacBigA%Pt(rows) + csrSR(cols)
!              JacBigA%Val(:,NbComp+1,nz) = JacBigA%Val(:,NbComp+1,nz) + divTemperatureNode(:,nums)

!              bigSm(NbComp+1,nums) = bigSm(NbComp+1,nums) - WIFws * SmTemperatureNode(nums)
!           end if
! #endif
       end do

       if(.not.something_is_injected) write(*,*) 'WARNING: nothing is injected in well', k, 'on proc', commRank+1
          
    end do

  end subroutine Jacobian_JacBigA_BigSm_wellinj


  ! loop of production well
  ! (qmol - qw) * (Pw - Pwmin) = qmol * (Pw-Pwmin) - qw * (Pw - Pwmin)
  ! where qw = sum_{s} sum_{i} q_{w,s,i}
  subroutine Jacobian_JacBigA_BigSm_wellprod

    integer :: k, rowk, colk, s, nums, rows, cols, m, mph, n, icp, nz
    double precision :: Pws, Ps, WIDws, WIFws, Ps_Pws
    double precision :: &
         dP_w(NbComp), dP_s(NbCompThermique,NbComp), &
         dP_ER_w, der_ER_s(NbCompThermique)
    logical :: something_is_produced

    nz = -1
    
    do k=1, NbWellProdLocal_Ncpus(commRank+1)
    something_is_produced = .false.
       ! A_kk, k is well
       rowk = k + NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1) &
            + NbCellLocal_Ncpus(commRank+1) + NbWellInjOwn_Ncpus(commRank+1)
       colk = k + NbNodeLocal_Ncpus(commRank+1) + NbFracLocal_Ncpus(commRank+1) &
            + NbCellLocal_Ncpus(commRank+1) + NbWellInjLocal_Ncpus(commRank+1)

       ! assembly equaitons of own wells
       if(k<=NbWellProdOwn_Ncpus(commRank+1)) then
          do n=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
             csrK( JacBigA%Num(n) ) = n - JacBigA%Pt(rowk)
          end do
       end if

       ! Pw_Pwmin = IncPressionWellProd(k) - DataWellProdLocal(k)%PressionMin ! P_w - P_w^{min}

       ! Two cases:
       !   if indwell=="p": well equation is Pw-Pwmin = 0
       !   if indwell=="f": well equation is qmol-qw = 0, where qw = sum_{s} sum_{i} q_{w,s,i}

       ! Step 1. well equation if indwell=="p"
       ! Step 2. loop of nodes of well k
       !         - well equation if indwell=="f"
       !         - reservoir equation

       ! Step 1. well equation: Pw - Pwmin = 0
       if( DataWellProdLocal(k)%IndWell == 'p') then
          something_is_produced = .true.
          if(k<=NbWellProdOwn_Ncpus(commRank+1)) then
             nz = JacBigA%Pt(rowk) + csrK(colk)
             JacBigA%Val(1,1,nz) = 1.d0
          end if
       end if

       ! Step 2.

       ! nodes of well k
       do s=NodebyWellProdLocal%Pt(k)+1, NodebyWellProdLocal%Pt(k+1)
          nums = NodebyWellProdLocal%Num(s) ! nums is node num

          Ps = IncNode(nums)%Pression       ! P_s
          Pws = PerfoWellProd(s)%Pression   ! P_{w,s}
          Ps_Pws = Ps - Pws

          WIDws = NodeDatabyWellProdLocal%Val(s)%WID ! WID_{w,s}
#ifdef _THERMIQUE_
          WIFws = NodeDatabyWellProdLocal%Val(s)%WIF ! WIF_{w,s}
#endif

          if(Ps_Pws > 0.d0) then ! if Ps_Pws < 0 then this term is zero
            something_is_produced = .true.
             ! derivative of
             !   sum_{Q_s \cap P_i} q_{w,s,i}
             !   sum_{Q_s \cap P_i} q_{w,s,e}
             dP_w(:) = 0.d0
             dP_s(:,:) = 0.d0
#ifdef _THERMIQUE_
             dP_ER_w = 0.d0
             der_ER_s(:) = 0.d0
#endif

             do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
                mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

                do icp=1, NbComp
                   if(MCP(icp,mph)==1) then ! \cap P_i

                      dP_w(icp) = dP_w(icp) - DensiteMolaireKrViscoCompNode(icp,m,nums) * WIDws

                      dP_s(:,icp) = dP_s(:,icp) + divDensiteMolaireKrViscoCompNode(:,icp,m,nums) * WIDws * Ps_Pws &
                           + DensiteMolaireKrViscoCompNode(icp,m,nums) * WIDws

                   end if
                end do

#ifdef _THERMIQUE_
                dP_ER_w = dP_ER_w - DensiteMolaireKrViscoEnthalpieNode(m,nums) * WIDws
                der_ER_s(:) = der_ER_s(:) + divDensiteMolaireKrViscoEnthalpieNode(:,m,nums) * WIDws * Ps_Pws &
                     + DensiteMolaireKrViscoEnthalpieNode(m,nums) * WIDws
#endif
             end do

             if(k<=NbWellProdOwn_Ncpus(commRank+1)) then ! own production well

                if( DataWellProdLocal(k)%IndWell == 'f') then

                   ! A_kk, k is own production well
                   nz = JacBigA%Pt(rowk) + csrK(colk)
                   do icp=1, NbComp
                      JacBigA%Val(1,1,nz) = JacBigA%Val(1,1,nz) - dP_w(icp) ! term -\sum_{i} q_{w,s,i}, derivative of P_w
                   end do

                   ! A_ks, k is own production well, s is node
                   nz = JacBigA%Pt(rowk) + csrK(nums)
                   do icp=1, NbComp
                      JacBigA%Val(:,1,nz) = JacBigA%Val(:,1,nz) - dP_s(:,icp) ! term -\sum_{i} q_{w,s,i}, derivative of node s
                   end do
                end if
             end if

             if(IdNodeLocal(nums)%Proc=="o") then ! node own, this node can not be in the boundary (why ?)

                rows = nums
                cols = nums
                do n=JacBigA%Pt(rows)+1, JacBigA%Pt(rows+1)
                   csrSR( JacBigA%Num(n) ) = n - JacBigA%Pt(rows)
                end do

                ! Ask, s is node, k is production well
                nz = JacBigA%Pt(rows) + csrSR(colk)
                do icp=1, NbComp
                   JacBigA%Val(1,icp,nz) = JacBigA%Val(1,icp,nz) + dP_w(icp) ! term q_{w,s,i}, derivative of P_w
                end do
#ifdef _THERMIQUE_
                JacBigA%Val(1,NbComp+1,nz) = JacBigA%Val(1,NbComp+1,nz) + dP_ER_w
#endif

                ! Ass, s is node
                nz = JacBigA%Pt(rows) + csrSR(cols)
                do icp=1, NbComp
                   JacBigA%Val(:,icp,nz) = JacBigA%Val(:,icp,nz) + dP_s(:,icp) ! term q_{w,s,i}, derivative of node s
                end do
#ifdef _THERMIQUE_
                JacBigA%Val(:,NbComp+1,nz) = JacBigA%Val(:,NbComp+1,nz) + der_ER_s(:)
#endif
             end if

          end if
       end do
       if(.not.something_is_produced) write(*,*) 'WARNING: nothing is produced from well', k, 'on proc', commRank+1
    end do

  end subroutine Jacobian_JacBigA_BigSm_wellprod

#ifdef _WIP_FREEFLOW_STRUCTURES_
  ! loop of node, index is nums
  ! 1.1 div 
  subroutine Jacobian_JacBigA_BigSm_FF_node
    ! div prims and Sm from FreeFlow term 
    double precision :: &
         divS3( NbIncTotalPrimMax, NbComp), & ! s for node, represent s in paper
         Sm3( NbComp), &
         divTFF( NbIncTotalPrimMax), &
         SmTFF
    integer :: nums, i, j, nz
   
    do nums=1, NbNodeOwn_Ncpus(commRank+1)

      if(IdFFNodeLocal(nums)) then ! loop over freeflow dof only

         ! compute the contribution of the freeflow
         call Jacobian_divMolarFreeFlow_node(nums,divS3,Sm3)
#ifdef _THERMIQUE_
         ! compute the thermal contribution of the freeflow
         call Jacobian_divThermalFreeFlow_node(nums,divTFF,SmTFF)
#endif

         ! the diagonal element (s,s) is JacBigA(nz)
         do i=JacBigA%Pt(nums)+1, JacBigA%Pt(nums+1)
            if( JacBigA%Num(i)==nums) then
               nz = i
               exit
            end if
         end do

         if( IdNodeLocal(nums)%P /= "d" ) then
            ! JacBigA%Val(:,:,nz) 
            do i=1, NbComp
               do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  JacBigA%Val(j,i,nz) = JacBigA%Val(j,i,nz) + divS3(j,i)
               end do
            end do
            ! Sm
            bigSm(1:NbComp,nums) = bigSm(1:NbComp,nums) - Sm3(1:NbComp)

         end if ! Dirichlet node

#ifdef _THERMIQUE_
         if( IdNodeLocal(nums)%T /= "d" ) then
            ! ps. divFourierFlux_s=0
            do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
               JacBigA%Val(j,NbComp+1,nz) = JacBigA%Val(j,NbComp+1,nz) + divTFF(j)
            end do
            ! Sm
            bigSm(NbComp+1,nums) = bigSm(NbComp+1,nums) - SmTFF
         end if ! Dirichlet node
#endif

      endif ! FreeFlow node

    enddo ! node nums


  end subroutine Jacobian_JacBigA_BigSm_FF_node


  ! Derivatives of the FreeFlow terms in the molar balance equations
  subroutine Jacobian_divMolarFreeFlow_node(nums, divS, Sm0)

    integer, intent(in) :: nums

    double precision, intent(out) :: &
         divS( NbIncTotalPrimMax, NbComp), &
         Sm0 ( NbComp)

    ! tmp
    integer :: m, mph, icp, j

    divS(:,:) = 0.d0
    Sm0(:) = 0.d0

    ! -> divS, node
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
      mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

      if(IncNode(nums)%FreeFlow_flowrate(mph)>=0.d0) then

         ! To understand better, change the order of the loop do m=.. and the loop do icp=..
         do icp=1, NbComp
            if(MCP(icp,mph)==1) then ! \cap P_i

               do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  divS(j,icp) = divS(j,icp) + SurfFreeFlowLocal(nums) * ( &
                        divFreeFlowMolarFlowrateCompNode(j,icp,m,nums) +  &
                        divFreeFlowHmCompNode(j,icp,m,nums) &
                        )
               end do

               ! Sm0
               Sm0(icp) = Sm0(icp) + SurfFreeFlowLocal(nums) * ( &
                     SmFreeFlowMolarFlowrateCompNode(icp,m,nums) + &
                     SmFreeFlowHmCompNode(icp,m,nums) &
                     )
            end if
         end do ! end of icp

      else ! IncNode(nums)%FreeFlow_flowrate(mph)<0.d0
      ! liq phase never enters in this loop because always FreeFlow_flowrate(liq)>=0.d0

         ! To understand better, change the order of the loop do m=.. and the loop do icp=..
         do icp=1, NbComp
            if(MCP(icp,mph)==1) then ! \cap P_i

               do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                  divS(j,icp) = divS(j,icp) + SurfFreeFlowLocal(nums) * ( & ! only gas phase because FreeFlow_flowrate(liq)>=0.d0
                        divFreeFlowMolarFlowrateNode(j,m,nums) * atm_comp(icp,mph) + & ! atm_comp is a constant, no derivative wrt atm_comp
                        divFreeFlowHmCompNode(j,icp,m,nums) &
                        )
               end do

               ! Sm0
               Sm0(icp) = Sm0(icp) + SurfFreeFlowLocal(nums) * ( &
                     SmFreeFlowMolarFlowrateNode(m,nums) * atm_comp(icp,mph) + &
                     SmFreeFlowHmCompNode(icp,m,nums) &
                     )
               end if
         end do ! end of icp
      endif

    enddo

  end subroutine Jacobian_divMolarFreeFlow_node

  ! Derivatives of the FreeFlow terms in the energy balance equation
  subroutine Jacobian_divThermalFreeFlow_node(nums, divS, Sm0)

    integer, intent(in) :: nums

    double precision, intent(out) :: &
         divS( NbIncTotalPrimMax), &
         Sm0

    !  tmp
    integer :: m, mph, j

    divS = 0.d0
    Sm0 = 0.d0

    ! -> divS, node
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
      mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

      if(IncNode(nums)%FreeFlow_flowrate(mph)>=0.d0) then

         do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
            divS(j) = divS(j) + SurfFreeFlowLocal(nums) * &
                      divFreeFlowMolarFlowrateEnthalpieNode(j,m,nums)
         enddo
         Sm0 = Sm0 + SurfFreeFlowLocal(nums) * SmFreeFlowMolarFlowrateEnthalpieNode(m,nums)
         
      else ! IncNode(nums)%FreeFlow_flowrate(mph)<0.d0
      ! liq phase never enters in this loop because always FreeFlow_flowrate(liq)>=0.d0
         do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
            divS(j) = divS(j) + SurfFreeFlowLocal(nums) * &
                      divFreeFlowMolarFlowrateNode(j,m,nums) * AtmEnthalpieNode(m,nums) ! AtmEnthalpieNode is a constant, no derivative
         enddo
         Sm0 = Sm0 + SurfFreeFlowLocal(nums) * &
               SmFreeFlowMolarFlowrateNode(m,nums) * AtmEnthalpieNode(m,nums)

      endif ! sign of flux
    enddo

    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
      divS(j) = divS(j) + SurfFreeFlowLocal(nums) * &
                divFreeFlowHTTemperatureNetRadiationNode(j,nums)
    enddo
    Sm0 = Sm0 + SurfFreeFlowLocal(nums) * SmFreeFlowHTTemperatureNetRadiationNode(nums)

  end subroutine Jacobian_divThermalFreeFlow_node
#endif

  ! term: \sum{P_i \cap Q_{k or s} } &
  !          (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI)
  !       = divK * div(X_k) + divS * div(X_s) + Sm
  !       where k is cell, s is node
  ! compute:
  !       divK, divS, Sm
  !
  ! if FluxDarcyKI>=0 then
  !    \sum{P_i \cap Q_{k} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
  ! else
  !    \sum{P_i \cap Q_{s} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
  !
  ! it is equivalent to:
  !    \sum{P_i \cap Q_{k} \cap FluxDarcyKI>=0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
  !    \sum{P_i \cap Q_{s} \cap FluxDarcyKI<0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
  subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellnode(k,s,nums, &
       divK, divS, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         Sm0 ( NbComp)

    ! tmp
    integer :: m, mph, icp, j

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    Sm0(:) = 0.d0

    ! -> divK, upwind is k, divS=0
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       ! if(commRank==0 .and. k==1 .and. s==1) then
       !    ! print*, divDensiteMolaireKrViscoCompCell(:,:,2,1)
       !    ! print*, FluxDarcyKI(mph,s,k)
       ! end  if


       if(FluxDarcyKI(mph,s,k)>=0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=...
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                   divK(j,icp) = divK(j,icp) + &
                        divDensiteMolaireKrViscoCompCell(j,icp,m,k) * FluxDarcyKI(mph,s,k)

                end do

                Sm0(icp) = Sm0(icp) + &
                     SmDensiteMolaireKrViscoCompCell(icp,m,k) * FluxDarcyKI(mph,s,k)
             end if
          end do ! end of icp

       end if
    end do ! end of Q_k


    ! -> divS, upwind is node, divK=0
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyKI(mph,s,k)<0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=..
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                   divS(j,icp) = divS(j,icp) + &
                        divDensiteMolaireKrViscoCompNode(j,icp,m,nums) * FluxDarcyKI(mph,s,k)
                end do

                ! Sm0
                ! if nums is dirichlet, sm is supposed to be null
                if( IdNodeLocal(nums)%P /= "d" ) then
                   Sm0(icp) = Sm0(icp) + &
                        SmDensiteMolaireKrViscoCompNode(icp,m,nums) * FluxDarcyKI(mph,s,k)
                end if

             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellnode



  ! term: \sum{P_i \cap Q_{k or s} } &
  !          (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI)
  !       = divK * div(X_k) + divS * div(X_s) + Sm
  !       where k is cell, s is frac
  ! compute:
  !       divK, divS, Sm
  !
  ! if FluxDarcyKI>=0 then
  !    \sum{P_i \cap Q_{k} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
  ! else
  !    \sum{P_i \cap Q_{s} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
  !
  ! it is equivalent to:
  !    \sum{P_i \cap Q_{k} \cap FluxDarcyKI>=0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divK
  !    \sum{P_i \cap Q_{s} \cap FluxDarcyKI<0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyKI) -> divS
  subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellfrac(k,s,nums, &
       divK, divS, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         Sm0( NbComp)

    ! tmp
    integer :: m, mph, icp, j, sf

    ! sf = s + NbNodeCell
    sf = s + NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    Sm0(:) = 0.d0

    ! -> divK, upwind is k, divS=0
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       if(FluxDarcyKI(mph,sf,k)>=0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=..
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
                   divK(j,icp) = divK(j,icp) + &
                        divDensiteMolaireKrViscoCompCell(j,icp,m,k) * FluxDarcyKI(mph,sf,k)

                end do

                Sm0(icp) = Sm0(icp) + &
                     SmDensiteMolaireKrViscoCompCell(icp,m,k) * FluxDarcyKI(mph,sf,k)
             end if
          end do ! end of icp

       end if
    end do ! end of Q_k


    ! -> divS, upwind is frac, divK=0
    do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

       if(FluxDarcyKI(mph,sf,k)<0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=..
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
                   divS(j,icp) = divS(j,icp) + &
                        divDensiteMolaireKrViscoCompFrac(j,icp,m,nums) * FluxDarcyKI(mph,sf,k)
                end do

                Sm0(icp) = Sm0(icp) + &
                     SmDensiteMolaireKrViscoCompFrac(icp,m,nums) * FluxDarcyKI(mph,sf,k)

             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_cellfrac


  ! term: \sum{P_i \cap Q_{k or s} } &
  !          (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyFI)
  !       = divK * div(X_k) + divS * div(X_s) + Sm
  !       where k is frac, s is node
  ! compute:
  !       divK, divS, Sm
  !
  ! if FluxDarcyFI>=0 then
  !    \sum{P_i \cap Q_{k} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyFI) -> divK
  ! else
  !    \sum{P_i \cap Q_{s} } (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyFI) -> divS
  !
  ! it is equivalent to:
  !    \sum{P_i \cap Q_{k} \cap FluxDarcyFI>=0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyFI) -> divK
  !    \sum{P_i \cap Q_{s} \cap FluxDarcyFI<0} (div (DensiteMolaire * PermRel / Viscosite * Comp ) * FluxDarcyFI) -> divS
  subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_fracnode(k,s,nums, &
       divK, divS, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         Sm0( NbComp)

    ! tmp
    integer :: m, mph, icp, j

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    Sm0(:) = 0.d0

    ! -> divK, upwind is k, divS=0
    do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is frac
       mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

       if(FluxDarcyFI(mph,s,k)>=0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=..
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
                   divK(j,icp) = divK(j,icp) + &
                        divDensiteMolaireKrViscoCompFrac(j,icp,m,k) * FluxDarcyFI(mph,s,k)

                end do

                Sm0(icp) = Sm0(icp) + &
                     SmDensiteMolaireKrViscoCompFrac(icp,m,k) * FluxDarcyFI(mph,s,k)
             end if
          end do ! end of icp

       end if
    end do ! end of Q_k


    ! -> divS, upwind is node, divK=0
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s, s is node
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyFI(mph,s,k)<0.d0) then

          ! To understand better, change the order of the loop do m=.. and the loop do icp=..
          do icp=1, NbComp
             if(MCP(icp,mph)==1) then ! \cap P_i

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
                   divS(j,icp) = divS(j,icp) + &
                        divDensiteMolaireKrViscoCompNode(j,icp,m,nums) * FluxDarcyFI(mph,s,k)
                end do

                if(IdNodeLocal(nums)%P /= "d") then
                   Sm0(icp) = Sm0(icp) + &
                        SmDensiteMolaireKrViscoCompNode(icp,m,nums) * FluxDarcyFI(mph,s,k)
                end if

             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_divDensiteMolaireKrViscoComp_DarcyFlux_fracnode



  ! term: \sum{P_i \cap (Q_k \cap Q_s)} &
  !           DensiteMolaire * PermRel / Viscosite * Comp * div(FluxDarcyKI)
  !       = ( divK * div(X_k)
  !         + divS * div(X_s)
  !         + \sum_r div(X_r) * div(X_r)
  !         + Sm
  ! compute
  !        divK, divS, divR, Sm
  subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellnode(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divK, divS, divR, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), & ! r represen
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         divR( NbIncTotalPrimMax, NbComp, NbNodeCellMax+NbFracCellMax), &
         Sm0 ( NbComp)

    ! tmp
    integer :: m, mph, j, icp, r, numr, rf
    integer :: NbNodeCell, NbFracCell

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    divR(:,:,:) = 0.d0
    Sm0(:) = 0.d0

    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    ! sum_{P_i \cap Q_{k} \cap FluxDarcyKI>=0}
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       if(FluxDarcyKI(mph,s,k)>=0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompCell(icp,m,k)
                end do

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompCell(icp,m,k)
                end do

                do r=1, NbNodeCell ! divR for r is node in dof(k)
                   numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) =  divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompCell(icp,m,k)
                   end do
                end do

                do r=1, NbFracCell ! divR for r is frac in dof(k)
                   numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))
                   rf = r + NbNodeCell

                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      divR(j,icp,rf) =  divR(j,icp,rf) &
                           + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoCompCell(icp,m,k)
                   end do
                end do

                ! Sm
                Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompCell(icp,m,k)

             end if
          end do ! end of icp

       end if
    end do ! end of Q_k

    ! sum_{P_i \cap Q_{s} \cap FluxDarcyKI<0}
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyKI(mph,s,k)<0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end do

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end do

                ! divR for r is node in dof(K)
                do r=1, NbNodeCell
                   numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) = divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                   end do
                end do

                ! divR for r is frac in dof(K)
                do r=1, NbFracCell
                   numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))
                   rf = r + NbNodeCell

                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      divR(j,icp,rf) = divR(j,icp,rf) &
                           + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                   end do
                end do

                ! Sm
                ! if nums is dirichlet, sm is supposed to be null
                if(IdNodeLocal(nums)%P /= "d") then
                   Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end if
             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellnode



  ! term: \sum{P_i \cap (Q_k \cap Q_s)} &
  !           DensiteMolaire * PermRel / Viscosite * Comp * div(FluxDarcyKI)
  !       = ( divK * div(X_k)
  !         + divS * div(X_s)
  !         + \sum_r div(X_r) * div(X_r)
  !         + Sm )
  ! compute:
  !        divK, divS, divR, Sm
  subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellfrac(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divK, divS, divR, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), & ! r represen
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         divR( NbIncTotalPrimMax, NbComp, NbNodeCellMax+NbFracCellMax), &
         Sm0 ( NbComp)

    ! tmp
    integer :: m, mph, j, icp, r, numr, rf, sf
    integer :: NbNodeCell, NbFracCell

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    divR(:,:,:) = 0.d0
    Sm0(:) = 0.d0

    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    sf = s + NbNodeCell

    ! sum_{P_i \cap Q_{k} \cap FluxDarcyKI>=0}
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       if(FluxDarcyKI(mph,sf,k)>=0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompCell(icp,m,k)
                end do

                do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompCell(icp,m,k)
                end do

                do r=1, NbNodeCell ! divR for r is node in dof(k)
                   numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) =  divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompCell(icp,m,k)
                   end do
                end do

                do r=1, NbFracCell ! divR for r is frac in dof(k)
                   numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))

                   rf = r + NbNodeCell

                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      divR(j,icp,rf) =  divR(j,icp,rf) &
                           + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoCompCell(icp,m,k)
                   end do
                end do

                ! Sm
                Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompCell(icp,m,k)

             end if
          end do ! end of icp

       end if
    end do ! end of Q_k

    ! sum_{P_i \cap Q_{s} \cap FluxDarcyKI<0}
    do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

       if(FluxDarcyKI(mph,sf,k)<0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompFrac(icp,m,nums)
                end do

                do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompFrac(icp,m,nums)
                end do

                ! divR for r is node in dof(K)
                do r=1, NbNodeCell
                   numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) = divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompFrac(icp,m,nums)
                   end do
                end do

                ! divR for r is frac in dof(K)
                do r=1, NbFracCell
                   numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r))

                   rf = r + NbNodeCell

                   do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                      divR(j,icp,rf) = divR(j,icp,rf) &
                           + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoCompFrac(icp,m,nums)
                   end do
                end do

                ! Sm
                Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompFrac(icp,m,nums)

             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_cellfrac


  subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_fracnode(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divK, divS, divR, Sm0)

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeFaceMax), & ! r represen
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divK( NbIncTotalPrimMax, NbComp), &
         divS( NbIncTotalPrimMax, NbComp), &
         divR( NbIncTotalPrimMax, NbComp, NbNodeFaceMax), &
         Sm0 ( NbComp)

    ! tmp
    integer :: m, mph, j, icp, r, numr, fk
    integer :: NbNodeFrac

    divK(:,:) = 0.d0
    divS(:,:) = 0.d0
    divR(:,:,:) = 0.d0
    Sm0(:) = 0.d0

    fk = FracToFaceLocal(k) ! face num
    NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

    ! sum_{P_i \cap Q_{k} \cap FluxDarcyKI>=0}
    do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

       if(FluxDarcyFI(mph,s,k)>=0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompFrac(icp,m,k)
                end do

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompFrac(icp,m,k)
                end do

                do r=1, NbNodeFrac ! divR for r is node in dof(k)
                   numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) =  divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompFrac(icp,m,k)
                   end do
                end do

                ! Sm
                Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompFrac(icp,m,k)

             end if
          end do ! end of icp

       end if
    end do ! end of Q_k

    ! sum_{P_i \cap Q_{s} \cap FluxDarcyKI<0}
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyFI(mph,s,k)<0.d0) then

          do icp=1, NbComp ! P_i
             if(MCP(icp,mph)==1) then

                do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divK
                   divK(j,icp) = divK(j,icp) &
                        + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end do

                do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
                   divS(j,icp) = divS(j,icp) &
                        + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end do

                ! divR for r is node in dof(K)
                do r=1, NbNodeFrac
                   numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

                   do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                      divR(j,icp,r) = divR(j,icp,r) &
                           + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                   end do
                end do

                ! Sm
                if(IdNodeLocal(nums)%P /= "d") then
                   Sm0(icp) = Sm0(icp) + SmDarcyFlux(mph) * DensiteMolaireKrViscoCompNode(icp,m,nums)
                end if

             end if
          end do ! end of icp

       end if

    end do ! end of Q_s

  end subroutine Jacobian_DensiteMolaireKrViscoComp_divDarcyFlux_fracnode



  subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellnode(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divEgK, divEgS, divEgR, SmEg)

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), & ! r represen
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divEgK( NbIncTotalPrimMax), &
         divEgS( NbIncTotalPrimMax), &
         divEgR( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), &
         SmEg

    ! tmp
    integer :: m, mph, j, r, numr, rf
    integer :: NbNodeCell, NbFracCell

    divEgK(:) = 0.d0
    divEgS(:) = 0.d0
    divEgR(:,:) = 0.d0
    SmEg = 0.d0

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    ! upwind is cell k, DensiteMolaireKrViscoEnthalpie(k)*DarcyFlux
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       if(FluxDarcyKI(mph,s,k)>=0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             divEgK(j) = divEgK(j) &
                  + divDensiteMolaireKrViscoEnthalpieCell(j,m,k) * FluxDarcyKI(mph,s,k) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
          end do

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
             divEgS(j) = divEgS(j) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
          end do

          ! divEgR
          do r=1, NbNodeCell ! divR for r is node in dof(k)
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) =  divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieCell(m,k)
             end do
          end do

          do r=1, NbFracCell ! divR for r is frac in dof(k)

             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divEgR(j,rf) =  divEgR(j,rf) &
                     + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoEnthalpieCell(m,k)
             end do
          end do

          ! Sm
          SmEg = SmEg &
               + SmDensiteMolaireKrViscoEnthalpieCell(m,k) * FluxDarcyKI(mph,s,k) &
               + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
       end if
    end do

    ! upwind is node s
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyKI(mph,s,k)<0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
             divEgK(j) = divEgK(j) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end do

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
             divEgS(j) = divEgS(j) &
                  + divDensiteMolaireKrViscoEnthalpieNode(j,m,nums) * FluxDarcyKI(mph,s,k) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end do

          ! divEgR, r is node in dof(k)
          do r=1, NbNodeCell
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) = divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
             end do
          end do

          ! divEgR, r is frac in dof(k)
          do r=1, NbFracCell
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divEgR(j,rf) = divEgR(j,rf) &
                     + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
             end do
          end do

          ! Sm
          ! if nums is Dirichlet, Sm is supposed to be null
          if(IdNodeLocal(nums)%T /= "d") then
             SmEg = SmEg &
                  + SmDensiteMolaireKrViscoEnthalpieNode(m,nums) * FluxDarcyKI(mph,s,k) &
                  + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end if

       end if
    end do


  end subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellnode



  subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellfrac(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divEgK, divEgS, divEgR, SmEg)

    ! nums: frac num

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), &
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divEgK( NbIncTotalPrimMax), &
         divEgS( NbIncTotalPrimMax), &
         divEgR( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), &
         SmEg

    ! tmp
    integer :: m, mph, j, r, numr, rf, sf
    integer :: NbNodeCell, NbFracCell

    divEgK(:) = 0.d0
    divEgS(:) = 0.d0
    divEgR(:,:) = 0.d0
    SmEg = 0.d0

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    sf = s + NbNodeCell

    ! upwind is cell k, DensiteMolaireKrViscoEnthalpie(k)*DarcyFlux
    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       if(FluxDarcyKI(mph,sf,k)>=0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             divEgK(j) = divEgK(j) &
                  + divDensiteMolaireKrViscoEnthalpieCell(j,m,k) * FluxDarcyKI(mph,sf,k) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
          end do

          ! if(commRank==1) then
          !    print*, k, s, nums, m, j, &
          !         divDarcyFlux_k(j,mph)*DensiteMolaireKrViscoEnthalpieCell(m,k)
          ! end if

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divS
             divEgS(j) = divEgS(j) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
          end do

          ! divEgR
          do r=1, NbNodeCell ! divR for r is node in dof(k)
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) =  divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieCell(m,k)
             end do
          end do

          do r=1, NbFracCell ! divR for r is frac in dof(k)

             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divEgR(j,rf) =  divEgR(j,rf) &
                     + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoEnthalpieCell(m,k)
             end do
          end do

          ! Sm
          SmEg = SmEg &
               + SmDensiteMolaireKrViscoEnthalpieCell(m,k) * FluxDarcyKI(mph,sf,k) &
               + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieCell(m,k)
       end if
    end do

    ! upwind is frac s
    do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

       if(FluxDarcyKI(mph,sf,k)<0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divK
             divEgK(j) = divEgK(j) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieFrac(m,nums)
          end do

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
             divEgS(j) = divEgS(j) &
                  + divDensiteMolaireKrViscoEnthalpieFrac(j,m,nums) * FluxDarcyKI(mph,sf,k) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieFrac(m,nums)
          end do

          ! divEgR, r is node in dof(k)
          do r=1, NbNodeCell
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) = divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieFrac(m,nums)
             end do
          end do

          ! divEgR, r is frac in dof(k)
          do r=1, NbFracCell
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num
             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divEgR(j,rf) = divEgR(j,rf) &
                     + divDarcyFlux_r(j,mph,rf) * DensiteMolaireKrViscoEnthalpieFrac(m,nums)
             end do
          end do

          ! Sm
          SmEg = SmEg &
               + SmDensiteMolaireKrViscoEnthalpieFrac(m,nums) * FluxDarcyKI(mph,sf,k) &
               + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieFrac(m,nums)
       end if
    end do


  end subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_cellfrac


  subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_fracnode(k,s,nums, &
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, SmDarcyFlux, &
       divEgK, divEgS, divEgR, SmEg)

    integer, intent(in) :: k, s, nums

    double precision, intent(in) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeFaceMax), & ! r represent nodes of dof(k), k is frac
         SmDarcyFlux( NbPhase)

    double precision, intent(out) :: &
         divEgK( NbIncTotalPrimMax), &
         divEgS( NbIncTotalPrimMax), &
         divEgR( NbIncTotalPrimMax, NbNodeFaceMax), &
         SmEg

    ! tmp
    integer :: m, mph, j, r, numr, fk
    integer :: NbNodeFrac

    divEgK(:) = 0.d0
    divEgS(:) = 0.d0
    divEgR(:,:) = 0.d0
    SmEg = 0.d0

    fk = FracToFaceLocal(k) ! fk is num face

    ! number of nodes in frac k
    NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

    ! upwind is frac k, DensiteMolaireKrViscoEnthalpieFrac(k)*DarcyFlux
    do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

       ! upwind is frac
       if(FluxDarcyFI(mph,s,k)>=0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
             divEgK(j) = divEgK(j) &
                  + divDensiteMolaireKrViscoEnthalpieFrac(j,m,k) * FluxDarcyFI(mph,s,k) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieFrac(m,k)
          end do

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divS
             divEgS(j) = divEgS(j) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieFrac(m,k)
          end do

          ! divEgR
          do r=1, NbNodeFrac ! divR for r is node in dof(k)
             numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) =  divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieFrac(m,k)
             end do
          end do

          ! Sm
          SmEg = SmEg &
               + SmDensiteMolaireKrViscoEnthalpieFrac(m,k) * FluxDarcyFI(mph,s,k) &
               + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieFrac(m,k)
       end if
    end do

    ! upwind is node s
    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if(FluxDarcyFI(mph,s,k)<0.d0) then

          ! divEgK
          do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divK
             divEgK(j) = divEgK(j) &
                  + divDarcyFlux_k(j,mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end do

          ! divEgS
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
             divEgS(j) = divEgS(j) &
                  + divDensiteMolaireKrViscoEnthalpieNode(j,m,nums) * FluxDarcyFI(mph,s,k) &
                  + divDarcyFlux_s(j,mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end do

          ! divEgR, r is node in dof(k)
          do r=1, NbNodeFrac
             numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divEgR(j,r) = divEgR(j,r) &
                     + divDarcyFlux_r(j,mph,r) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
             end do
          end do

          ! Sm
          if(IdNodeLocal(nums)%T /= "d") then
             SmEg = SmEg &
                  + SmDensiteMolaireKrViscoEnthalpieNode(m,nums) * FluxDarcyFI(mph,s,k) &
                  + SmDarcyFlux(mph) * DensiteMolaireKrViscoEnthalpieNode(m,nums)
          end if

       end if
    end do

  end subroutine Jacobian_divDensiteMolaireKrViscoEnthalpieDarcyFlux_fracnode



  ! div: div prim
  ! div (V_{k,s}^alpha ) = divDarcyFlux_k * div(X_k)
  !                      + divDarcyFlux_s * div(X_s)
  !                      + \sum_{r \in V_k} divDarcyFlux_r * div(X_r)
  !                      + SmDarcyFlux
  subroutine Jacobian_divDarcyFlux_cellnode(k,s,nums,&
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
       SmDarcyFlux)

    ! k: cell num
    ! s: node num in cell k
    ! nums: node num in mesh

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), & ! r represent s' in paper
         SmDarcyFlux( NbPhase)

    double precision :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    integer :: NbNodeCell, NbFracCell
    logical :: Id_Qks(NbPhase)

    integer :: r, numr, rf, j, m, mph
    double precision :: sum_aks, sum_aksgz

    divDarcyFlux_r(:,:,:) = 0.d0
    SmDarcyFlux(:) = 0.d0

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    ! sum_aks = \sum_r a_{k,s}^r
    ! sum_aksgz = \sum_r a_{k,s}^r * g * (z_k-z_r)
    sum_aks = 0.d0
    sum_aksgz = 0.d0

    do r=1, NbNodeCell ! r is node
       numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r) ! numr is face num here

       sum_aks = sum_aks + TkLocal_Darcy(k)%pt(s,r)
       sum_aksgz = sum_aksgz + TkLocal_Darcy(k)%pt(s,r) * (XCellLocal(3,k)-XNodeLocal(3,numr))
    end do

    do r=1, NbFracCell ! r is frac
       numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r) ! numr is face num here
       rf = r + NbNodeCell

       sum_aks = sum_aks + TkLocal_Darcy(k)%pt(s,rf)
       sum_aksgz = sum_aksgz + TkLocal_Darcy(k)%pt(s,rf) * (XCellLocal(3,k)-XFaceLocal(3,numr))
    end do

    sum_aksgz = sum_aksgz * gravity

    ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
    call Jacobian_divrho_cellnode(k,s,nums, &
         divrho_k, Smrho_k, &
         divrho_s, Smrho_s)

    Id_Qks(:) = .false.

    ! if( k==1 .and. s==1 .and. commRank==1) then
    !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)

    !       print*, sum_aksgz/20.d0, divrho_k(j,1), divrho_k(j,2)

    !       ! print*, IncNode(nums)%Saturation(2)
    !       ! print*, IncCell(k)%Saturation
    !       ! print*, divrho_s(j,1)!*sum_aksgz!*DensiteMolaireKrViscoCompCell(1,1,1)
    !       ! print*, sum_aksgz*DensiteMolaireKrViscoCompCell(1,1,1)
    !       ! print*, sum_aksgz, DensiteMolaireKrViscoCompCell(1,1,1)
    !       ! print*, divDensiteMassiqueNode(j,1,nums)
    !       ! print*, divSaturationNode(j,1,nums)
    !       ! print*, DensiteMassiqueCell(1,k)
    !       ! print*, IncNode(nums)%Saturation(2)
    !    end do
    ! end if

    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       Id_Qks(mph) = .true.

       ! divDarcyFlux_k
       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
          divDarcyFlux_k(j,mph) = &
               sum_aks * divPressionCell(j,k) &        ! \sum a_{ks}^{s'} P_k
               + sum_aks * divPressionCapCell(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
               + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from k
       SmDarcyFlux(mph) = SmDarcyFlux(mph)&
            + sum_aks * SmPressionCell(k) &
            + sum_aksgz * Smrho_k(mph)    ! SmPressionCap=0

       ! divDarcyFlux_s
       do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
          divDarcyFlux_s(j,mph) = &
               sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from s
       SmDarcyFlux(mph) = SmDarcyFlux(mph) &
            + sum_aksgz * Smrho_s(mph)

       ! divDarcyFlux_r, r represent s' in paper, r is node
       do r=1, NbNodeCell
          numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

          ! P(mph) = PressionNode + PressionCap(mph)
          do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
             divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                  - TkLocal_Darcy(k)%pt(s,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                  - TkLocal_Darcy(k)%pt(s,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
          end do

          ! SmDarcyFlux from r (node)
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               - TkLocal_Darcy(k)%pt(s,r) * SmPressionNode(numr)    ! -a_{ks}^{s'} * Sm

       end do ! end of r node

       ! divDarcyFlux_r, r represent s' in paper, r is frac
       do r=1, NbFracCell
          numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here

          rf = r + NbNodeCell

          do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
             divDarcyFlux_r(j,mph,rf) = divDarcyFlux_r(j,mph,rf) &
                  - TkLocal_Darcy(k)%pt(s,rf) * divPressionFrac(j,numr) &        ! a_{ks}^{s'} -P_s'
                  - TkLocal_Darcy(k)%pt(s,rf) * divPressionCapFrac(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
          end do

          ! SmDarcyFlux from r (frac)
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               - TkLocal_Darcy(k)%pt(s,rf) * SmPressionFrac(numr)    ! -a_{ks}^{s'} * Sm

       end do ! end of r frac
    end do

    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if( Id_Qks(mph) .eqv. .false.) then ! this phase is not in Q_k

          ! divDarcyFlux_k
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             divDarcyFlux_k(j,mph) = &
                  sum_aks * divPressionCell(j,k) &        ! \sum a_{ks}^{s'} P_k
                  + sum_aks * divPressionCapCell(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
                  + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from k
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aks * SmPressionCell(k) &
               + sum_aksgz * Smrho_k(mph)   ! SmPressionCap=0

          ! divDarcyFlux_s
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
             divDarcyFlux_s(j,mph) = &
                  sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from s
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aksgz * Smrho_s(mph)

          ! divDarcyFlux_r, r represent s' in paper, r is node
          do r=1, NbNodeCell
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r) ! numr is frac num here

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                     - TkLocal_Darcy(k)%pt(s,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                     - TkLocal_Darcy(k)%pt(s,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
             end do

             ! SmDarcyFlux_r from r (node)
             SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                  - TkLocal_Darcy(k)%pt(s,r) * SmPressionNode(numr) ! -a_{ks}^{s'} * Sm

          end do ! end of r node

          ! divDarcyFlux_r, r represent s' in paper, r is frac
          do r=1, NbFracCell
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here

             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divDarcyFlux_r(j,mph,rf) = divDarcyFlux_r(j,mph,rf) &
                     - TkLocal_Darcy(k)%pt(s,rf) * divPressionFrac(j,numr) &        ! a_{ks}^{s'} -P_s'
                     - TkLocal_Darcy(k)%pt(s,rf) * divPressionCapFrac(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
             end do

             ! SmDarcyFlux from r (frac)
             SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                  - TkLocal_Darcy(k)%pt(s,rf) * SmPressionFrac(numr) ! -a_{ks}^{s'} * Sm

          end do ! end of r frac

       end if ! end of Id_Qki(mph)
    end do ! end of Q_s

    ! if( k==1 .and. s==1 .and. commRank==1) then
    !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
    !       print*, sum_aksgz, divDarcyFlux_k(j,1), divDarcyFlux_k(j,2)
    !    end do
    ! end if

  end subroutine Jacobian_divDarcyFlux_cellnode



  ! div: div prim
  ! div (V_{k,s}^alpha ) = divDarcyFlux_k * div(X_k)
  !                      + divDarcyFlux_s * div(X_s)
  !                      + \sum_r divDarcyFlux_r * div(X_r)
  !                      + SmDarcyFlux
  subroutine Jacobian_divDarcyFlux_cellfrac(k,s,nums,&
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
       SmDarcyFlux)

    ! k: cell num
    ! s: frac num in cell k
    ! nums: frac num

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeCellMax+NbFracCellMax), & ! r represent s' in paper
         SmDarcyFlux( NbPhase)

    double precision :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    integer :: NbNodeCell, NbFracCell
    logical :: Id_Qks(NbPhase)

    integer :: r, numr, j, m, mph, sf, rf
    double precision :: sum_aks, sum_aksgz

    divDarcyFlux_r(:,:,:) = 0.d0
    SmDarcyFlux(:) = 0.d0

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    sf = s + NbNodeCell ! TkLocal_Darcy(k)%pt(sf,s') for frac s

    ! sum_aks = \sum_r a_{k,s}^r
    ! sum_aksgz = \sum_r a_{k,s}^r * g * (z_k-z_r)
    sum_aks = 0.d0
    sum_aksgz = 0.d0

    do r=1, NbNodeCell ! r is node
       numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r) ! numr is face num here

       sum_aks = sum_aks + TkLocal_Darcy(k)%pt(sf,r)
       sum_aksgz = sum_aksgz + TkLocal_Darcy(k)%pt(sf,r) * (XCellLocal(3,k)-XNodeLocal(3,numr))
    end do

    do r=1, NbFracCell ! r is frac
       numr = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r) ! numr is face num here
       rf = r + NbNodeCell

       sum_aks = sum_aks + TkLocal_Darcy(k)%pt(sf,rf)
       sum_aksgz = sum_aksgz + TkLocal_Darcy(k)%pt(sf,rf) * (XCellLocal(3,k)-XFaceLocal(3,numr))
    end do

    sum_aksgz = sum_aksgz * gravity

    ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
    call Jacobian_divrho_cellfrac(k,s,nums, &
         divrho_k, Smrho_k, &
         divrho_s, Smrho_s)

    Id_Qks(:) = .false.

    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       Id_Qks(mph) = .true.

       ! divDarcyFlux_k
       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
          divDarcyFlux_k(j,mph) = &
               sum_aks * divPressionCell(j,k) &        ! \sum a_{ks}^{s'} P_k
               + sum_aks * divPressionCapCell(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
               + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from k
       SmDarcyFlux(mph) = SmDarcyFlux(mph) &
            + sum_aks * SmPressionCell(k) &
            + sum_aksgz * Smrho_k(mph)    ! SmPressionCap=0

       ! divDarcyFlux_s
       do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
          divDarcyFlux_s(j,mph) = &
               sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from s
       SmDarcyFlux(mph) = SmDarcyFlux(mph) &
            + sum_aksgz * Smrho_s(mph)

       ! divDarcyFlux_r, r represent s' in paper, r is node
       do r=1, NbNodeCell
          numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

          do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
             divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                  - TkLocal_Darcy(k)%pt(sf,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                  - TkLocal_Darcy(k)%pt(sf,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
          end do

          ! SmDarcyFlux from r (node)
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               - TkLocal_Darcy(k)%pt(sf,r) * SmPressionNode(numr)    ! -a_{ks}^{s'} * Sm

       end do ! end of r node

       ! divDarcyFlux_r, r represent s' in paper, r is frac
       do r=1, NbFracCell
          numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here

          rf = r + NbNodeCell

          do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
             divDarcyFlux_r(j,mph,rf) = divDarcyFlux_r(j,mph,rf) &
                  - TkLocal_Darcy(k)%pt(sf,rf) * divPressionFrac(j,numr) &        ! a_{ks}^{s'} -P_s'
                  - TkLocal_Darcy(k)%pt(sf,rf) * divPressionCapFrac(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
          end do

          ! SmDarcyFlux from r (frac)
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               - TkLocal_Darcy(k)%pt(sf,rf) * SmPressionFrac(numr)    ! -a_{ks}^{s'} * Sm

       end do ! end of r frac

    end do

    do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

       if( Id_Qks(mph) .eqv. .false.) then ! this phase is not in Q_k

          ! divDarcyFlux_k
          do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
             divDarcyFlux_k(j,mph) = &
                  sum_aks * divPressionCell(j,k) &        ! \sum a_{ks}^{s'} P_k
                  + sum_aks * divPressionCapCell(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
                  + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from k
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aks * SmPressionCell(k)  &
               + sum_aksgz * Smrho_k(mph)   ! SmPressionCap=0

          ! divDarcyFlux_s
          do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
             divDarcyFlux_s(j,mph) = &
                  sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from s
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aksgz * Smrho_s(mph)

          ! divDarcyFlux_r, r represent s' in paper, r is node
          do r=1, NbNodeCell
             numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r) ! numr is frac num here

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                     - TkLocal_Darcy(k)%pt(sf,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                     - TkLocal_Darcy(k)%pt(sf,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
             end do

             ! SmDarcyFlux from r (node)
             SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                  - TkLocal_Darcy(k)%pt(sf,r) * SmPressionNode(numr) ! -a_{ks}^{s'} * Sm

          end do ! end of r node

          ! divDarcyFlux_r, r represent s' in paper, r is frac
          do r=1, NbFracCell
             numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here

             rf = r + NbNodeCell

             do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
                divDarcyFlux_r(j,mph,rf) = divDarcyFlux_r(j,mph,rf) &
                     - TkLocal_Darcy(k)%pt(sf,rf) * divPressionFrac(j,numr) &        ! a_{ks}^{s'} -P_s'
                     - TkLocal_Darcy(k)%pt(sf,rf) * divPressionCapFrac(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
             end do

             ! SmDarcyFlux from r (frac)
             SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                  - TkLocal_Darcy(k)%pt(sf,rf) * SmPressionFrac(numr) ! -a_{ks}^{s'} * Sm

          end do ! end of r frac

       end if ! end of Id_Qki(mph)
    end do ! end of Q_s

    ! if( k==1 .and. s==1 .and. commRank==1) then
    !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
    !       print*, sum_aks, divDarcyFlux_k(j,1), divDarcyFlux_k(j,2)
    !    end do
    ! end if

  end subroutine Jacobian_divDarcyFlux_cellfrac



  ! div: div prim
  ! div (V_{k,s}^alpha ) = divDarcyFlux_k * div(X_k)
  !                      + divDarcyFlux_s * div(X_s)
  !                      + \sum_r divDarcyFlux_r * div(X_r)
  !                      + SmDarcyFlux
  subroutine Jacobian_divDarcyFlux_fracnode(k,s,nums,&
       divDarcyFlux_k, divDarcyFlux_s, divDarcyFlux_r, &
       SmDarcyFlux)

    ! k: frac num
    ! s: node num in frac k
    ! nums: node num in mesh

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divDarcyFlux_k( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_s( NbIncTotalPrimMax, NbPhase), &
         divDarcyFlux_r( NbIncTotalPrimMax, NbPhase, NbNodeFaceMax), & ! r represent s' in paper                                !
         SmDarcyFlux( NbPhase)

    double precision :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    integer :: NbNodeFrac
    logical :: Id_Qks(NbPhase)

    integer :: r, numr, j, m, mph, fk
    double precision :: sum_aks, sum_aksgz

    divDarcyFlux_r(:,:,:) = 0.d0
    SmDarcyFlux(:) = 0.d0

    fk = FracToFaceLocal(k) ! fk is num face

    ! number of nodes in frac k
    NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

    ! sum_aks = \sum_r a_{k,s}^r
    ! sum_aksgz = \sum_r a_{k,s}^r * g * (z_k-z_r)
    sum_aks = 0.d0
    sum_aksgz = 0.d0

    do r=1, NbNodeFrac ! r is node
       numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

       sum_aks = sum_aks + TkFracLocal_Darcy(k)%pt(s,r)
       sum_aksgz = sum_aksgz + TkFracLocal_Darcy(k)%pt(s,r) * (XFaceLocal(3,fk)-XNodeLocal(3,numr))
    end do

    sum_aksgz = sum_aksgz * gravity

    ! div ( rho_{k,s}^alpha ) for all phase Q_k \cup Q_s
    call Jacobian_divrho_fracnode(k,s,nums, &
         divrho_k, Smrho_k, &
         divrho_s, Smrho_s)

    Id_Qks(:) = .false.

    do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k, k is frac
       mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

       Id_Qks(mph) = .true.

       ! divDarcyFlux_k
       do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
          divDarcyFlux_k(j,mph) = &
               sum_aks * divPressionFrac(j,k) &        ! \sum a_{ks}^{s'} P_k
               + sum_aks * divPressionCapFrac(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
               + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from k
       SmDarcyFlux(mph) = SmDarcyFlux(mph) &
            + sum_aks * SmPressionFrac(k) &
            + sum_aksgz * Smrho_k(mph)    ! SmPressionCap=0

       ! divDarcyFlux_s
       do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
          divDarcyFlux_s(j,mph) = &
               sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
       end do

       ! SmDarcyFlux from s
       SmDarcyFlux(mph) = SmDarcyFlux(mph) &
            + sum_aksgz * Smrho_s(mph)

       ! divDarcyFlux_r, r represent s' in paper, r is node
       do r=1, NbNodeFrac
          numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

          do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
             divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                  - TkFracLocal_Darcy(k)%pt(s,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                  - TkFracLocal_Darcy(k)%pt(s,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
          end do

          ! SmDarcyFlux from r
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               - TkFracLocal_Darcy(k)%pt(s,r) * SmPressionNode(numr)    ! -a_{ks}^{s'} * Sm

       end do ! end of r node

    end do

    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       if( Id_Qks(mph) .eqv. .false.) then ! this phase is not in Q_k

          ! divDarcyFlux_k
          do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
             divDarcyFlux_k(j,mph) = &
                  sum_aks * divPressionFrac(j,k) &        ! \sum a_{ks}^{s'} P_k
                  + sum_aks * divPressionCapFrac(j,mph,k) & ! \sum a_{ks}^{s'} PressionCap_k
                  + sum_aksgz * divrho_k(j,mph)           ! \sum a_{ks}^{s'} divrho_k * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from k
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aks * SmPressionFrac(k) &
               + sum_aksgz * Smrho_k(mph)   ! SmPressionCap=0

          ! divDarcyFlux_s
          do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
             divDarcyFlux_s(j,mph) = &
                  sum_aksgz * divrho_s(j,mph)             ! \sum a_{ks}^{s'} divrho_s * g * (z_k-z_s')
          end do

          ! SmDarcyFlux from s
          SmDarcyFlux(mph) = SmDarcyFlux(mph) &
               + sum_aksgz * Smrho_s(mph)

          ! divDarcyFlux_r, r represent s' in paper, r is node
          do r=1, NbNodeFrac
             numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r) ! numr is frac num here

             do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
                divDarcyFlux_r(j,mph,r) = divDarcyFlux_r(j,mph,r) &
                     - TkFracLocal_Darcy(k)%pt(s,r) * divPressionNode(j,numr) &        ! a_{ks}^{s'} -P_s'
                     - TkFracLocal_Darcy(k)%pt(s,r) * divPressionCapNode(j,mph,numr)     ! a_{ks}^{s'} -PressionCap_s'
             end do

             ! SmDarcyFlux from r (node)
             SmDarcyFlux(mph) = SmDarcyFlux(mph) &
                  - TkFracLocal_Darcy(k)%pt(s,r) * SmPressionNode(numr) ! -a_{ks}^{s'} * Sm

          end do ! end of r node

       end if ! end of Id_Qki(mph)
    end do ! end of Q_s

  end subroutine Jacobian_divDarcyFlux_fracnode



  ! div ( rho_{k,s}^alpha )
  !   = divrho_k * div(X_k) + divrho_s * div(X_s)
  !   + Smrho_k + Smrho_s
  subroutine Jacobian_divrho_cellnode(k,s,nums, &
       divrho_k, Smrho_k, &
       divrho_s, Smrho_s)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    ! double precision :: Satki
    integer :: j, m, mph, tmp_compt(NbPhase)

    divrho_k(:,:) = 0.d0
    divrho_s(:,:) = 0.d0
    Smrho_k(:) = 0.d0
    Smrho_s(:) = 0.d0
    tmp_compt(:) = 0

    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       ! Satki = IncCell(k)%Saturation(mph) + IncNode(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha
       tmp_compt(mph) = tmp_compt(mph) + 1

       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
          divrho_k(j,mph) = divDensiteMassiqueCell(j,mph,k)
       end do

       Smrho_k(mph) = SmDensiteMassiqueCell(mph,k)

       ! if( abs(Satki)<eps) then ! Satki == 0

       !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
       !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueCell(j,mph,k)
       !    end do

       !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueCell(mph,k)

       !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
       !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueNode(j,mph,nums)
       !    end do

       !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueNode(mph,nums)

       ! else

       !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
       !       divrho_k(j,mph) = &
       !            divSaturationCell(j,m,k) / Satki * DensiteMassiqueCell(mph,k) &
       !            - divSaturationCell(j,m,k) / (Satki**2) &
       !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
       !            + IncCell(k)%Saturation(mph) / Satki * divDensiteMassiqueCell(j,mph,k) &
       !                          !
       !            - divSaturationCell(j,m,k) / (Satki**2) &
       !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums)
       !    end do

       !    Smrho_k(mph) = &
       !         IncCell(k)%Saturation(mph) / Satki * SmDensiteMassiqueCell(mph,k) ! SmSaturation=0

       !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
       !       divrho_s(j,mph) = &
       !            - divSaturationNode(j,m,nums) / (Satki**2) &
       !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
       !            + divSaturationNode(j,m,nums) / Satki * DensiteMassiqueNode(mph,nums) &
       !            - divSaturationNode(j,m,nums) / (Satki**2) &
       !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums) &
       !            + IncNode(nums)%Saturation(mph) / Satki * divDensiteMassiqueNode(j,mph,nums)
       !    end do

       !    Smrho_s(mph) = &
       !         IncNode(nums)%Saturation(mph) / Satki * SmDensiteMassiqueNode(mph,nums)

       ! end if ! end of Id_Qks(mph)

    end do ! end of Q_k


    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       tmp_compt(mph) = tmp_compt(mph) + 1

          ! Satki = IncCell(k)%Saturation(mph) + IncNode(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha

       do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
          divrho_s(j,mph) = divDensiteMassiqueNode(j,mph,nums)
       end do

       Smrho_s(mph) = SmDensiteMassiqueNode(mph,nums)

          ! if( abs(Satki)<eps) then ! Satki == 0

          !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
          !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueCell(j,mph,k)
          !    end do

          !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueCell(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
          !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueNode(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueNode(mph,nums)

          ! else

          !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
          !       divrho_k(j,mph) = &
          !            + divSaturationCell(j,m,k) / Satki * DensiteMassiqueCell(mph,k) &
          !            - divSaturationCell(j,m,k) / (Satki**2) &
          !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
          !            + IncCell(k)%Saturation(mph) / Satki * divDensiteMassiqueCell(j,mph,k) &
          !                       !
          !            - divSaturationCell(j,m,k) / (Satki**2) &
          !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums)
          !    end do

          !    Smrho_k(mph) = &
          !         IncCell(k)%Saturation(mph) / Satki * SmDensiteMassiqueCell(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
          !       divrho_s(j,mph) = &
          !            - divSaturationNode(j,m,nums) / (Satki**2) &
          !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
          !            + divSaturationNode(j,m,nums) / Satki * DensiteMassiqueNode(mph,nums) &
          !            - divSaturationNode(j,m,nums) / (Satki**2) &
          !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums) &
          !            + IncNode(nums)%Saturation(mph) / Satki * divDensiteMassiqueNode(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = &
          !         IncNode(nums)%Saturation(mph) / Satki * SmDensiteMassiqueNode(mph,nums)

          ! end if ! end of Id_Qks(mph)

    end do

    do m=1, NbPhase
       divrho_k(:,m) = divrho_k(:,m) / max(tmp_compt(m), 1)
       divrho_s(:,m) = divrho_s(:,m) / max(tmp_compt(m), 1)
       Smrho_k(m) = Smrho_k(m) / max(tmp_compt(m), 1)
       Smrho_s(m) = Smrho_s(m) / max(tmp_compt(m), 1)
    enddo


    ! if(k==1 .and. s==1 .and. commRank==0) then
    !    print*, "ph 1", divrho_s(:,1)
    !    print*, "ph 2", divrho_s(:,2)
    ! end if

  end subroutine Jacobian_divrho_cellnode


  ! div ( rho_{k,s}^alpha )
  !   = divrho_k * div(X_k) + divrho_s * div(X_s)
  !   + Smrho_k + Smrho_s
  subroutine Jacobian_divrho_cellfrac(k,s,nums, &
       divrho_k, Smrho_k, &
       divrho_s, Smrho_s)

    ! k: cell num
    ! s: frac num in cell k
    ! nums: frac num in mesh

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    ! double precision :: Satki
    integer :: j, m, mph, tmp_compt(NbPhase)

    divrho_k(:,:) = 0.d0
    divrho_s(:,:) = 0.d0
    Smrho_k(:) = 0.d0
    Smrho_s(:) = 0.d0
    tmp_compt = 0

    do m=1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncCell(k)%ic)

       !Satki = IncCell(k)%Saturation(mph) + IncFrac(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha
       tmp_compt(mph) = tmp_compt(mph) + 1

       do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
          divrho_k(j,mph) = divDensiteMassiqueCell(j,mph,k)
       end do

       Smrho_k(mph) = SmDensiteMassiqueCell(mph,k)

       ! if( abs(Satki)<eps) then ! Satki == 0

       !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
       !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueCell(j,mph,k)
       !    end do

       !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueCell(mph,k)

       !    do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divrho_s
       !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueFrac(j,mph,nums)
       !    end do

       !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueFrac(mph,nums)

       ! else

       !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
       !       divrho_k(j,mph) = &
       !            divSaturationCell(j,m,k) / Satki * DensiteMassiqueCell(mph,k) &
       !            - divSaturationCell(j,m,k) / (Satki**2) &
       !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
       !            + IncCell(k)%Saturation(mph) / Satki * divDensiteMassiqueCell(j,mph,k) &
       !                          !
       !            - divSaturationCell(j,m,k) / (Satki**2) &
       !            * IncFrac(nums)%Saturation(mph) * DensiteMassiqueFrac(mph,nums)
       !    end do

       !    Smrho_k(mph) = &
       !         IncCell(k)%Saturation(mph) / Satki * SmDensiteMassiqueCell(mph,k) ! SmSaturation=0

       !    do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
       !       divrho_s(j,mph) = &
       !            - divSaturationFrac(j,m,nums) / (Satki**2) &
       !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
       !            + divSaturationFrac(j,m,nums) / Satki * DensiteMassiqueFrac(mph,nums) &
       !            - divSaturationFrac(j,m,nums) / (Satki**2) &
       !            * IncFrac(nums)%Saturation(mph) * DensiteMassiqueFrac(mph,nums) &
       !            + IncFrac(nums)%Saturation(mph) / Satki * divDensiteMassiqueFrac(j,mph,nums)
       !    end do

       !    Smrho_s(mph) = &
       !         IncFrac(nums)%Saturation(mph) / Satki * SmDensiteMassiqueFrac(mph,nums)

       ! end if ! end of Id_Qks(mph)

    end do ! end of Q_k


    do m=1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

       ! Satki = IncCell(k)%Saturation(mph) + IncFrac(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha
       tmp_compt(mph) = tmp_compt(mph) + 1

       do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divrho_s
          divrho_s(j,mph) = divDensiteMassiqueFrac(j,mph,nums)
       end do

       Smrho_s(mph) = SmDensiteMassiqueFrac(mph,nums)

          ! if( abs(Satki)<eps) then ! Satki == 0

          !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic) ! divrho_k
          !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueCell(j,mph,k)
          !    end do

          !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueCell(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic) ! divrho_s
          !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueFrac(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueFrac(mph,nums)

          ! else

          !    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
          !       divrho_k(j,mph) = &
          !            + divSaturationCell(j,m,k) / Satki * DensiteMassiqueCell(mph,k) &
          !            - divSaturationCell(j,m,k) / (Satki**2) &
          !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
          !            + IncCell(k)%Saturation(mph) / Satki * divDensiteMassiqueCell(j,mph,k) &
          !                       !
          !            - divSaturationCell(j,m,k) / (Satki**2) &
          !            * IncFrac(nums)%Saturation(mph) * DensiteMassiqueFrac(mph,nums)
          !    end do

          !    Smrho_k(mph) = &
          !         IncCell(k)%Saturation(mph) / Satki * SmDensiteMassiqueCell(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncFrac(nums)%ic)
          !       divrho_s(j,mph) = &
          !            - divSaturationFrac(j,m,nums) / (Satki**2) &
          !            * IncCell(k)%Saturation(mph) * DensiteMassiqueCell(mph,k) &
          !            + divSaturationFrac(j,m,nums) / Satki * DensiteMassiqueFrac(mph,nums) &
          !            - divSaturationFrac(j,m,nums) / (Satki**2) &
          !            * IncFrac(nums)%Saturation(mph) * DensiteMassiqueFrac(mph,nums) &
          !            + IncFrac(nums)%Saturation(mph) / Satki * divDensiteMassiqueFrac(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = &
          !         IncFrac(nums)%Saturation(mph) / Satki * SmDensiteMassiqueFrac(mph,nums)

          ! end if ! end of Id_Qks(mph)

    end do

    do m=1, NbPhase
      divrho_k(:,m) = divrho_k(:,m) / max(tmp_compt(m), 1)
      divrho_s(:,m) = divrho_s(:,m) / max(tmp_compt(m), 1)
      Smrho_k(m) = Smrho_k(m) / max(tmp_compt(m), 1)
      Smrho_s(m) = Smrho_s(m) / max(tmp_compt(m), 1)
    enddo

    ! if(k==1 .and. s==1 .and. commRank==0) then
    !    print*, "ph 1", divrho_s(:,1)
    !    print*, "ph 2", divrho_s(:,2)
    ! end if

  end subroutine Jacobian_divrho_cellfrac



  ! div ( rho_{k,s}^alpha )
  !   = divrho_k * div(X_k) + divrho_s * div(X_s)
  !   + Smrho_k + Smrho_s
  subroutine Jacobian_divrho_fracnode(k, s, nums, &
       divrho_k, Smrho_k, &
       divrho_s, Smrho_s)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divrho_k( NbIncTotalPrimMax, NbPhase), &
         divrho_s( NbIncTotalPrimMax, NbPhase), &
         Smrho_k( NbPhase), &
         Smrho_s( NbPhase)

    ! double precision :: Satki
    integer :: j, m, mph, tmp_compt(NbPhase)

    divrho_k(:,:) = 0.d0
    divrho_s(:,:) = 0.d0
    Smrho_k(:) = 0.d0
    Smrho_s(:) = 0.d0
    tmp_compt(:) = 0


    do m=1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
       mph = NumPhasePresente_ctx(m,IncFrac(k)%ic)

       ! Satki = IncFrac(k)%Saturation(mph) + IncNode(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha
       tmp_compt(mph) = tmp_compt(mph) + 1

       do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divrho_k
          divrho_k(j,mph) = divDensiteMassiqueFrac(j,mph,k)
       end do

       Smrho_k(mph) = SmDensiteMassiqueFrac(mph,k)

       ! if( abs(Satki)<eps) then ! Satki == 0

       !    do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divrho_k
       !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueFrac(j,mph,k)
       !    end do

       !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueFrac(mph,k)

       !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
       !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueNode(j,mph,nums)
       !    end do

       !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueNode(mph,nums)

       ! else

       !    do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
       !       divrho_k(j,mph) = &
       !            divSaturationFrac(j,m,k) / Satki * DensiteMassiqueFrac(mph,k) &
       !            - divSaturationFrac(j,m,k) / (Satki**2) &
       !            * IncFrac(k)%Saturation(mph) * DensiteMassiqueFrac(mph,k) &
       !            + IncFrac(k)%Saturation(mph) / Satki * divDensiteMassiqueFrac(j,mph,k) &
       !                          !
       !            - divSaturationFrac(j,m,k) / (Satki**2) &
       !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums)
       !    end do

       !    Smrho_k(mph) = &
       !         IncFrac(k)%Saturation(mph) / Satki * SmDensiteMassiqueFrac(mph,k) ! SmSaturation=0

       !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
       !       divrho_s(j,mph) = &
       !            - divSaturationNode(j,m,nums) / (Satki**2) &
       !            * IncFrac(k)%Saturation(mph) * DensiteMassiqueFrac(mph,k) &
       !            + divSaturationNode(j,m,nums) / Satki * DensiteMassiqueNode(mph,nums) &
       !            - divSaturationNode(j,m,nums) / (Satki**2) &
       !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums) &
       !            + IncNode(nums)%Saturation(mph) / Satki * divDensiteMassiqueNode(j,mph,nums)

       !    end do

       !    Smrho_s(mph) = &
       !         IncNode(nums)%Saturation(mph) / Satki * SmDensiteMassiqueNode(mph,nums)

       ! end if ! end of Id_Qks(mph)

    end do ! end of Q_k


    do m=1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
       mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

       ! Satki = IncFrac(k)%Saturation(mph) + IncNode(nums)%Saturation(mph) ! S_k^alpha+S_i^alpha
       tmp_compt(mph) = tmp_compt(mph) + 1

       do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
          divrho_s(j,mph) = divDensiteMassiqueNode(j,mph,nums)
       end do

       Smrho_s(mph) = SmDensiteMassiqueNode(mph,nums)

          ! if( abs(Satki)<eps) then ! Satki == 0

          !    do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic) ! divrho_k
          !       divrho_k(j,mph) = 0.5d0 * divDensiteMassiqueFrac(j,mph,k)
          !    end do

          !    Smrho_k(mph) = 0.5d0 * SmDensiteMassiqueFrac(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic) ! divrho_s
          !       divrho_s(j,mph) = 0.5d0 * divDensiteMassiqueNode(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = 0.5d0 * SmDensiteMassiqueNode(mph,nums)

          ! else

          !    do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
          !       divrho_k(j,mph) = &
          !            + divSaturationFrac(j,m,k) / Satki * DensiteMassiqueFrac(mph,k) &
          !            - divSaturationFrac(j,m,k) / (Satki**2) &
          !            * IncFrac(k)%Saturation(mph) * DensiteMassiqueFrac(mph,k) &
          !            + IncFrac(k)%Saturation(mph) / Satki * divDensiteMassiqueFrac(j,mph,k) &
          !                       !
          !            - divSaturationFrac(j,m,k) / (Satki**2) &
          !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums)
          !    end do

          !    Smrho_k(mph) = &
          !         IncFrac(k)%Saturation(mph) / Satki * SmDensiteMassiqueFrac(mph,k)

          !    do j=1, NbIncTotalPrim_ctx(IncNode(nums)%ic)
          !       divrho_s(j,mph) = &
          !            - divSaturationNode(j,m,nums) / (Satki**2) &
          !            * IncFrac(k)%Saturation(mph) * DensiteMassiqueFrac(mph,k) &
          !            + divSaturationNode(j,m,nums) / Satki * DensiteMassiqueNode(mph,nums) &
          !            - divSaturationNode(j,m,nums) / (Satki**2) &
          !            * IncNode(nums)%Saturation(mph) * DensiteMassiqueNode(mph,nums) &
          !            + IncNode(nums)%Saturation(mph) / Satki * divDensiteMassiqueNode(j,mph,nums)
          !    end do

          !    Smrho_s(mph) = &
          !         IncNode(nums)%Saturation(mph) / Satki * SmDensiteMassiqueNode(mph,nums)


          ! end if ! end of Id_Qks(mph)

    end do

    do m=1, NbPhase
      divrho_k(:,m) = divrho_k(:,m) / max(tmp_compt(m), 1)
      divrho_s(:,m) = divrho_s(:,m) / max(tmp_compt(m), 1)
      Smrho_k(m) = Smrho_k(m) / max(tmp_compt(m), 1)
      Smrho_s(m) = Smrho_s(m) / max(tmp_compt(m), 1)
    enddo

    ! if(k==1 .and. s==1 .and. commRank==0) then
    !    print*, "ph 1", divrho_s(:,1)
    !    print*, "ph 2", divrho_s(:,2)
    ! end if

  end subroutine Jacobian_divrho_fracnode



  ! div Flux Fourier
  subroutine Jacobian_divFourierFlux_cellnode(k,s,nums, &
       divFourierFlux_k, & ! sum_{s'} a_{k,s}^s' T_k
       divFourierFlux_r, & ! sum_{s'} a_{k,s}^s' -T_s'
       SmFourierFlux )

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divFourierFlux_k( NbIncTotalPrimMax), &
         divFourierFlux_r( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), & ! r represent s' in paper
         SmFourierFlux

    ! tmp
    integer :: j, r, numr, rf
    integer :: NbNodeCell, NbFracCell

    double precision :: sum_aks

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    ! sum_aks = \sum_{r \in dof(k)} a_{k,s}^r
    sum_aks = 0.d0
    do r=1, NbNodeCell+NbFracCell
       sum_aks = sum_aks + TkLocal_Fourier(k)%pt(s,r)
    end do

    ! divFourierfFlux_k
    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
       divFourierFlux_k(j) = sum_aks * divTemperatureCell(j,k)
    end do

    SmFourierFlux = sum_aks * SmTemperatureCell(k)

    ! divFourierFlux_r, r represent s' in paper, r is node
    do r=1, NbNodeCell
       numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

       do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
          divFourierFlux_r(j,r) = - TkLocal_Fourier(k)%pt(s,r) * divTemperatureNode(j,numr)
       end do

       SmFourierFlux = SmFourierFlux &
            - TkLocal_Fourier(k)%pt(s,r) * SmTemperatureNode(numr)
    end do

    ! divFourierFlux_r, r represent s' in paper, r is frac
    do r=1, NbFracCell
       numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here
       rf = r + NbNodeCell

       do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
          divFourierFlux_r(j,rf) = - TkLocal_Fourier(k)%pt(s,rf) * divTemperatureFrac(j,numr)
       end do

       SmFourierFlux = SmFourierFlux &
            - TkLocal_Fourier(k)%pt(s,rf) * SmTemperatureFrac(numr)
    end do

  end subroutine Jacobian_divFourierFlux_cellnode


  ! div Flux Fourier
  subroutine Jacobian_divFourierFlux_cellfrac(k,s,nums, &
       divFourierFlux_k, & ! sum_{s'} a_{k,s}^s' T_k
       divFourierFlux_r, & ! sum_{s'} a_{k,s}^s' -T_s'
       SmFourierFlux)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divFourierFlux_k( NbIncTotalPrimMax), &
         divFourierFlux_r( NbIncTotalPrimMax, NbNodeCellMax+NbFracCellMax), & ! r represent s' in paper, dof(k)
         SmFourierFlux

    ! tmp
    integer :: j, r, numr, sf, rf
    integer :: NbNodeCell, NbFracCell

    double precision :: sum_aks

    ! number of nodes/fracs in cell k
    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

    sf = s + NbNodeCell

    ! sum_aks = \sum_{r \in dof(k)} a_{k,s}^r
    sum_aks = 0.d0
    do r=1, NbNodeCell+NbFracCell
       sum_aks = sum_aks + TkLocal_Fourier(k)%pt(sf,r)
    end do

    ! divFourierfFlux_k
    do j=1, NbIncTotalPrim_ctx(IncCell(k)%ic)
       divFourierFlux_k(j) = sum_aks * divTemperatureCell(j,k)
    end do

    SmFourierFlux = sum_aks * SmTemperatureCell(k)

    ! divFourierFlux_r, r is node
    do r=1, NbNodeCell
       numr = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+r)

       do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
          divFourierFlux_r(j,r) = - TkLocal_Fourier(k)%pt(sf,r) * divTemperatureNode(j,numr)
       end do

       SmFourierFlux = SmFourierFlux &
            - TkLocal_Fourier(k)%pt(sf,r) * SmTemperatureNode(numr)
    end do

    ! divFourierFlux_r, r is frac
    do r=1, NbFracCell
       numr = FaceToFracLocal( FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+r)) ! numr is frac num here
       rf = r + NbNodeCell

       do j=1, NbIncTotalPrim_ctx(IncFrac(numr)%ic)
          divFourierFlux_r(j,rf) = - TkLocal_Fourier(k)%pt(sf,rf) * divTemperatureFrac(j,numr)
       end do

       SmFourierFlux = SmFourierFlux &
            - TkLocal_Fourier(k)%pt(sf,rf) * SmTemperatureFrac(numr)
    end do

  end subroutine Jacobian_divFourierFlux_cellfrac


  ! div Flux Fourier
  subroutine Jacobian_divFourierFlux_fracnode(k,s,nums, &
       divFourierFlux_k, & ! sum_{s'} a_{k,s}^s' T_k
       divFourierFlux_r, & ! sum_{s'} a_{k,s}^s' -T_s'
       SmFourierFlux)

    integer, intent(in) :: k, s, nums

    double precision, intent(out) :: &
         divFourierFlux_k( NbIncTotalPrimMax), &
         divFourierFlux_r( NbIncTotalPrimMax, NbNodeFaceMax), & ! r represent s' in paper
         SmFourierFlux

    ! tmp
    integer :: j, r, numr, fk
    integer :: NbNodeFrac

    double precision :: sum_aks

    fk = FracToFaceLocal(k) ! fk is num face

    ! number of nodes in frac k
    NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

    ! sum_aks = \sum_r a_{k,s}^r
    sum_aks = 0.d0
    do r=1, NbNodeFrac ! r is node
       sum_aks = sum_aks + TkFracLocal_Fourier(k)%pt(s,r)
    end do

    ! divFourierfFlux_k
    do j=1, NbIncTotalPrim_ctx(IncFrac(k)%ic)
       divFourierFlux_k(j) = sum_aks * divTemperatureFrac(j,k)
    end do

    SmFourierFlux = sum_aks * SmTemperatureFrac(k)

    ! divFourierFlux_r, r represent s' in paper, r is node
    do r=1, NbNodeFrac
       numr = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+r)

       do j=1, NbIncTotalPrim_ctx(IncNode(numr)%ic)
          divFourierFlux_r(j,r) = - TkFracLocal_Fourier(k)%pt(s,r) * divTemperatureNode(j,numr)
       end do

       SmFourierFlux = SmFourierFlux &
            - TkFracLocal_Fourier(k)%pt(s,r) * SmTemperatureNode(numr)
    end do

  end subroutine Jacobian_divFourierFlux_fracnode


  ! Regularization of Jacobian
  ! operation on JacBigA, bigSm
  subroutine Jacobian_Regularization

    integer :: k, rowk, colk

    ! rows of node own
    do k=1, NbNodeOwn_Ncpus(commRank+1)

       rowk = k ! row of k
       colk = rowk ! col of diag element in rowk

       ! regularization for node
       call Jacobian_Regularization_row(k, rowk, colk, "n")
    end do

    ! rows of frac own
    do k=1, NbFracOwn_Ncpus(commRank+1)

       rowk = k + NbNodeOwn_Ncpus(commRank+1) ! row of k
       colk = k + NbNodeLocal_Ncpus(commRank+1) ! col of diag element in rowk

       ! regularization for frac
       call Jacobian_Regularization_row(k, rowk, colk, "f")
    end do

    ! rows of cell
    do k=1, NbCellLocal_Ncpus(commRank+1)

       rowk = k + NbNodeOwn_Ncpus(commRank+1) &
            + NbFracOwn_Ncpus(commRank+1) ! row of k

       colk = k + NbNodeLocal_Ncpus(commRank+1) &
            + NbFracLocal_Ncpus(commRank+1)! col of diag element in rowk

       ! regularization for cell
       call Jacobian_Regularization_row(k, rowk, colk, "c")
    end do


  end subroutine Jacobian_Regularization


  ! sub subroutine of Jacobian_Regularization
  ! used for regularization for node/frac/cell
  subroutine Jacobian_Regularization_row(k, rowk, colk, cv)

    ! row and col of diag block element in JacBigA
    integer, intent(in) :: k, rowk, colk
    character, intent(in) :: cv

    integer :: i, icp, iph, j, nz, nzj
    integer :: errcode, Ierr
    double precision :: sumcol

    ! look for the diag: JacBigA(:,:,nz)
    do i=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
       if( JacBigA%Num(i)==colk) then
          nz = i
          exit
       end if
    end do

    do i=1, NbCompThermique

       ! look for col nul in JacBigA(:,:,nz)
       ! ps. index order of matrix JacBigA(:,:,nz) is (col,row)
       sumcol = 0.d0 ! sum of col i
       do j=1, NbCompThermique
          sumcol = sumcol + abs( JacBigA%Val(i,j,nz))
       end do

       ! warning
       ! TODO: CHECKME: the following is a magic number
       ! eps*1e-3 for small permeability in the pressure equation
       ! TODO: find an automatic scaling
       if(sumcol<eps*1e-3) then ! col i is null

          ! look for component C_{i}^alpha corresponding to the col i

          if(cv .eq.'n') then ! node
             j = NumIncTotalPrimNode(i,k)
             icp = 0
             iph = 0

             if (j<=NbIncPTC_ctx(IncNode(k)%ic)) then
                icp = NumIncPTC2NumIncComp_comp_ctx(j,IncNode(k)%ic)
                iph = NumIncPTC2NumIncComp_phase_ctx(j,IncNode(k)%ic)
             endif
          else if(cv .eq. 'f') then ! fracs
             j = NumIncTotalPrimFrac(i,k)
             icp = NumIncPTC2NumIncComp_comp_ctx(j,IncFrac(k)%ic)
             iph = NumIncPTC2NumIncComp_phase_ctx(j,IncFrac(k)%ic)
          else if(cv .eq. 'c') then ! cell
             j = NumIncTotalPrimCell(i,k)
             icp = NumIncPTC2NumIncComp_comp_ctx(j,IncCell(k)%ic)
             iph = NumIncPTC2NumIncComp_phase_ctx(j,IncCell(k)%ic)
          else
             print*, ""
             print*, "Regularization error: cv should be node/frac/cell"
             errcode = 51
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if(icp .eq. 0) then
             print*, " "
             ! print*, JacBigA%Val(:,1,nz)
             ! print*, JacBigA%Val(:,2,nz)
             print*, "Regularization error: icp=0"
             write(*,'(A,4I8)') "    row/col/i/commRank=", rowk, colk, i, commRank
             errcode = 52
             write(*,*) "j ", j, "cv ", cv, " ic ", IncCell(k)%ic
             write(*,*) "sumcol ", sumcol
             write(*,*) "i", i
             ! write(*,*) "x ", XFaceLocal(:,FracToFaceLocal(k))

             ! write(*,*) "Pression ", IncFrac(k)%Pression
             ! write(*,*) "Cwater: ", IncFrac(k)%Comp(:,1)
             ! write(*,*) "Coil  : ", IncFrac(k)%Comp(:,2)
             ! write(*,*) "Satu  : ", IncFrac(k)%Saturation(:)

             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          endif

          ! row icp in JacBigA is set to be zero
          do j=JacBigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)
             nzj = JacBigA%Num(j)
             JacBigA%Val(:,icp,nzj) = 0.d0
          end do

          ! JacBigA(i,icp,nz)=1
          JacBigA%Val(i,icp,nz) = 1.d0

          ! bigSm(i,rowk) = 0
          bigSm(i,rowk) = 0.d0
       end if
    end do

  end subroutine Jacobian_Regularization_row


  ! Alignment of Jacobian: diag method
  ! operation on JacA, Sm
  subroutine Jacobian_Alignment_diag

    integer :: k, rowk, colk

    ! rows of node own
    do k=1, NbNodeOwn_Ncpus(commRank+1)

       rowk = k ! row of k
       colk = rowk ! col of diag element in rowk

       call Jacobian_Alignment_diag_row(rowk, colk)
    end do

    ! rows of frac own
    do k=1, NbFracOwn_Ncpus(commRank+1)

       rowk = k + NbNodeOwn_Ncpus(commRank+1)   ! row of k
       colk = k + NbNodeLocal_Ncpus(commRank+1) ! col of diag element in rowk

       call Jacobian_Alignment_diag_row(rowk, colk)
    end do

  end subroutine Jacobian_Alignment_diag


  ! sub subroutine of Jacobian_Alignment: diag method
  ! used for Alignment for node/frac
  subroutine Jacobian_Alignment_diag_row(rowk, colk)

    integer, intent(in) :: rowk, colk
    integer :: i, nz, errcode, Ierr

    ! optimal size of the WORK array, pre-queried for machine
    integer, parameter :: lwork &
         = 320 ! Mac

    double precision, dimension(lwork) :: work
    integer :: ipival(NbCompThermique), info

    double precision, dimension(NbCompThermique, NbCompThermique) :: &
         AA, BB
    double precision, dimension(NbCompThermique) :: &
         Smk

    ! look for the diag: JacA(:,:,nz)
    do i=JacA%Pt(rowk)+1, JacA%Pt(rowk+1)
       if( JacA%Num(i)==colk) then
          nz = i
          exit
       end if
    end do

    ! BB = inv(JacA%Val(:,:,nz))
    ! ps. the index order of JacA%Val(:,:,nz) is (col, row)
    ! so the index order of BB is also (col, row)

    BB = JacA%Val(:,:,nz)
    call dgetrf(NbCompThermique, NbCompThermique, &
         BB, NbCompThermique, ipival, info)
    if(info/=0) then
       print*, "dgetrf error", info, "in Alignment, rowk/colk = ", rowk, colk
       print*, "shape of BB", shape(BB)
       print*, BB
       print *, "Try to locate error..."
       call Jacobian_JacBigA_locate_frac_row(rowk)
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    call dgetri(NbCompThermique, BB, NbCompThermique, &
         ipival, work, lwork, info)
    ! print*, lwork(1)
    if(info /= 0) then
       print*, "dgetri error", info, "in Alignment, rowk/colk = ", rowk, colk
       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    ! JacA%Val(:,:,i) = JacA%Val(:,:,i)*BB, row k
    do i=JacA%Pt(rowk)+1, JacA%Pt(rowk+1)

       AA(:,:) = JacA%Val(:,:,i)
       call dgemm('N', 'N', NbCompThermique, NbCompThermique, NbCompThermique, &
            1.d0, AA, NbCompThermique, BB, NbCompThermique, 0.d0, &
            JacA%Val(:,:,i), NbCompThermique)
    end do

    ! Sm(:,rowk) = BB * Sm(:,rowk), rowk
    ! transpose of BB is necessary since the index of BB is (col, row)
    Smk(:) = Sm(:,rowk)
    call dgemv('T', NbCompThermique, NbCompThermique, 1.d0, &
         BB, NbCompThermique, Smk, 1, &
         0.d0, Sm(:,rowk), 1)

  end subroutine Jacobian_Alignment_diag_row


  ! Alignment of Jacobian: manually method
  ! operation on JacA, Sm
  subroutine Jacobian_Alignment_man

    integer :: k, rowk

    ! rows of node own
    do k=1, NbNodeOwn_Ncpus(commRank+1)


       ! WARNING: l'alignement avec alignemat ne s'applique qu'aux eqs de conservation
       ! donc pas aux noeuds DIR DIR
       ! TODO: le cas DIR Neu ou Neu Dir reste a faire

       if ( (IdNodeLocal(k)%P.ne."d").and.(IdNodeLocal(k)%T.ne."d") ) then
       rowk = k ! row of k
       call Jacobian_Alignment_man_row(k, rowk, IncNode(k)%ic)

       else if ( (IdNodeLocal(k)%P.eq."d").and.(IdNodeLocal(k)%T.ne."d") ) then
          call CommonMPI_abort('in manual alignment Jacobian &
            mix Dir/Neu node are not implemented')
       else if ( (IdNodeLocal(k)%T.eq."d").and.(IdNodeLocal(k)%P.ne."d") ) then
          call CommonMPI_abort('in manual alignment Jacobian &
            mix Dir/Neu node are not implemented')
       endif
    end do

    ! rows of frac own
    do k=1, NbFracOwn_Ncpus(commRank+1)

       rowk = k + NbNodeOwn_Ncpus(commRank+1) ! row of k
       call Jacobian_Alignment_man_row(k, rowk, IncFrac(k)%ic)
    end do

  end subroutine Jacobian_Alignment_man


  ! sub subroutine of Jacobian_Alignment: manually method
  ! used for Alignment for node/frac
  subroutine Jacobian_Alignment_man_row(k, rowk, ic)

    integer, intent(in) :: k, rowk, ic
    integer :: i
!!$    integer :: j, nz 

    double precision, dimension(NbCompThermique, NbCompThermique) :: &
         AA, BB
    double precision, dimension(NbCompThermique) :: &
         Smk

    ! the index order of JacA%Val(:,:,nz) is (col, row)
    ! the index order of aligmethod(:,:,ic) is also (col, row)


!!$    do i=JacA%Pt(rowk)+1, JacA%Pt(rowk+1)
!!$       if (JacA%Num(i)==rowk) then
!!$          nz = i
!!$       endif
!!$    enddo
!!$   
!!$    write(*,*)
!!$    BB(:,:) = aligmat(:,:,ic)        
!!$    write(*,*)
!!$    write(*,*)' ic',ic
!!$    do i=1,NbCompThermique       
!!$    write(*,*)' ligne i ',i
!!$       write(*,*)' Aligne ',(BB(i,j),j=1,NbCompThermique)
!!$    enddo      
!!$    
!!$    AA(:,:) = JacA%Val(:,:,nz)
!!$    write(*,*)
!!$    do i=1,NbCompThermique       
!!$    write(*,*)' ligne i ',i
!!$       write(*,*)' AA avant ',(AA(i,j),j=1,NbCompThermique)
!!$    enddo


    
    BB(:,:) = aligmat(:,:,ic)

    ! JacA%Val(:,:,i) = JacA%Val(:,:,i) * aligmethod(:,:,ic)
    ! since all the matrix are transpose
    do i=JacA%Pt(rowk)+1, JacA%Pt(rowk+1)

       AA(:,:) = JacA%Val(:,:,i)
       call dgemm('N', 'N', NbCompThermique, NbCompThermique, NbCompThermique, &
            1.d0, AA, NbCompThermique, BB, NbCompThermique, 0.d0, &
            JacA%Val(:,:,i), NbCompThermique)
    end do

    ! Sm(:,rowk) = BB * Sm(:,rowk), rowk
    ! transpose of aligmethod(:,:,ic) is necessary
    ! since the index order of BB is (col, row)
    Smk(:) = Sm(:,rowk)
    call dgemv('T', NbCompThermique, NbCompThermique, 1.d0, &
         BB, NbCompThermique, Smk, 1, &
         0.d0, Sm(:,rowk), 1)


!!$    AA(:,:) = JacA%Val(:,:,nz)
!!$    write(*,*)
!!$    do i=1,NbCompThermique       
!!$    write(*,*)' ligne i ',i
!!$       write(*,*)' AA apres ',(AA(i,j),j=1,NbCompThermique)
!!$    enddo


    


    
  end subroutine Jacobian_Alignment_man_row


  ! Schur complement
  ! (A1, A2; A3, A4) ->  A1-A3*A2/A4
  ! BigSm =(b1;b2)   ->  b1-A3*(b2/14), b1, b2 are vectors

  ! node/frac,  cell,  well
  !  |   A1,     A2,   Anw  | node own/frac
  !  |   A3,     A4,   0    | cell
  !  |   Awn,    0,    Aww  | well own
  ! size of 11 is (NbUnkOwn, NbUnkLocal)
  ! size of 14 is (NbCellLocal, NbCellLocal)
  ! A4 is a diagnal matrix

  ! Some index of code
  !  row i in A1
  !  row lmi in A3
  !  col lmj in A3
  !  col n in A1
  !  the element considered in A4 is (lmj,lmi), %Num(mij)

  subroutine Jacobian_Schur

    integer :: k, rowk, mk
    integer :: i, mi, rowmi, colmi, kmi, mj, lmj, mij, n, ln, j, lj, ni, info
    integer :: row1, row2, rowi
    integer :: errcode, Ierr
    integer :: NbNodeFracOwn, NbNodeFracOwnCellLocal, NbNodeFracLocal, NbNodeFracCellLocal

    double precision :: DB(NbCompThermique, NbCompThermique)

    NbNodeFracLocal = NbNodeLocal_Ncpus(commRank+1) + NbFracLocal_Ncpus(commRank+1)
    NbNodeFracCellLocal = NbNodeFracLocal + NbCellLocal_Ncpus(commRank+1)

    NbNodeFracOwn = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)
    NbNodeFracOwnCellLocal = NbNodeFracOwn + NbCellLocal_Ncpus(commRank+1)

    if(.not. allocated(ipiv)) then
       allocate(ipiv(NbCompThermique, NbCellLocal_Ncpus(commRank+1)) )
    end if

    ! compute LU for each block element of A4
    do k=1, NbCellLocal_Ncpus(commRank+1)
       rowk = k + NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)

       mk = JacBigA%Pt(rowk+1) ! last element in row is cell

       ! diag block is JacBigA%Val(:,:,mk)
       call dgetrf(NbCompThermique, NbCompThermique, &
            JacBigA%Val(:,:,mk), NbCompThermique, ipiv(:,k), info)

       if(info /= 0) then
          write(0,'(A,I5)', advance='no') "dgetrf error in Schur complement, k =", k
          write(0,'(A,I5)') ",  info =", info
          call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
       end if
    end do

    ! loop of rows for node own not dir and frac own
    do i=1, NbNodeFracOwn  ! rows in the part A1

       ! cycle if the row is dir node
#ifdef _THERMIQUE_
    if(i<=NbNodeOwn_Ncpus(commRank+1)) then
        if(IdNodeLocal(i)%P=="d" .and. (IdNodeLocal(i)%T=="d")) then
            cycle
        end if
    end if
#else
    if(i<=NbNodeOwn_Ncpus(commRank+1)) then
        if(IdNodeLocal(i)%P=="d") then
            cycle
        end if
    end if
#endif

       do mi=JacBigA%Pt(i)+1, JacBigA%Pt(i+1) ! loop of non zeros in row i
          colmi = JacBigA%Num(mi) ! %Num(mi)

          if ((colmi>NbNodeFracLocal) .and. (colmi<=NbNodeFracCellLocal)) then ! cols in the part (A2;A4), part cell

             kmi = colmi - NbNodeFracLocal ! cell kmi is in column colmi
             rowmi = kmi + NbNodeFracOwn   ! cell kmi is in row    rowmi

             mij = JacBigA%Pt(rowmi+1) ! last nonzero element of row colmi since A4 is diagnal

             ! DB = JacBigA%Val(mij)^-1*JacBigA%Val(mi)
             DB(:,:) = JacBigA%Val(:,:,mi)
             call dgetrs('N', NbCompThermique, NbCompThermique, &
                  JacBigA%Val(:,:,mij), NbCompThermique, ipiv(:,kmi), &
                  DB, NbCompThermique, info)

             if(info /= 0) then
                write(0,'(A,I5)', advance='no') "dgetrs error in Schur complement to compute DB, mi =", mi
                write(0,'(A,I5)') "  info =", info
                call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
             end if

             ! second member
             ! BigSm(:,i) = BigSm(:,i) - BigSm(:,rowmi)*JacBigA%Val(mi)/JacBigA%Val(mij)
             call dgemv('T', NbCompThermique, NbCompThermique, -1.d0, &
                  DB, NbCompThermique, BigSm(:,rowmi), 1, &
                  1.d0, BigSm(:,i), 1)

             ! matrix
             do mj=JacBigA%Pt(rowmi)+1, JacBigA%Pt(rowmi+1) ! row rowmi is in (A3,A4)
                lmj = JacBigA%Num(mj) ! %Num(mj)

                if(lmj<=NbNodeFracLocal) then ! cols in the part (A1;A3)

                   ! schur compement for A1(i,:)
                   do n=JacBigA%Pt(i)+1, JacBigA%Pt(i+1)
                      ln = JacBigA%Num(n)

                      if(ln==lmj) then

                         ! JacBigA%Val(n) = JacBigA%Val(n) - JacBigA%Val(mj)*JacBigA%Val(mi)/JacBigA%Val(mij)
                         call dgemm('N', 'N', NbCompThermique, NbCompThermique, NbCompThermique, &
                              -1.d0, JacBigA%Val(:,:,mj), NbCompThermique, DB, NbCompThermique, 1.d0, &
                              JacBigA%Val(:,:,n), NbCompThermique)
                      end if
                   end do

                end if
             end do

          end if
       end do
    end do

    ! JacA is (A1, Anw; Awn, Aww)
    do i=1, NbNodeOwn_Ncpus(commRank+1) ! node own

       ni = CellbyNodeOwn%Pt(i+1) - CellbyNodeOwn%Pt(i)

       do j=1, JacA%Pt(i+1)-JacA%Pt(i)
          lj = JacA%Num( JacA%Pt(i)+j) ! col

          ! A1, row node, col node/frac
          if(lj<=NbNodeFracLocal) then
             JacA%Val(:,:,JacA%Pt(i)+j) = JacBigA%Val(:,:,JacBigA%Pt(i)+j)

          ! Anw, row node, col well
          else
             JacA%val(:,:,JacA%Pt(i)+j) = JacBigA%Val(:,:,JacBigA%Pt(i)+ni+j) ! +ni to jump column cell
          end if
       end do
    end do


    ! A1, row frac, col node/frac
    do i=1, NbFracOwn_Ncpus(commRank+1) ! frac own
       rowi = i + NbNodeOwn_Ncpus(commRank+1)

       do j=1, JacA%Pt(rowi+1)-JacA%Pt(rowi)
          lj = JacA%Num( JacA%Pt(rowi)+j) ! col

          JacA%Val(:,:,JacA%Pt(rowi)+j) = JacBigA%Val(:,:,JacBigA%Pt(rowi)+j)
       end do
    end do

    ! Awn and Aww
    do i=1, NbWellInjOwn_Ncpus(commRank+1) + NbWellProdOwn_Ncpus(commRank+1)

       row1 = i + NbNodeFracOwn ! row in JacA
       row2 = i + NbNodeFracOwnCellLocal ! row in JacBigA

       do j=1, JacA%Pt(row1+1)-JacA%Pt(row1) ! = JacBigA%Pt(row2+1)-JacBigA%Pt(row2)
          JacA%Val(:,:,JacA%Pt(row1)+j) = JacBigA%Val(:,:,JacBigA%Pt(row2)+j)
       end do
    end do

    ! Sm(:,:)
    do i=1, NbNodeFracOwn ! node and frac
       Sm(:,i) = BigSm(:,i)
    end do

    do i=1, NbWellInjOwn_Ncpus(commRank+1)+NbWellProdOwn_Ncpus(commRank+1) ! well

       row1 = i + NbNodeFracOwn ! row in Sm
       row2 = i + NbNodeFracOwnCellLocal ! row in BigSm

       Sm(:,row1) = BigSm(:,row2)
    end do

  end subroutine Jacobian_Schur


  ! Non-zero sturcture of Jacobian
  subroutine Jacobian_StrucJacBigA

    !            node, frac, cell, wellinj, wellprod
    !           | a11, a12, a13, a14, a15 | node own
    ! JacBigA = | a21, a22, a23, 0,   0   | frac own
    !           | a31, a32, a33, 0,   0   | cell (own and ghost)
    !           ! a41, 0,   0,   a44, 0   | wellinj own
    !           ! a51, 0,   0,   0,   a55 | wellprod own

    !            node, frac, wellinj, wellprod
    ! JacA =    | A11, A12, A13, A14 | node own
    !           | A21, A22, 0,   0   | frac own
    !           | A31, 0,   A33, 0   | wellinj own,
    !           | A41, 0,   0,   A44 | wellprod own
    ! the nonzero structure of Aij is same as aij, i,j=1,2 in JacBigA
    ! A31=a31, A33=a33, A41=a41, A44=a44

    ! The non zero structure of aij is based on the connectivities
    ! a11: NodebyNodeOwn if not dirichlet
    ! a12: FracbyNodeOwn if not dir
    ! a13: CellbyNodeOwn if not dir
    ! a14: WellInjbyNodeOwn if not dir
    ! a15: WellProdbyNodeOwn if not dir
    ! Rq: if node own is Dir then only one non zero/row for a1*

    ! a21: NodebyFracOwn
    ! a22: FracbyFracOwn
    ! a23: CellbyFracOwn

    ! a31: NodebyCellLocal
    ! a32: FracbyCellLocal
    ! a33: diag block

    ! a41: NodebyWellInjLocal
    ! a44: diag
    ! a51: NodebyWellProdLocal
    ! a55: diag

    ! Four steps in this subroutine
    !   Number of non zeros each line: NbNnzbyline
    !   bigA%Nb, bigA%Pt using NbNnzbyline
    !   bigA%Num, non zero structure of bigA
    !   arrange bigA%Num such that in each row, the cols is in inscreasing order

    integer, dimension(:), allocatable :: &
         nbNnzbyline ! number of non zeros each line

    integer :: i, j, k, jf, Nz, start
    integer :: tmp

    integer :: &
         nbNodeOwn, nbFracOwn, nbWellInjOwn, nbWellProdOwn, &
         nbNodeLocal, nbFracLocal, nbCellLocal, nbWellInjLocal, nbWellProdLocal

    nbNodeOwn = NbNodeOwn_Ncpus(commRank+1)
    nbFracOwn = NbFracOwn_Ncpus(commRank+1)
    nbWellInjOwn = NbWellInjOwn_Ncpus(commRank+1)
    nbWellProdOwn = NbWellProdOwn_Ncpus(commRank+1)

    nbNodeLocal = NbNodeLocal_Ncpus(commRank+1)
    nbFracLocal = NbFracLocal_Ncpus(commRank+1)
    nbCellLocal = NbCellLocal_Ncpus(commRank+1)
    nbWellInjLocal = NbWellInjLocal_Ncpus(commRank+1)
    nbWellProdLocal = NbWellProdLocal_Ncpus(commRank+1)

    JacBigA%Nb = nbNodeOwn + nbFracOwn + nbCellLocal + nbWellInjOwn + nbWellProdOwn
    allocate(nbNnzbyLine(JacBigA%Nb))
    allocate(JacBigA%Pt(JacBigA%Nb+1))
    allocate(BigSm(NbCompThermique,JacBigA%Nb))

    ! a1* in JacbigA
    do i=1, nbNodeOwn

       ! Darcy dir and T dir
       ! only in this case, one non-zero (block) this row, it is Id
#ifdef _THERMIQUE_
       if( (IdNodeLocal(i)%P=="d") .and. (IdNodeLocal(i)%T=="d")) then
#else
       if( (IdNodeLocal(i)%P=="d")) then
#endif
          nbNnzbyLine(i) = 1
       else
          nbNnzbyLine(i) = NodebyNodeOwn%Pt(i+1) - NodebyNodeOwn%Pt(i) &
               + FracbyNodeOwn%Pt(i+1) - FracbyNodeOwn%Pt(i) &
               + CellbyNodeOwn%Pt(i+1) - CellbyNodeOwn%Pt(i) &
               + WellInjbyNodeOwn%Pt(i+1) - WellInjbyNodeOwn%Pt(i) &
               + WellProdbyNodeOwn%Pt(i+1) - WellProdbyNodeOwn%Pt(i)
       end if
    end do

    ! if(commRank==1) then
    !    print*, nbNnzbyLine(1:NbNodeOwn_Ncpus(commRank+1))
    ! end if

    ! a2* in JacbigA, rq: frac is not dir face
    start = nbNodeOwn
    do i=1, nbFracOwn
       nbNnzbyLine(start+i) = NodebyFracOwn%Pt(i+1) - NodebyFracOwn%Pt(i) &
            + FracbyFracOwn%Pt(i+1) - FracbyFracOwn%Pt(i) &
            + CellbyFracOwn%Pt(i+1) - CellbyFracOwn%Pt(i)
    end do

    ! a3* in JacbigA
    start = start + nbFracOwn
    do i=1, nbCellLocal
       nbNnzbyLine(start+i) = &
            NodebyCellLocal%Pt(i+1) - NodebyCellLocal%Pt(i)  &
            + FracbyCellLocal%Pt(i+1) - FracbyCellLocal%Pt(i) + 1
    end do

    ! a4* in JacbigA
    start = start + nbCellLocal
    do i=1, nbWellInjOwn
        nbNnzbyLine(start+i) = &
            NodebyWellInjLocal%Pt(i+1) - NodebyWellInjLocal%Pt(i) + 1 ! a14 and a44, (+1 sicne a44 is diag)
    end do

    ! a5* in JacbigA
    start = start + nbWellInjOwn
    ! print*, 'DEBUG HERE', NbWellProdOwn_Ncpus(commRank+1)
    do i=1, nbWellProdOwn
        nbNnzbyLine(start+i) = &
            NodebyWellProdLocal%Pt(i+1) - NodebyWellProdLocal%Pt(i) + 1 ! a55 is diag
    end do

    ! JacBigA%Pt
    Nz = 0
    JacBigA%Pt(1) = 0
    do i=1, JacBigA%Nb
       Nz = Nz + nbNnzbyLine(i)
       JacBigA%Pt(i+1) = Nz
    end do

    deallocate(nbNnzbyLine)

    ! JacBigA%Num
    ! Rq: for Fracby*, *=CellLocal/FracOwn/NodeOwn
    !   jf=Fracby*(jf) is face number, not frac number,
    !   FaceToFrac(jf) is frac number

    allocate(JacBigA%Num(Nz))
    allocate(JacBigA%Val(NbCompThermique, NbCompThermique, Nz))

    start = 0

    do i=1, nbNodeOwn

       ! dir Darcy and dir T
#ifdef _THERMIQUE_
       if((IdNodeLocal(i)%P=="d") .and. (IdNodeLocal(i)%T=="d")) then
#else
       if((IdNodeLocal(i)%P=="d")) then
#endif

          JacBigA%Num(start+1) = i ! node=(node own, node ghost)
          start = start + 1

       else ! one of Darcy or T is not dir
          ! a11(i,:)
          do j=1, NodebyNodeOwn%Pt(i+1)-NodebyNodeOwn%Pt(i)
             JacBigA%Num(start+j) = NodebyNodeOwn%Num(j+NodebyNodeOwn%Pt(i))
          end do
          start = start + NodebyNodeOwn%Pt(i+1)-NodebyNodeOwn%Pt(i)

          ! a12(i,:)
          do j=1, FracbyNodeOwn%Pt(i+1)-FracbyNodeOwn%Pt(i)
             jf = FracbyNodeOwn%Num(j+FracbyNodeOwn%Pt(i)) ! jf is face num, need to transform to frac num
             JacBigA%Num(start+j) = FaceToFracLocal(jf) + nbNodeLocal ! col
          end do
          start = start + FracbyNodeOwn%Pt(i+1)-FracbyNodeOwn%Pt(i)

          ! a13(i,:)
          do j=1, CellbyNodeOwn%Pt(i+1)-CellbyNodeOwn%Pt(i)
             JacBigA%Num(start+j) = CellbyNodeOwn%Num(j+CellbyNodeOwn%Pt(i)) + nbNodeLocal + nbFracLocal ! col
          end do
          start = start + CellbyNodeOwn%Pt(i+1)-CellbyNodeOwn%Pt(i)

          ! a14(i,:)
          do j=1, WellInjbyNodeOwn%Pt(i+1)-WellInjbyNodeOwn%Pt(i)
             JacBigA%Num(start+j) = WellInjbyNodeOwn%Num(j+WellInjbyNodeOwn%Pt(i)) &
                  + nbNodeLocal + nbFracLocal + nbCellLocal
          end do
          start = start + WellInjbyNodeOwn%Pt(i+1)-WellInjbyNodeOwn%Pt(i)

          ! a15(i;:)
          do j=1, WellProdbyNodeOwn%Pt(i+1)-WellProdbyNodeOwn%Pt(i)
             JacBigA%Num(start+j) = WellProdbyNodeOwn%Num(j+WellProdbyNodeOwn%Pt(i)) &
                  + nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal
          end do
          start = start + WellProdbyNodeOwn%Pt(i+1)-WellProdbyNodeOwn%Pt(i)
       end if
    end do

    do i=1, nbFracOwn

       ! a21(i,:)
       do j=1, NodebyFracOwn%Pt(i+1)-NodebyFracOwn%Pt(i)
          JacBigA%Num(start+j) = NodebyFracOwn%Num(j+NodebyFracOwn%Pt(i))
       end do
       start = start + NodebyFracOwn%Pt(i+1)-NodebyFracOwn%Pt(i)

       ! a22(i,:)
       do j=1, FracbyFracOwn%Pt(i+1)-FracbyFracOwn%Pt(i)
          jf = FracbyFracOwn%Num(j+FracbyFracOwn%Pt(i))
          JacBigA%Num(start+j) =  FaceToFracLocal(jf) + nbNodeLocal ! col
       end do
       start = start + FracbyFracOwn%Pt(i+1)-FracbyFracOwn%Pt(i)

       ! a23(i,:)
       do j=1, CellbyFracOwn%Pt(i+1)-CellbyFracOwn%Pt(i)
          JacBigA%Num(start+j) = CellbyFracOwn%Num(j+CellbyFracOwn%Pt(i)) + nbNodeLocal + nbFracLocal ! col
       end do
       start = start + CellbyFracOwn%Pt(i+1)-CellbyFracOwn%Pt(i)

       ! a24 = a25 = 0
    end do

    do i=1, nbCellLocal

       ! a31(i,:)
       do j=1, NodebyCellLocal%Pt(i+1)-NodebyCellLocal%Pt(i)
          JacBigA%Num(start+j) = NodebyCellLocal%Num(j+NodebyCellLocal%Pt(i))
       end do
       start = start + NodebyCellLocal%Pt(i+1)-NodebyCellLocal%Pt(i)

       ! a32(i,:)
       do j=1, FracbyCellLocal%Pt(i+1)-FracbyCellLocal%Pt(i)
          jf = FracbyCellLocal%Num(j+FracbyCellLocal%Pt(i))
          JacBigA%Num(start+j) = FaceToFracLocal(jf) + nbNodeLocal ! col
       end do
       start = start + FracbyCellLocal%Pt(i+1)-FracbyCellLocal%Pt(i)

       ! a33(i,:)
       JacBigA%Num(start+1) = i + nbNodeLocal + nbFracLocal ! col
       start = start + 1

       ! a34 = a35 = 0
    end do

    do i=1, nbWellInjOwn

       ! a41(i,:)
       do j=1, NodebyWellInjLocal%Pt(i+1) - NodebyWellInjLocal%Pt(i)
          JacBigA%Num(start+j) = NodebyWellInjLocal%Num(j+NodebyWellInjLocal%Pt(i))
       end do
       start = start + NodebyWellInjLocal%Pt(i+1) - NodebyWellInjLocal%Pt(i)

       ! a44(i,:)
       JacBigA%Num(start+1) = i + nbNodeLocal + nbFracLocal + nbCellLocal
       start = start + 1
    end do

    do i=1, nbWellProdOwn

       ! a51(i,:)
       do j=1, NodebyWellProdLocal%Pt(i+1) - NodebyWellProdLocal%Pt(i)
          JacBigA%Num(start+j) = NodebyWellProdLocal%Num(j+NodebyWellProdLocal%Pt(i))
       end do
       start = start + NodebyWellProdLocal%Pt(i+1) - NodebyWellProdLocal%Pt(i)

       ! a55(i,:)
       JacBigA%Num(start+1) = i + nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal
       start = start + 1
    end do

    ! Sort JacbigA%Num
    do i=1,JacBigA%Nb

       do k=1, JacBigA%Pt(i+1)-JacBigA%Pt(i)
          do j=JacBigA%Pt(i)+1, JacBigA%Pt(i+1)-k

             if( JacBigA%Num(j) > JacBigA%Num(j+1) ) then
                tmp = JacBigA%Num(j)
                JacBigA%Num(j) = JacBigA%Num(j+1)
                JacBigA%Num(j+1) = tmp
             end if

          end do
       end do
    end do

  end subroutine Jacobian_StrucJacBigA


  ! Non zero sturcture of JacA: JacBigA after Schur
  subroutine Jacobian_StrucJacA

    !         node, frac, wellinj, wellprod
    ! JacA = | A11, A12, A13, A14 | node own
    !        | A21, A22, 0,   0   | frac own
    !        | A31, 0,   A33, 0   | wellinj own
    !        | A41, 0,   0,   A44 | wellprod own
    ! the nonzero structure of Aij is same as aij in JacBigA, i,j=1,2

    ! The non zero structure of aij is based on the connectivities
    ! A11: NodebyNodeOwn if not dir
    ! A12: FracbyNodeOwn if not dir
    ! Rq: if node own is Dir then only one non zero per row for A1*

    ! A21: NodebyFracOwn
    ! A22: FracbyFracOwn
    ! A13: WellInjbyNodeOwn
    ! A14: WellProdbyNodeOwn

    ! A31: NodebyWellInjLocal
    ! A41: NodebyWellProdLocal
    ! A33, A44: diag

    ! Four steps in this subroutine
    !   Number of non zeros each line: NbNnzbyline
    !   bigA%Nb, bigA%Pt using NbNnzbyline
    !   bigA%Num, non zero structure of bigA
    !   arrange bigA%Num such that in each row, the cols is in inscreasing order

    integer, dimension(:), allocatable :: &
         nbNnzbyline ! number of non zeros each line

    integer :: i, j, k, jf, Nz, start
    integer :: tmp

    integer :: &
         nbNodeOwn, nbFracOwn, nbWellInjOwn, nbWellProdOwn, &
         nbNodeLocal, nbFracLocal, nbCellLocal, nbWellInjLocal, nbWellProdLocal

    nbNodeOwn = NbNodeOwn_Ncpus(commRank+1)
    nbFracOwn = NbFracOwn_Ncpus(commRank+1)
    nbWellInjOwn = NbWellInjOwn_Ncpus(commRank+1)
    nbWellProdOwn = NbWellProdOwn_Ncpus(commRank+1)

    nbNodeLocal = NbNodeLocal_Ncpus(commRank+1)
    nbFracLocal = NbFracLocal_Ncpus(commRank+1)
    nbCellLocal = NbCellLocal_Ncpus(commRank+1)
    nbWellInjLocal = NbWellInjLocal_Ncpus(commRank+1)
    nbWellProdLocal = NbWellProdLocal_Ncpus(commRank+1)

    ! JacA%Nb
    JacA%Nb = nbNodeOwn + nbFracOwn + nbWellInjOwn + nbWellProdOwn

    allocate(nbNnzbyLine(JacA%Nb))
    allocate(JacA%Pt(JacA%Nb+1))
    allocate(Sm(NbCompThermique, JacA%Nb))
    allocate(csrK(nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal + nbWellProdLocal))
    csrK(:) = 0
    allocate(csrSR(nbNodeLocal + nbFracLocal + nbCellLocal + nbWellInjLocal + nbWellProdLocal))
    csrSR(:) = 0

    ! A1* in JacA
    do i=1, nbNodeOwn

       ! Darcy dir and T dir
#ifdef _THERMIQUE_
       if( (IdNodeLocal(i)%P=="d") .and. (IdNodeLocal(i)%T=="d")) then !
#else
       if( (IdNodeLocal(i)%P=="d")) then !
#endif
          nbNnzbyLine(i) = 1 ! one non zero (block) this row
       else
          nbNnzbyLine(i) = NodebyNodeOwn%Pt(i+1) - NodebyNodeOwn%Pt(i) &
               + FracbyNodeOwn%Pt(i+1) - FracbyNodeOwn%Pt(i) &
               + WellInjbyNodeOwn%Pt(i+1) - WellInjbyNodeOwn%Pt(i) &
               + WellProdbyNodeOwn%Pt(i+1) - WellProdbyNodeOwn%Pt(i)
       end if
    end do

    ! A2* in JacA, rq: frac is not dir face
    start = nbNodeOwn
    do i=1,nbFracOwn
       nbNnzbyLine(start+i) = NodebyFracOwn%Pt(i+1) - NodebyFracOwn%Pt(i) &
            + FracbyFracOwn%Pt(i+1) - FracbyFracOwn%Pt(i)
    end do

    ! A3* in JacA
    start = start + nbFracOwn
    do i=1, nbWellInjOwn
       nbNnzbyLine(start+i) = &
            NodebyWellInjLocal%Pt(i+1) - NodebyWellInjLocal%Pt(i) + 1 ! A44 is di
    end do

    ! A4* in JacA
    start = start + nbWellInjOwn
    do i=1, nbWellProdOwn
       nbNnzbyLine(start+i) = &
            NodebyWellProdLocal%Pt(i+1) - NodebyWellProdLocal%Pt(i) + 1 ! A55 is diag
    end do


    ! JacA%Pt
    Nz = 0
    JacA%Pt(1) = 0
    do i=1, JacA%Nb
       Nz = Nz + nbNnzbyLine(i)
       JacA%Pt(i+1) = Nz
    end do

    deallocate(nbNnzbyLine)

    ! JacA%Num
    ! Rq: for Fracby*, *=CellLocal/FracOwn/NodeOwn
    !   jf=Fracby*(jf) is face number, not frac number,
    !   FaceToFrac(jf) is frac number

    allocate( JacA%Num(Nz) )
    JacA%Num(:) = 0

    start = 0

    do i=1, nbNodeOwn

       ! Darcy and T are both dir
#ifdef _THERMIQUE_
       if((IdNodeLocal(i)%P=="d") .and. (IdNodeLocal(i)%T=="d")) then
#else
       if(IdNodeLocal(i)%P=="d") then
#endif
          JacA%Num(start+1) = i ! node=(node own, node ghost)
          start = start + 1

       else ! one of Darcy and T is not dir

          ! A11(i,:)
          do j=1, NodebyNodeOwn%Pt(i+1)-NodebyNodeOwn%Pt(i)
             JacA%Num(start+j) = NodebyNodeOwn%Num(j+NodebyNodeOwn%Pt(i))
          end do
          start = start + NodebyNodeOwn%Pt(i+1)-NodebyNodeOwn%Pt(i)

          ! A12(i,:)
          do j=1, FracbyNodeOwn%Pt(i+1)-FracbyNodeOwn%Pt(i)
             jf = FracbyNodeOwn%Num(j+FracbyNodeOwn%Pt(i)) ! jf is face num, need to transform to frac num
             JacA%Num(start+j) = FaceToFracLocal(jf) + nbNodeLocal ! col
          end do
          start = start + FracbyNodeOwn%Pt(i+1)-FracbyNodeOwn%Pt(i)

          ! A13(i,:)
          do j=1, WellInjbyNodeOwn%Pt(i+1)-WellInjbyNodeOwn%Pt(i)
             JacA%Num(start+j) = WellInjbyNodeOwn%Num(j+WellInjbyNodeOwn%Pt(i)) &
                  + nbNodeLocal + nbFracLocal
          end do
          start = start + WellInjbyNodeOwn%Pt(i+1)-WellInjbyNodeOwn%Pt(i)

          ! A14(i,:)
          do j=1, WellProdbyNodeOwn%Pt(i+1)-WellProdbyNodeOwn%Pt(i)
             JacA%Num(start+j) = WellProdbyNodeOwn%Num(j+WellProdbyNodeOwn%Pt(i)) &
                  + nbNodeLocal + nbFracLocal + nbWellInjLocal
          end do
          start = start + WellProdbyNodeOwn%Pt(i+1)-WellProdbyNodeOwn%Pt(i)
       end if
    end do

    do i=1, nbFracOwn

       ! A21(i,:)
       do j=1, NodebyFracOwn%Pt(i+1)-NodebyFracOwn%Pt(i)
          JacA%Num(start+j) = NodebyFracOwn%Num(j+NodebyFracOwn%Pt(i))
       end do
       start = start + NodebyFracOwn%Pt(i+1)-NodebyFracOwn%Pt(i)

       ! A22(i,:)
       do j=1, FracbyFracOwn%Pt(i+1)-FracbyFracOwn%Pt(i)
          jf = FracbyFracOwn%Num(j+FracbyFracOwn%Pt(i))
          JacA%Num(start+j) =  FaceToFracLocal(jf) + nbNodeLocal ! col
       end do
       start = start + FracbyFracOwn%Pt(i+1)-FracbyFracOwn%Pt(i)
    end do

    do i=1, nbWellInjOwn

       ! A31(i,:)
       do j=1, NodebyWellInjLocal%Pt(i+1)-NodebyWellInjLocal%Pt(i)
          JacA%Num(start+j) = NodebyWellInjLocal%Num(j+NodebyWellInjLocal%Pt(i))
       end do
       start = start + NodebyWellInjLocal%Pt(i+1)-NodebyWellInjLocal%Pt(i)

       ! A33(i,:)
       JacA%Num(start+1) = i + nbNodeLocal + nbFracLocal
       start = start + 1
    end do

    do i=1, nbWellProdOwn

       ! A41(i,:)
       do j=1, NodebyWellProdLocal%Pt(i+1)-NodebyWellProdLocal%Pt(i)
          JacA%Num(start+j) = NodebyWellProdLocal%Num(j+NodebyWellProdLocal%Pt(i))
       end do
       start = start + NodebyWellProdLocal%Pt(i+1)-NodebyWellProdLocal%Pt(i)

       ! A44(i,:)
       JacA%Num(start+1) = i + nbNodeLocal + nbFracLocal + nbWellInjLocal
       start = start + 1
    end do

    ! Sort JacA%Num
    do i=1,JacA%Nb

       do k=1, JacA%Pt(i+1)-JacA%Pt(i)
          do j=JacA%Pt(i)+1, JacA%Pt(i+1)-k

             if( JacA%Num(j) > JacA%Num(j+1) ) then
                tmp = JacA%Num(j)
                JacA%Num(j) = JacA%Num(j+1)
                JacA%Num(j+1) = tmp
             end if

          end do
       end do

    end do

    ! allocate JacA%Val
    allocate( JacA%Val( NbCompThermique, NbCompThermique, Nz)) ! number of non zero
!    JacA%Val(:,:,:) = 0.d0

  end subroutine Jacobian_StrucJacA


  ! free
  subroutine Jacobian_free

    deallocate(JacBigA%Pt)
    deallocate(JacBigA%Num)
    deallocate(JacBigA%Val)

    deallocate(JacA%Pt)
    deallocate(JacA%Num)
    deallocate(JacA%Val)

    deallocate(bigSm)
    deallocate(Sm)

    deallocate(csrK)
    deallocate(csrSR)

  end subroutine Jacobian_free


  ! num transformation, used for TkLocal(k)%(i,j), where i, j is node or frac
  ! from i, j to its num (local)
  ! k: loop of cell
  subroutine Jacobian_RowCol_KSR(k, &
       nbNodeCell, nbFracCell, &
       rowk,  colk, &
       rowSR, colSR)

    integer, intent(in) :: k, &
         nbNodeCell, nbFracCell

    integer, intent(out) :: &
         rowk, colk, &
         rowSR( NbNodeCellMax+NbFracCellMax), &
         colSR( NbNodeCellMax+NbFracCellMax)

    integer :: i, in

    ! rowk
    rowk = k &
         + NbNodeOwn_Ncpus(commRank+1) &
         + NbFracOwn_Ncpus(commRank+1)

    ! colk
    colk = k &
         + NbNodeLocal_Ncpus(commRank+1) &
         + NbFracLocal_Ncpus(commRank+1)

    ! colSR, nodes
    do i=1, nbNodeCell
       colSR(i) = NodebyCellLocal%Num(i + NodebyCellLocal%Pt(k))
    end do

    ! colSR, frac
    do i=1, nbFracCell


       ! in is face number, need to be transformed to number of frac using FaceToFracLocal
       in = FracbyCellLocal%Num(i + FracbyCellLocal%Pt(k))

       colSR(nbNodeCell+i) = FaceToFracLocal(in) &
            + NbNodeLocal_Ncpus(commRank+1)
    end do

    rowSR(:) = 0

    ! rowSR, nodes
    do i=1, nbNodeCell

       ! node
       ! Rq. if i is dir Darcy, rowSR(i) will not be used for Darcy
       !     if i is dir Fourier, rowSR(i) will not be used for Fourier
       rowSR(i) = NodebyCellLocal%Num(i + NodebyCellLocal%Pt(k))
    end do

    ! rowSR, frac
    do i=1,nbFracCell

       in = FracbyCellLocal%Num(i + FracbyCellLocal%Pt(k))
       if(in <= NbFaceOwn_Ncpus(commRank+1) ) then ! own, rq: in is face number, not frac number
          rowSR(nbNodeCell+i) = FaceToFracLocal(in) &
               + NbNodeOwn_Ncpus(commRank+1)
       end if
    end do

  end subroutine Jacobian_RowCol_KSR


  ! num transformation, used for TkFracLocal(k)%(i,j), where i, j is node
  ! from i, j to its num (local)
  ! k: loop of frac
  subroutine Jacobian_RowCol_FR(k, &
       nbNodeFrac, &
       rowk, colk, &
       rowSR, colSR)

    integer, intent(in) :: k, nbNodeFrac

    integer, intent(out) :: &
         rowk, colk, &
         rowSR( NbNodeFaceMax), &
         colSR( NbNodeFaceMax)

    integer :: i, fk

    fk = FracToFaceLocal(k) ! this frac is which face

    ! rowk
    rowk = k + NbNodeOwn_Ncpus(commRank+1)

    ! colk
    colk = k + NbNodeLocal_Ncpus(commRank+1)

    ! rowR
    rowSR(:) = 0
    do i=1, nbNodeFrac
       rowSR(i) = NodebyFaceLocal%Num(i + NodebyFaceLocal%Pt(fk))
    end do

    ! colSR = rowR
    colSR(:) = rowSR(:)
    ! colSR(:)  = 0
    ! do i=1, nbNodeFrac
    !    colSR(i) = NodebyFaceLocal%Num(i + NodebyFaceLocal%Pt(fk))
    ! end do

  end subroutine Jacobian_RowCol_FR

  subroutine Jacobian_GetSolCell_C(increments_pointers) &
     bind(C, name="Jacobian_GetSolCell")

    type(Newton_increments_pointers), intent(in), value :: increments_pointers
    type(Newton_increments) :: increments
    
    call Newton_pointers_to_values(increments_pointers, increments)
    call Jacobian_GetSolCell( &
       increments%nodes, increments%fractures, increments%cells &
    )

  end subroutine Jacobian_GetSolCell_C

  ! A3 * IncrementNodeFrac + A4 * IncrementCell = bigSm (cell part)
  ! A4 is a block diag matrix,
  ! its block element has been done LU factorization
  subroutine Jacobian_GetSolCell( &
       NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell &
    )

    real(c_double), dimension(:,:), intent(in) :: &
         NewtonIncreNode, &
         NewtonIncreFrac

    real(c_double), dimension(:,:), intent(out) :: &
         NewtonIncreCell

    integer :: k, rowk, s, cols, nums, i, j, nz
    integer :: info, errcode, Ierr
    integer :: NbNodeFracOwn, NbNodeLocal
    double precision :: tmpInc(NbCompThermique)

    NbNodeFracOwn = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)
    NbNodeLocal = NbNodeLocal_Ncpus(commRank+1)

    do k=1, NbCellLocal_Ncpus(commRank+1)

       rowk = k + NbNodeFracOwn

       tmpInc(:) = 0.d0

       do s=JacbigA%Pt(rowk)+1, JacBigA%Pt(rowk+1)-1  ! A3
          cols = JacbigA%Num(s)

          nz = JacBigA%Pt(rowk+1) ! A4

          ! cols is in part node
          if(cols<=NbNodeLocal) then

             do i=1, NbCompThermique
                do j=1, NbCompThermique

                   tmpInc(i) = tmpInc(i) &
                        + JacBigA%Val(j,i,s) * NewtonIncreNode(j,cols)
                end do
             end do

          else
             ! cols is in part frac
             ! this col corresponts to frac nums
             nums = cols - NbNodeLocal

             do i=1, NbCompThermique
                do j=1, NbCompThermique

                   tmpInc(i) = tmpInc(i) &
                        + JacBigA%Val(j,i,s) * NewtonIncreFrac(j,nums)
                end do
             end do
          end if
       end do

       NewtonIncreCell(1:NbCompThermique,k) = &
            bigSm(1:NbCompThermique,rowk) - tmpInc(1:NbCompThermique)

       call dgetrs('T', NbCompThermique, 1, &
            JacBigA%Val(:,:,nz), NbCompThermique, ipiv(:,k), &
            NewtonIncreCell(1:NbCompThermique,k), NbCompThermique, info)

       if(info /= 0) then
          write(0,'(A,I5)', advance='no') "dgetrs error in inverse Schur complement k =", k
          write(0,'(A,I5)') ",  info =", info
          call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
       end if

       ! if(commRank==1 .and. k==1) then
       !    print*, ""
       !    write(*,'(ES22.13)') NewtonIncreCell(1:NbCompThermique,k)
       ! end if

    end do

  end subroutine Jacobian_GetSolCell

end module Jacobian


! nz = JacBigA%Pt(rowk) + csrK(colk)
! cols = colSR(s)
! nz = JacBigA%Pt(rowk) + csrK(cols)
! if(k==1 .and. s==1 .and. commRank==0) then

!    ! print*, NumIncTotalPrimCell(:,1)
!    ! print*, NumIncTotalPrimNode(:,1)

!    ! do i=1, NbComp
!    !    do j=1, NbIncTotalPrim_ctx(3)
!    !       print*, JacBigA%Val(j,i,nz)+divS1(j,i) + divS2(j,i) + divR2(j,i,s)
!    !    end do
!    !    print*, ""
!    ! end do

!    do j=1, NbIncTotalPrim_ctx(3)
!       ! print*, JacBigA%Val(j,NbComp+1,nz) + divFourierFlux_k(j)+ divEgK(j)
!       print*, JacBigA%Val(j,NbComp+1,nz) + divEgR(j,s) + divFourierFlux_r(j,s) + divEgS(j)
!    end do

! end if

! do i=JacA%Pt(rowk)+1, JacA%Pt(rowk+1)
!    if( JacA%Num(i)==colk) then
!       JacA%Val(:,:,i) = 0.d0
!       do j=1, 5
!          JacA%Val(j,j,i) = 2.d0
!       end do
!    else
!       JacA%Val(:,:,i) = 0.d0
!       do j=1, 5
!          JacA%Val(j,j,i) = 6.d0
!       end do
!    end if
! end do
