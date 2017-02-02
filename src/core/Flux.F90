module Flux

  use MeshSchema

  use IncCV
  use NumbyContext

  use LoisThermoHydro
  use VAGFrac

  implicit none

  ! flux Darcy V_{k,s}^alpha, k is cell/frac
  double precision, allocatable, dimension(:,:,:), public :: &
      FluxDarcyKI, &
      FluxDarcyFI

  ! flux Fourier
  double precision, allocatable, dimension(:,:), protected :: &
      FluxFourierKI, &
      FluxFourierFI

  public :: &
      Flux_allocate, &
      Flux_free,     &
      Flux_DarcyFlux_Cell, &
      Flux_DarcyFlux_Frac, &
      Flux_FourierFlux_Cell, &
      Flux_FourierFlux_Frac

contains

  ! allocate
  subroutine Flux_allocate

    ! flux
    allocate( FluxDarcyKI &
        (NbPhase, NbNodeCellMax+NbFracCellMax, NbCellLocal_Ncpus(commRank+1)) )
    allocate( FluxDarcyFI &
        (NbPhase, NbNodeFaceMax, NbFracLocal_Ncpus(commRank+1)) )

#ifdef _THERMIQUE_

    allocate( FluxFourierKI &
        (NbNodeCellMax+NbFracCellMax, NbCellLocal_Ncpus(commRank+1)) )
    allocate( FluxFourierFI &
        (NbNodeFaceMax, NbFracLocal_Ncpus(commRank+1)))

#endif

  end subroutine Flux_allocate

  ! free
  subroutine Flux_free

    deallocate(FluxDarcyKI)
    deallocate(FluxDarcyFI)

#ifdef _THERMIQUE_

    deallocate(FluxFourierKI)
    deallocate(FluxFourierFI)
#endif

  end subroutine Flux_free

  ! Structure of this subroutine:
  ! loop of cell k
  !   a. loop of node i of cell k
  !       1. compute rho_ki_alpha
  !          loops of Q_k and Q_i

  !       2. loop of node j
  !          loops of Q_k and Q_i

  !       3. loop of frac j
  !          loops of Q_k and Q_i

  !   b. loop of frac i of cell k
  !       1. compute rho_ki_alpha
  !          loops of Q_k and Q_i

  !       2. loop of node j
  !          loops of Q_k and Q_i

  !       3. loop of frac j
  !          loops of Q_k and Q_i
  subroutine Flux_DarcyFlux_Cell

    integer :: k, i, j, fj, fi
    integer :: numi, numj, nph_i, nph_k, numph_i, numph_k
    integer :: NbNodeCell, NbFracCell

    double precision :: rho_ki_alpha(NbPhase)
    logical :: Id_Qki(NbPhase)
    double precision :: Tkij, Pkj, zkj, Satki

    FluxDarcyKI(:,:,:) = 0.d0

    ! FluxDarcyKI
    do k=1, NbCellLocal_Ncpus(commRank+1) ! loop of cell

      ! number of nodes/fracs in cell k
      NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
      NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

      ! i is node. Two parts:
      !   1. compute rho_ki^alpha
      !   2. loop of j
      do i=1, NbNodeCell

        numi = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+i) ! num of node i

        ! compute rho_ki^alpha: loop of Q_k and loop of Q_i

        rho_ki_alpha(:) = 0.d0

        Id_Qki(:) = .false.

        do nph_k = 1, NbPhasePresente_ctx( IncCell(k)%ic) ! phases present: Q_k
          numph_k = NumPhasePresente_ctx(nph_k, IncCell(k)%ic)

          Id_Qki(numph_k) = .true.

          ! Satki = IncCell(k)%Saturation(numph_k) + IncNode(numi)%Saturation(numph_k) ! S_k^alpha+S_i^alpha

          rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
              0.5d0 * (DensiteMassiqueCell(numph_k, k) + DensiteMassiqueNode(numph_k, numi))

          ! if(abs(Satki)<eps) then
          !    rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
          !         0.5d0 * (DensiteMassiqueCell(numph_k, k) + DensiteMassiqueNode(numph_k, numi))
          ! else
          !    rho_ki_alpha(numph_k) = & 
          !         (IncCell(k)%Saturation(numph_k) * DensiteMassiqueCell(numph_k, k) & 
          !         + IncNode(numi)%Saturation(numph_k) * DensiteMassiqueNode(numph_k, numi))/Satki
          ! end if
        end do

        do nph_i = 1, NbPhasePresente_ctx( IncNode(numi)%ic) ! phases present: Q_i
          numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)

          if(Id_Qki(numph_i) .eqv. .false.) then ! not computed in the loop of Q_k

            ! Satki = IncCell(k)%Saturation(numph_i) + IncNode(numi)%Saturation(numph_i) ! S_k^alpha+S_i^alpha

            rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
                0.5d0 * (DensiteMassiqueCell(numph_i, k) + DensiteMassiqueNode(numph_i, numi))

            ! if(abs(Satki)<eps) then
            !    rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
            !         0.5d0 * (DensiteMassiqueCell(numph_i, k) + DensiteMassiqueNode(numph_i, numi))
            ! else
            !    rho_ki_alpha(numph_i) = & 
            !         (IncCell(k)%Saturation(numph_i) * DensiteMassiqueCell(numph_i, k) & 
            !         + IncNode(numi)%Saturation(numph_i) * DensiteMassiqueNode(numph_i, numi))/Satki
            ! end if
          end if
        end do

        ! i is node, j is node
        do j=1, NbNodeCell

          numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+j) ! num of node j

          Tkij = TkLocal_Darcy(k)%Pt(i,j) ! a_{k,s}^s'
          Pkj = IncCell(k)%Pression - IncNode(numj)%Pression ! P_k - P_s'
          zkj  = Gravite*(XCellLocal(3,k) - XNodeLocal(3,numj)) ! g*(z_k - z_s')

          ! loop of Q_k
          Id_Qki(:) = .false.

          do nph_k=1, NbPhasePresente_ctx(IncCell(k)%ic)
            numph_k = NumPhasePresente_ctx(nph_k, IncCell(k)%ic)

            Id_Qki(numph_k) = .true. ! this phase is in Q_k, not need to consider in the loop of Q_i

            FluxDarcyKI(numph_k,i,k) = FluxDarcyKI(numph_k,i,k) &
                + (Pkj + PressionCapCell(nph_k,k) - PressionCapNode(nph_k,numj) &
                + rho_ki_alpha(numph_k)*zkj) * Tkij
          end do

          ! loop of Q_i
          do nph_i=1, NbPhasePresente_ctx(IncNode(numi)%ic)
            numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)

            if( Id_Qki(numph_i) .eqv. .false.) then ! this phase is not in Q_k

              FluxDarcyKI(numph_i,i,k) = FluxDarcyKI(numph_i,i,k) &
                  + (Pkj + PressionCapCell(nph_i,k) - PressionCapNode(nph_i,numj) &
                  + rho_ki_alpha(numph_i)*zkj) * Tkij
            end if
          end do

        end do

        ! i is node, j is frac
        do j=1, NbFracCell

          fj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+j) 
          numj = FaceToFracLocal(fj) ! fj is face number, numj is frac number

          Tkij = TkLocal_Darcy(k)%Pt(i,j+NbNodeCell) ! a_{k,s}^s'
          Pkj = IncCell(k)%Pression - IncFrac(numj)%Pression ! P_k - P_s'
          zkj  = Gravite*(XCellLocal(3,k) - XFaceLocal(3,fj)) ! g*(z_k - z_s')

          ! loop of Q_k
          Id_Qki(:) = .false.

          do nph_k=1, NbPhasePresente_ctx(IncCell(k)%ic)
            numph_k = NumPhasePresente_ctx(nph_k, IncCell(k)%ic)

            Id_Qki(numph_k) = .true. ! this phase is in Q_k, not need to consider in the loop of Q_i

            FluxDarcyKI(numph_k,i,k) = FluxDarcyKI(numph_k,i,k) &
                + (Pkj + PressionCapCell(nph_k,k) - PressionCapFrac(nph_k,numj) &
                + rho_ki_alpha(numph_k)*zkj) * Tkij
          end do

          ! loop of Q_i, ps. not Q_j !
          do nph_i=1, NbPhasePresente_ctx(IncNode(numi)%ic)
            numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)

            if( Id_Qki(numph_i) .eqv. .false.) then ! this phase is not in Q_k

              FluxDarcyKI(numph_i,i,k) = FluxDarcyKI(numph_i,i,k) &
                  + (Pkj + PressionCapCell(nph_i,k) - PressionCapFrac(nph_i,numj) &
                  + rho_ki_alpha(numph_i)*zkj) * Tkij
            end if
          end do

        end do ! end of j
      end do


      ! i is frac. Two parts:
      !   1. compute rho_ki^alpha
      !   2. loop of j
      do i=1, NbFracCell

        fi = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+i) ! fi is face number
        numi = FaceToFracLocal(fi) ! numi is frac number

        rho_ki_alpha(:) = 0.d0

        Id_Qki(:) = .false.

        do nph_k = 1, NbPhasePresente_ctx( IncCell(k)%ic) ! phases present: Q_k
          numph_k = NumPhasePresente_ctx( nph_k, IncCell(k)%ic)

          Id_Qki(numph_k) = .true.

          ! Satki = IncCell(k)%Saturation(numph_k) + IncFrac(numi)%Saturation(numph_k) ! S_k^alpha+S_i^alpha

          rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
              0.5d0 * (DensiteMassiqueCell(numph_k, k) + DensiteMassiqueFrac(numph_k, numi))

          ! if(abs(Satki)<eps) then
          !   rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
          !         0.5d0 * (DensiteMassiqueCell(numph_k, k) + DensiteMassiqueFrac(numph_k, numi))
          ! else
          !    rho_ki_alpha(numph_k) = & 
          !         (IncCell(k)%Saturation(numph_k) * DensiteMassiqueCell(numph_k, k) & 
          !         + IncFrac(numi)%Saturation(numph_k) * DensiteMassiqueFrac(numph_k, numi))/Satki
          ! end if
        end do

        do nph_i = 1, NbPhasePresente_ctx( IncFrac(numi)%ic) ! phases present: Q_i
          numph_i = NumPhasePresente_ctx( nph_i, IncFrac(numi)%ic)

          if(Id_Qki(numph_i) .eqv. .false.) then

            ! Satki = IncCell(k)%Saturation(numph_i) + IncFrac(numi)%Saturation(numph_i) ! S_k^alpha+S_i^alpha

            rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
                0.5d0 * (DensiteMassiqueCell(numph_i, k) + DensiteMassiqueFrac(numph_i, numi))

            ! if(abs(Satki)<eps) then
            !    rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
            !         0.5d0 * (DensiteMassiqueCell(numph_i, k) + DensiteMassiqueFrac(numph_i, numi))
            ! else
            !    rho_ki_alpha(numph_i) = & 
            !         (IncCell(k)%Saturation(numph_i) * DensiteMassiqueCell(numph_i, k) & 
            !         + IncFrac(numi)%Saturation(numph_i) * DensiteMassiqueFrac(numph_i, numi))/Satki
            ! end if

          end if
        end do


        ! i is frac, j is node
        do j=1, NbNodeCell

          numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+j)

          Tkij = TkLocal_Darcy(k)%Pt(i+NbNodeCell,j)            ! a_{k,s}^s'
          Pkj  = IncCell(k)%Pression - IncNode(numj)%Pression   ! P_k - P_s'
          zkj  = Gravite*(XCellLocal(3,k) - XNodeLocal(3,numj)) ! g*(z_k - z_s')

          ! loop of Q_k
          Id_Qki(:) = .false.

          do nph_k=1, NbPhasePresente_ctx(IncCell(k)%ic)
            numph_k = NumPhasePresente_ctx(nph_k, IncCell(k)%ic)

            Id_Qki(numph_k) = .true. ! this phase is in Q_k, not need to consider in the loop of Q_i

            FluxDarcyKI(numph_k,i+NbNodeCell,k) = FluxDarcyKI(numph_k,i+NbNodeCell,k) &
                + (Pkj + PressionCapCell(nph_k,k) - PressionCapNode(nph_k,numj) &
                + rho_ki_alpha(numph_k)*zkj) * Tkij
          end do

          ! loop of Q_i, Ps. not Q_j !
          do nph_i=1, NbPhasePresente_ctx(IncFrac(numi)%ic)
            numph_i = NumPhasePresente_ctx(nph_i, IncFrac(numi)%ic)

            if( Id_Qki(numph_i) .eqv. .false.) then ! this phase is not in Q_k

              FluxDarcyKI(numph_i,i+NbNodeCell,k) = FluxDarcyKI(numph_i,i+NbNodeCell,k) &
                  + (Pkj + PressionCapCell(nph_i,k) - PressionCapNode(nph_i,numj) &
                  + rho_ki_alpha(numph_i)*zkj) * Tkij
            end if
          end do

          ! if(commRank==1 .and. k==1) then
          !    print*, "j node", j, Pkj, zkj, (Pkj+rho_ki_alpha(2)*zkj), Tkij
          ! end if

        end do ! end of j

        ! i is frac, j is frac
        do j=1, NbFracCell

          fj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+j) 
          numj = FaceToFracLocal(fj) ! fj is face number, numj is frac number

          Tkij = TkLocal_Darcy(k)%Pt(i+NbNodeCell,j+NbNodeCell) ! a_{k,s}^s'
          Pkj = IncCell(k)%Pression - IncFrac(numj)%Pression    ! P_k - P_s'
          zkj  = Gravite*(XCellLocal(3,k) - XFaceLocal(3,fj))   ! g*(z_k - z_s')

          ! loop of Q_k
          Id_Qki(:) = .false.

          do nph_k=1, NbPhasePresente_ctx(IncCell(k)%ic)
            numph_k = NumPhasePresente_ctx(nph_k, IncCell(k)%ic)

            Id_Qki(numph_k) = .true. ! this phase is in Q_k, not need to consider in the loop of Q_i

            FluxDarcyKI(numph_k,i+NbNodeCell,k) = FluxDarcyKI(numph_k,i+NbNodeCell,k) &
                + (Pkj + PressionCapCell(nph_k,k) - PressionCapFrac(nph_k,numj) &
                + rho_ki_alpha(numph_k)*zkj) * Tkij
          end do

          ! loop of Q_i, Ps. not Q_j !
          do nph_i=1, NbPhasePresente_ctx(IncFrac(numi)%ic)
            numph_i = NumPhasePresente_ctx(nph_i, IncFrac(numi)%ic)

            if( Id_Qki(numph_i) .eqv. .false.) then ! this phase is not in Q_k

              FluxDarcyKI(numph_i,i+NbNodeCell,k) = FluxDarcyKI(numph_i,i+NbNodeCell,k) &
                  + (Pkj + PressionCapCell(nph_i,k) - PressionCapFrac(nph_i,numj) &
                  + rho_ki_alpha(numph_i)*zkj) * Tkij
            end if
          end do

        end do ! fin loop of j

      end do ! fin loop of i in cell k

    end do ! fin loop cell k

    ! if(commRank==1) then
    !    do i=1, 9
    !       print*, i, FluxDarcyKI(:,i,1)
    !    end do
    ! end if

  end subroutine Flux_DarcyFlux_Cell


  ! Structure of this subroutine:
  ! loop of frac k
  !   a. loop of node i of frac k
  !       1. compute rho_ki_alpha
  !          loops of Q_k and Q_i

  !       2. loop of node j
  !          loops of Q_k and Q_i
  subroutine Flux_DarcyFlux_Frac

    integer :: k, fk, i, j, numi, numj
    integer :: nph_i, numph_i, nph_k, numph_k
    integer :: NbNodeFrac

    double precision :: Pkj, Satki, Tkij, zkj

    double precision :: rho_ki_alpha(NbPhase)
    logical :: Id_Qki(NbPhase)

    FluxDarcyFI(:,:,:) = 0

    ! loop of frac
    do k=1, NbFracLocal_Ncpus(commRank+1)

      fk = FracToFaceLocal(k) ! fk is face number

      ! number of nodes in a frac
      NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

      do i=1, NbNodeFrac

        numi = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+i) ! node number of i

        ! compute rho_ki^alpha: loop of Q_k and loop of Q_i
        rho_ki_alpha(:) = 0.d0

        Id_Qki(:) = .false.

        do nph_k=1, NbPhasePresente_ctx( IncFrac(k)%ic) ! phases present: Q_k
          numph_k = NumPhasePresente_ctx(nph_k, IncFrac(k)%ic)

          Id_Qki(numph_k) = .true.

          ! Satki = IncFrac(k)%Saturation(numph_k) + IncNode(numi)%Saturation(numph_k) ! S_k^alpha+S_i^alpha

          rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
              0.5d0 * (DensiteMassiqueFrac(numph_k, k) + DensiteMassiqueNode(numph_k, numi) )

          ! if(abs(Satki)<eps) then
          !    rho_ki_alpha(numph_k) = & ! Attention: numph_k used for densitemassique
          !         0.5d0 * (DensiteMassiqueFrac(numph_k, k) + DensiteMassiqueNode(numph_k, numi) )
          ! else
          !    rho_ki_alpha(numph_k) = & 
          !         (IncFrac(k)%Saturation(numph_k) * DensiteMassiqueFrac(numph_k, k) & 
          !         + IncNode(numi)%Saturation(numph_k) * DensiteMassiqueNode(numph_k, numi))/Satki
          ! end if
        end do

        do nph_i=1,  NbPhasePresente_ctx( IncNode(numi)%ic) ! phases present: Q_i
          numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)

          if( Id_Qki(numph_i) .eqv. .false.) then ! not computed in loop of Q_k

            ! Satki = IncFrac(k)%Saturation(numph_i) + IncNode(numi)%Saturation(numph_i) ! S_k^alpha+S_i^alpha

            rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
                0.5d0 * (DensiteMassiqueFrac(numph_i, k) + DensiteMassiqueNode(numph_i, numi))

            ! if(abs(Satki)<eps) then
            !    rho_ki_alpha(numph_i) = & ! Attention: numph_k used for densitemassique
            !         0.5d0 * (DensiteMassiqueFrac(numph_i, k) + DensiteMassiqueNode(numph_i, numi))
            ! else
            !    rho_ki_alpha(numph_i) = & 
            !         (IncFrac(k)%Saturation(numph_i) * DensiteMassiqueFrac(numph_i, k) & 
            !         + IncNode(numi)%Saturation(numph_i) * DensiteMassiqueNode(numph_i, numi))/Satki
            ! end if

          end if
        end do

        ! if(commRank==1 .and. k==1) then
        !    print*, i, rho_ki_alpha(:)
        ! end if


        do j=1, NbNodeFrac
          numj = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk)+j)

          Tkij = TkFracLocal_Darcy(k)%Pt(i,j) ! a_{k,s}^s'
          Pkj = IncFrac(k)%Pression - IncNode(numj)%Pression ! P_k - P_s'
          zkj  = Gravite*(XFaceLocal(3,fk) - XNodeLocal(3,numj)) ! g*(z_k - z_s')

          ! loop of Q_k
          Id_Qki(:) = .false.

          do nph_k=1, NbPhasePresente_ctx(IncFrac(k)%ic)
            numph_k = NumPhasePresente_ctx(nph_k, IncFrac(k)%ic)

            Id_Qki(numph_k) = .true. ! this phase is in Q_k, not need to consider in the loop of Q_i

            FluxDarcyFI(numph_k,i,k) = FluxDarcyFI(numph_k,i,k) &
                + (Pkj + PressionCapFrac(nph_k,k) - PressionCapNode(nph_k,numj) &
                + rho_ki_alpha(numph_k)*zkj) * Tkij

            ! if(commRank==1 .and. k==1) then
            !    print*, k, i, j, numph_k, FluxDarcyFI(numph_k,i,k), (Pkj+rho_ki_alpha(numph_k)*zkj), Tkij
            ! end if

          end do

          ! loop of Q_i
          do nph_i=1, NbPhasePresente_ctx(IncNode(numi)%ic)
            numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)

            if(Id_Qki(numph_i) .eqv. .false.) then ! this phase is not in Q_k

              FluxDarcyFI(numph_i,i,k) = FluxDarcyFI(numph_i,i,k) &
                  + (Pkj + PressionCapFrac(nph_i,k) - PressionCapNode(nph_i,numj) &
                  + rho_ki_alpha(numph_i)*zkj) * Tkij
            end if
          end do

        end do
      end do

    end do ! end loop of frac k


  end subroutine Flux_DarcyFlux_Frac


  subroutine Flux_FourierFlux_Cell

    integer :: k, i, j, fj
    integer :: numj
    integer :: NbNodeCell, NbFracCell

    FluxFourierKI(:,:) = 0.d0

    ! FluxFourierKI
    do k=1, NbCellLocal_Ncpus(commRank+1) ! loop of cell

      ! number of nodes/fracs in cell k
      NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
      NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

      ! i is node
      do i=1, NbNodeCell

        ! i is node, j is node
        do j=1, NbNodeCell

          numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+j)

          FluxFourierKI(i,k) = FluxFourierKI(i,k) &
              + TkLocal_Fourier(k)%Pt(i,j) * &
              (IncCell(k)%Temperature - IncNode(numj)%Temperature)
        end do

        ! i is node, j is frac
        do j=1, NbFracCell

          fj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+j) 
          numj = FaceToFracLocal(fj) ! fj is face number

          FluxFourierKI(i,k) = FluxFourierKI(i,k) &
              + TkLocal_Fourier(k)%Pt(i,j+NbNodeCell)&
              *(IncCell(k)%Temperature - IncFrac(numj)%Temperature)
        end do
      end do

      ! i is frac
      do i=1, NbFracCell

        ! i is frac, j is node
        do j=1, NbNodeCell

          numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k)+j)

          FluxFourierKI(i+NbNodeCell,k) = FluxFourierKI(i+NbNodeCell,k) &
              + TkLocal_Fourier(k)%Pt(i+NbNodeCell,j) &
              *(IncCell(k)%Temperature - IncNode(numj)%Temperature)
        end do

        ! i is frac, j is frac
        do j=1, NbFracCell

          fj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k)+j) 
          numj = FaceToFracLocal(fj) ! fj is face number

          FluxFourierKI(i+NbNodeCell,k) = FluxFourierKI(i+NbNodeCell,k) &
              + TkLocal_Fourier(k)%Pt(i+NbNodeCell,j+NbNodeCell) &
              *(IncCell(k)%Temperature - IncFrac(numj)%Temperature)
        end do
      end do

    end do ! fin loop cell: k

  end subroutine Flux_FourierFlux_Cell


  subroutine Flux_FourierFlux_Frac

    integer :: k, fk, i, j, numi, numj
    integer :: NbNodeFrac

    FluxFourierFI(:,:) = 0

    ! loop of frac
    do k=1, NbFracLocal_Ncpus(commRank+1)

      fk = FracToFaceLocal(k) ! fk is face number

      ! number of nodes in a frac
      NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

      do i=1, NbNodeFrac

        do j=1, NbNodeFrac
          numj = NodebyFaceLocal%Num( NodebyFaceLocal%Pt(fk)+j)

          FluxFourierFI(i,k) = FluxFourierFI(i,k) &
              + TkFracLocal_Fourier(k)%pt(i,j) &
              *(IncFrac(k)%Temperature - IncNode(numj)%Temperature)
        end do
      end do

    end do ! end loop of frac k


  end subroutine Flux_FourierFlux_Frac

end module Flux