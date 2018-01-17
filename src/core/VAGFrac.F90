! Three levels node numerotation
!    num (local) : num in this proc
!    num (cell) : num in a cell
!    num (face) : num in a face

! Unknown = (node, frac)
!    in a face, size is nbNodeFace+1, only 1 frac in a face frac
!    in a cell, size is nbNodeCell+nbFracCell

module VAGFrac

  use CommonMPI
  use MeshSchema
  use DefModel

  implicit none

  type Type_Trans
    double precision, dimension(:,:), pointer :: pt
  end type Type_Trans

  ! Outputs: Trans...
  type(Type_Trans), dimension(:), allocatable, protected :: &
      TkLocal_Darcy ,    & ! transmissivites darcy on matrix
      TkFracLocal_Darcy, & ! transmissivites darcy on frac
      TkLocal_Fourier ,  & ! transmissivites fourier on matrix
      TkFracLocal_Fourier  ! transmissivites fourier on frac

  ! Volume
  ! = sum_{M_s} alpha_{k,s} vol_K ...
  double precision, allocatable, dimension(:), protected :: &
      VolDarcyCell, &
      VolDarcyFrac, &
      VolDarcyNode

  ! Porous volume Darcy
  ! = sum_{M_s} alpha_{k,s} phi * vol_K ...
  double precision, allocatable, dimension(:), protected :: &
      PoroVolDarcyCell, &
      PoroVolDarcyFrac, &
      PoroVolDarcyNode

#ifdef _THERMIQUE_

  ! Two porous volume Fourier
  ! = sum_{M_s} alpha_{k,s} phi * vol_K ...
  ! = sum_{M_s} alpha_{k,s} (1-phi) * vol_K ...
  double precision, allocatable, dimension(:), protected :: &
      PoroVolFourierCell, & 
      PoroVolFourierFrac, &
      PoroVolFourierNode

  double precision, allocatable, dimension(:), protected :: &
      Poro_1volFourierCell, & 
      Poro_1volFourierFrac, &
      Poro_1volFourierNode  

  double precision, allocatable, dimension(:), protected :: &
      CellThermalSourceVol, & 
      FracThermalSourceVol, & 
      NodeThermalSourceVol
#endif


  integer, allocatable, dimension(:), private :: &
      UnkFaceToUnkCell ! from node in num (face) to node in num (cell)

  public :: &
      VAGFrac_TransDarcy,     &
      VAGFrac_VolsDarcy,   &
      VAGFrac_Free

#ifdef _THERMIQUE_

  public :: &
      VAGFrac_TransFourier,   &
      VAGFrac_VolsFourier
#endif

  private:: &
      VAGFrac_Trans,     & ! comput trans
      VAGFrac_TransFrac, & ! comput transfrac
      VAGFrac_UnkFaceToUnkCell, &
      VAGFrac_VolTetra, &
      VAGFrac_VecNormalT, & ! vec normal
      VAGFrac_InvA

contains


  ! compute Trans and TransFrac for Darcy
  subroutine VAGFrac_TransDarcy(Perm, PermFrac)

    double precision, dimension(:,:,:), intent(in) :: Perm
    double precision, dimension(:), intent(in) :: PermFrac

    call VAGFrac_Trans(Perm, TkLocal_Darcy)
    call VAGFrac_TransFrac(PermFrac, TkFracLocal_Darcy)

  end subroutine VAGFrac_TransDarcy


#ifdef _THERMIQUE_

  ! compute Trans and TransFrac for Fourier
  subroutine VAGFrac_TransFourier(CondThermic, CondThermicFrac)

    double precision, dimension(:,:,:), intent(in) :: CondThermic
    double precision, dimension(:), intent(in) :: CondThermicFrac

    call VAGFrac_Trans(CondThermic, TkLocal_Fourier)
    call VAGFrac_TransFrac(CondThermicFrac, TkFracLocal_Fourier)

  end subroutine VAGFrac_TransFourier
#endif

  ! called by VAGFrac_TransDarcy/Fourier
  subroutine VAGFrac_Trans(lamda, TkLocal)

    ! lamda is permeability for darcy
    ! lamda is thermal conductivity for fourier
    double precision, dimension(:,:,:),intent(in) :: lamda

    ! output
    type(Type_Trans), dimension(:), allocatable, intent(out) :: TkLocal

    ! calcul du schema en ne remplissant que les noeuds des faces pour les gk
    double precision, allocatable, dimension(:,:) :: GkT

    ! tmp
    integer :: &
        in1, in2, inf, & ! num (face) in Unknown in a face
        n1, n2 ! num (local) in Unknown in a cell

    integer :: i,j,k,m,mm,nn,ipt,& 
        in,k1,k2,kn1,kn2

    double precision, dimension(3) :: &
        xk, xf, x1, x2, xt, & ! cordinate
        v12k, v1fk, vf2k, v   ! normal directive

    double precision :: &
        ss

    double precision :: surf, volT

    integer :: nbNodeCell, nbNodeFace, nbFracCell, nbFracFace, nbUnkFace

    ! tmp values to simply notations
    integer :: nbCellLocal, nbFaceLocal, nbNodeLocal, nbFracLocal

    nbCellLocal = NbCellLocal_Ncpus(commRank+1)
    nbFaceLocal = NbFaceLocal_Ncpus(commRank+1)
    nbNodeLocal = NbNodeLocal_Ncpus(commRank+1)
    nbFracLocal = NbFracLocal_Ncpus(commRank+1)

    ! TkLocal
    allocate( TkLocal(nbCellLocal) )

    ! GkT
    allocate( GkT(3,nbNodeFaceMax+1) ) ! +1 for frac

    allocate( UnkFaceToUnkCell(nbNodeFaceMax+1) )     

    ! three levels loops
    !     mailles: k
    !     face   : i
    !     node   : 

    ! nbNodeCell: number of nodes in cell k; scalar, defined within the loop
    ! nbFracCell: number of face frac in cell k; scalar, defined within the loop
    ! nbNodeFace: number of nodes in face i; scalar, defined within the loop

    ! boucle sur les mailles k own du proc 
    do k=1,nbCellLocal

      ! cordinate center of cell
      xk(:) = XCellLocal(:,k)

      ! number of node, frac in cell k
      nbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
      nbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

      ! allocate TkLocal(k) and initial
      allocate( TkLocal(k)%pt(nbNodeCell+nbFracCell,nbNodeCell+nbFracCell) )
      TkLocal(k)%pt(:,:) = 0.d0

      ! ! cell permeability tensor (cf PhysicsFunction)
      ! rt = RockTypebyCellLocal(k)

      ! loop of face i in mesh k

      ! In a face or in a cell, Unk = (node,frac)
      ! Two num for Unk, in face and in cell, UnkFacetoUnkCell is used to trasform
      do j = FacebyCellLocal%Pt(k)+1, FacebyCellLocal%Pt(k+1)

        i = FacebyCellLocal%Num(j) ! num (local) of face

        ! number of nodes, frac, unknown in face i
        nbNodeFace = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)
        nbFracFace = 1
        nbUnkFace = nbNodeFace + nbFracFace

        ! Use:
        !    NodebyFaceLocal
        !    NodebyCellLocal
        !    FracbyCellLocal
        ! output:
        !    UnkFaceToUnkCell: from Unk (num in face) to Unk (num in cell)
        ! Rq: = 0 if face i is not frac
        call VAGFrac_UnkFaceToUnkCell(k,i)

        inf = nbNodeFace + 1 ! num (face) of frac in a face, in Unk=(node, frac)

        ! isobarycentre de la face
        xf(:) = 0.d0

        do m = NodebyFaceLocal%Pt(i)+1,NodebyFaceLocal%Pt(i+1)
          xf(:) = xf(:) + XNodeLocal(:,NodebyFaceLocal%Num(m))
        enddo
        xf(:) = xf(:)/dble(NodebyFaceLocal%Pt(i+1)-NodebyFaceLocal%Pt(i))

        ! loop of nodes of face i
        do m = NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1)

          ! edge with nodes n1, n2
          ! num (face) of n1 and n2 are in1 and in2            
          in1 = m - NodebyFaceLocal%Pt(i)
          n1 = NodebyFaceLocal%Num(m)

          if (m==NodebyFaceLocal%Pt(i+1)) then 
            n2 = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(i)+1)
            in2 = 1
          else 
            n2 = NodebyFaceLocal%Num(m+1)
            in2 = in1 + 1
          endif

          x1(:) = XNodeLocal(:,n1)
          x2(:) = XNodeLocal(:,n2)

          ! vol of tetra
          call VAGFrac_VolTetra(x1,x2,xf,xk,volT)

          ! centre du tetra 
          xt(:) = (x1(:)+x2(:)+xk(:)+xf(:))/4.d0

          ! vecteur normal a 12k pointant vers X et de norme la surface de 12k 
          call VAGFrac_VecNormalT(x1,x2,xk,xt,v,surf)
          v12k(:) = v(:)*surf

          ! vecteur normal a 1fk pointant vers X et de norme la surface de 1fk 
          call VAGFrac_VecNormalT(x1,xf,xk,xt,v,surf)
          v1fk(:) = v(:)*surf

          ! vecteur normal a f2k pointant vers X et de norme la surface de f2k 
          call VAGFrac_VecNormalT(xf,x2,xk,xt,v,surf)
          vf2k(:) = v(:)*surf

          ! Set GkT=0
          GkT(:,:) = 0.d0

          ! Grad T_12fk = 1/(3*volT) ( v12k (uf-uk) + v1fk (u1-uk) + vf2k (u2-uk) ) 
          ! where uf = 1/NbNodebyFace \sum_(n=node of face) u_n

          if (IdFaceLocal(i) /= -2) then ! this face i (i is num (local) of face) is not a frac

            ss = 1.d0/(3.d0*volT*dble(nbNodeFace))

            ! boucle sur les noeuds de la face pour le terme 1/(3*volT) V12k (uk-uf)
            do ipt=NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1)

              in = ipt - NodebyFaceLocal%Pt(i) ! loop for nodes in num (face)
              GkT(:,in) = v12k(:)*ss
            end do

            ss = 1.d0/(3.d0*volT)

            GkT(:,in1) = GkT(:,in1) + vf2k(:)*ss
            GkT(:,in2) = GkT(:,in2) + v1fk(:)*ss

            ! this face i (i is num (local) of face) is a frac
          else if (IdFaceLocal(i) == -2) then

            ss = 1.d0/(3.d0*volT)

            GkT(:,in1) = vf2k(:)*ss
            GkT(:,in2) = v1fk(:)*ss
            GkT(:,inf) = v12k(:)*ss
          end if

          ! ! Test Gradient in debug mode
          ! call VAGFrac_TestGkT(i, xk, xf, GkT)

          ! calcul des Tkij node
          do k1=1, nbUnkFace
            kn1 = UnkFaceToUnkCell(k1)

            do k2=1, nbUnkFace
              kn2 = UnkFaceToUnkCell(k2)

              if((kn1 /= 0) .and. (kn2 /= 0)) then ! kn1, kn2=0 if i is not frac
                do mm=1,3
                  do nn=1,3
                    TkLocal(k)%pt(kn1,kn2) = TkLocal(k)%pt(kn1,kn2) &
                        + volT*lamda(mm,nn,k)*GkT(mm,k1)*GkT(nn,k2)
                  enddo
                enddo

              end if ! end of if kn1/=0...

            end do ! end of k2
          end do ! end of k1

        enddo ! fin boucle sur les noeuds de la face, index: m          
      enddo ! fin boucle sur les faces de la maille, index: j

    enddo ! fin boucle sur les mailles, index: k

    deallocate(GkT)
    deallocate(UnkFaceToUnkCell)

  end subroutine VAGFrac_Trans


  ! TransFrac
  subroutine VAGFrac_TransFrac(lamdafrac, TkFracLocal)

    ! lamda is permeability for darcy
    ! lamda is thermal conductivity for fourier
    double precision, dimension(:),intent(in) :: lamdafrac

    ! output
    type(Type_Trans), dimension(:), allocatable, intent(out) :: TkFracLocal

    double precision, dimension(3,nbNodeFaceMax) :: GfT

    double precision, dimension(3) :: &
        xf, x1, x2, xt,  & ! cordinate
        v                      ! normal directive

    double precision :: &
        Surf12f     ! surface of a triangle with nodes 1,2 and center of face

    double precision :: &
        AA(3,3), BB(3,3)

    integer :: &
        n1, n2, & ! num (local) of a node
        in1, in2  ! num (face) of a node

    integer :: &
        i, ifrac,    & ! i: loop of face frac, ifrac: num (local) of i
        mm, m, &
        i1, i2

    ! tmp values to simply notations
    integer :: nbCellLocal, nbFracLocal, nbNodeFace

    ! number of frac/cell
    nbFracLocal = NbFracLocal_Ncpus(commRank+1)
    nbCellLocal = NbCellLocal_Ncpus(commRank+1)

    ! TkFracLocal
    allocate( TkFracLocal(nbFracLocal))

    ! GfT = 0
    GfT(:,:) = 0.d0

    ! boucle sur les face frac     
    do ifrac = 1, nbFracLocal

      i = FracToFaceLocal(ifrac)

      ! num of nodes in face i
      nbNodeFace = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)

      ! allocate TkFraclocal and initial
      allocate( TkFracLocal(ifrac)%pt(nbNodeFace, nbNodeFace))
      TkFracLocal(ifrac)%pt(:,:) = 0.d0

      ! isobarycentre de la face 
      xf(:) = XFaceLocal(:,i)

      ! init GfT as zero
      GfT(:,:) = 0.d0

      do m = NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1)

        ! edge of nodes n1, n2
        ! num (face) of n1 and n2 are in1 and in2            
        in1 = m - NodebyFaceLocal%Pt(i)
        n1 = NodebyFaceLocal%Num(m)

        if (m==NodebyFaceLocal%Pt(i+1)) then 
          n2 = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(i)+1)
          in2 = 1
        else 
          n2 = NodebyFaceLocal%Num(m+1)
          in2 = in1 + 1
        endif

        x1(:) = XNodeLocal(:,n1)
        x2(:) = XNodeLocal(:,n2)          
        xt(:) = (x1(:)+x2(:)+xf(:))/3.d0

        call VAGFrac_VecNormalT(x1,x2,xf,xt,v,Surf12f)

        ! Gradient tangentiel 
        AA(1,:) = -x1(:) + xf(:)
        AA(2,:) = -x2(:) + xf(:)
        AA(3,:) = v(:)

        ! BB = inv(AA)
        call VAGFrac_InvA(AA, BB)

        GfT(:,in1) = BB(:,1) 
        GfT(:,in2) = BB(:,2)

        ! ! Test Gradient tangentiel dans mode de debug
        ! call VAGFrac_TestGfT(i, x1, x2, xf, v, GfT)

        ! TkFraclocal
        ! Rq: GfT(:,j)=0 if j is not in1, in2
        do i1=1, nbNodeFace
          if((i1==in1) .or. (i1==in2)) then

            do i2=1, nbNodeFace
              if((i2==in1) .or. (i2==in2)) then

                do mm=1,3
                  TkFracLocal(ifrac)%pt(i1,i2) = TkFracLocal(ifrac)%pt(i1,i2) &
                      + Thickness*Surf12f*lamdafrac(ifrac)*GfT(mm,i1)*GfT(mm,i2)
                end do

              end if ! end if of i2=in1..
            end do ! end of i2

          end if ! end if of i1=in1...
        end do ! end of i1

        ! Reset GfT=0 for next iteration
        GfT(:,in1) = 0.d0
        GfT(:,in2) = 0.d0

      end do ! end of loop edge in face       

    end do ! end of loop face

  end subroutine VAGFrac_TransFrac


  ! shift the fraction omega of the cell volume to the nodes according to the
  ! label
  SUBROUTINE VAGFrac_SplitCellVolume( &
      NbCellLocal, &
      CellLabel, &
      NbNodeLocal, &
      NodeLabel, &
      IsVolumeNode, &
      omega, &
      NodebyCellLocal, &
      CellVolume, &
      NodeVolume)

    INTEGER, INTENT(IN) :: NbCellLocal
    INTEGER, INTENT(IN) :: CellLabel(NbCellLocal)
    INTEGER, INTENT(IN) :: NbNodeLocal
    INTEGER, INTENT(IN) :: NodeLabel(NbNodeLocal)
    LOGICAL, INTENT(IN) :: IsVolumeNode(NbNodeLocal)
    DOUBLE PRECISION, INTENT(IN) :: omega
    TYPE(CSR), INTENT(IN) :: NodebyCellLocal

    DOUBLE PRECISION, INTENT(INOUT) :: CellVolume(nbCellLocal)
    DOUBLE PRECISION, INTENT(INOUT) :: NodeVolume(nbNodeLocal)

    INTEGER :: k, ptnumi, numi
    INTEGER :: NbVolume, NbInternalVolume

    DOUBLE PRECISION :: SplitVolume

    DO k=1, nbCellLocal
      NbVolume = 0
      NbInternalVolume = 0

      ! loop of nodes in cell
      DO ptnumi = NodebyCellLocal%Pt(k)+1, NodebyCellLocal%Pt(k+1)
        numi = NodebyCellLocal%Num(ptnumi)
        
        IF(IsVolumeNode(numi))THEN
          NbVolume = NbVolume + 1

          IF(CellLabel(k) == NodeLabel(numi))THEN
            NbInternalVolume = NbInternalVolume + 1
          ENDIF
        ENDIF
      ENDDO

      SplitVolume = omega * CellVolume(k) * NbVolume / NbInternalVolume

      ! loop of nodes in cell
      do ptnumi = NodebyCellLocal%Pt(k)+1, NodebyCellLocal%Pt(k+1)
        numi = NodebyCellLocal%Num(ptnumi)

        IF(IsVolumeNode(numi) .AND. CellLabel(k) == NodeLabel(numi))THEN
          CellVolume(k) = CellVolume(k) - SplitVolume
          NodeVolume(numi) = NodeVolume(numi) + SplitVolume
        end if
      end do
    ENDDO
  END SUBROUTINE VAGFrac_SplitCellVolume


  ! Compute vols darcy:
  !   VolDarcy and PoroVolDarcy
  subroutine VAGFrac_VolsDarcy

    integer :: k, i, ifrac, numi, ptnumi
    integer :: Ierr, errcode
    integer :: NbVolume, NbInternalVolume

    double precision :: s1, s2

    ! VolDarcy
    allocate(VolDarcyCell(NbCellLocal_Ncpus(commRank+1)) )
    allocate(VolDarcyFrac(NbFracLocal_Ncpus(commRank+1)) )
    allocate(VolDarcyNode(NbNodeLocal_Ncpus(commRank+1)) )

    ! PoroVolDarcy
    allocate(PoroVolDarcyCell(NbCellLocal_Ncpus(commRank+1)) )
    allocate(PoroVolDarcyFrac(NbFracLocal_Ncpus(commRank+1)) )
    allocate(PoroVolDarcyNode(NbNodeLocal_Ncpus(commRank+1)) )

    ! Darcy volume
    VolDarcyCell = VolCellLocal
    VolDarcyFrac = Thickness * SurfFracLocal
    VolDarcyNode = 0.d0

    CALL VAGFrac_SplitCellVolume( &
      NbCellLocal_Ncpus(commRank+1), &
      CellRocktypeLocal(1,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(1,:), &
      IdNodeLocal(:)%P /= "d" .AND. IdNodeLocal(:)%Frac == "n", &
      omegaDarcyCell, &
      NodebyCellLocal, &
      VolDarcyCell, &
      VolDarcyNode)

    CALL VAGFrac_SplitCellVolume( &
      NbFracLocal_Ncpus(commRank+1), &
      FracRocktypeLocal(1,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(1,:), &
      IdNodeLocal(:)%P /= "d", &
      omegaDarcyFrac, &
      NodebyFractureLocal, &
      VolDarcyFrac, &
      VolDarcyNode)

    ! Porosity
    PoroVolDarcyCell = PorositeCellLocal * VolCellLocal
    PoroVolDarcyFrac = PorositeFracLocal * Thickness * SurfFracLocal
    PoroVolDarcyNode = 0.d0

    CALL VAGFrac_SplitCellVolume( &
      NbCellLocal_Ncpus(commRank+1), &
      CellRocktypeLocal(1,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(1,:), &
      IdNodeLocal(:)%P /= "d" .AND. IdNodeLocal(:)%Frac == "n", &
      omegaDarcyCell, &
      NodebyCellLocal, &
      PoroVolDarcyCell, &
      PoroVolDarcyNode)

    CALL VAGFrac_SplitCellVolume( &
      NbFracLocal_Ncpus(commRank+1), &
      FracRocktypeLocal(1,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(1,:), &
      IdNodeLocal(:)%P /= "d", &
      omegaDarcyFrac, &
      NodebyFractureLocal, &
      PoroVolDarcyFrac, &
      PoroVolDarcyNode)


    ! check if vol is positive
    do k=1, NbCellLocal_Ncpus(commRank+1)
      if(VolDarcyCell(k)<eps) then
        print*, "vol darcy cell < 0: ", k

        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

    do k=1, NbFracLocal_Ncpus(commRank+1)
      if(VolDarcyFrac(k)<eps) then

        print*, "vol darcy frac < 0: ", k
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

    ! if(commRank==1) then
    !    do k=1, 2
    !       print*, VolCellLocal(k), VolDarcyCell(k)
    !    end do
    ! end if

  end subroutine VAGFrac_VolsDarcy


#ifdef _THERMIQUE_

  ! compute porosite fourier
  subroutine VAGFrac_VolsFourier

    integer :: k, i, j, ifrac, numj, numi
    integer :: Ierr, errcode
    integer :: NbNodeCell

    double precision :: s

    allocate(PoroVolFourierCell(NbCellLocal_Ncpus(commRank+1)))
    allocate(PoroVolFourierFrac(NbFracLocal_Ncpus(commRank+1)))
    allocate(PoroVolFourierNode(NbNodeLocal_Ncpus(commRank+1)))

    allocate(Poro_1volFourierCell(NbCellLocal_Ncpus(commRank+1)))
    allocate(Poro_1volFourierFrac(NbFracLocal_Ncpus(commRank+1)))
    allocate(Poro_1volFourierNode(NbNodeLocal_Ncpus(commRank+1)))

    allocate(CellThermalSourceVol(NbCellLocal_Ncpus(commRank+1)))
    allocate(FracThermalSourceVol(NbFracLocal_Ncpus(commRank+1)))
    allocate(NodeThermalSourceVol(NbNodeLocal_Ncpus(commRank+1)))

    PoroVolFourierCell = PorositeCellLocal * VolCellLocal
    PoroVolFourierFrac = PorositeFracLocal * Thickness * SurfFracLocal
    PoroVolFourierNode = 0.d0

    Poro_1VolFourierCell = (1 - PorositeCellLocal) * VolCellLocal
    Poro_1VolFourierFrac = (1 - PorositeFracLocal) * Thickness * SurfFracLocal
    Poro_1VolFourierNode = 0.d0

    CellThermalSourceVol = CellThermalSourceLocal * VolCellLocal
    FracThermalSourceVol = FracThermalSourceLocal * Thickness * SurfFracLocal
    NodeThermalSourceVol = 0.d0

    ! Fourier volume
    CALL VAGFrac_SplitCellVolume( &
      NbCellLocal_Ncpus(commRank+1), &
      CellRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d" .AND. IdNodeLocal(:)%Frac == "n", &
      omegaFourierCell, &
      NodebyCellLocal, &
      PoroVolFourierCell, &
      PoroVolFourierNode)

    CALL VAGFrac_SplitCellVolume( &
      NbFracLocal_Ncpus(commRank+1), &
      FracRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d", &
      omegaFourierFrac, &
      NodebyFractureLocal, &
      PoroVolFourierFrac, &
      PoroVolFourierNode)

    CALL VAGFrac_SplitCellVolume( &
      NbCellLocal_Ncpus(commRank+1), &
      CellRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d" .AND. IdNodeLocal(:)%Frac == "n", &
      omegaFourierCell, &
      NodebyCellLocal, &
      Poro_1VolFourierCell, &
      Poro_1VolFourierNode)

    CALL VAGFrac_SplitCellVolume( &
      NbFracLocal_Ncpus(commRank+1), &
      FracRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d", &
      omegaFourierFrac, &
      NodebyFractureLocal, &
      Poro_1VolFourierFrac, &
      Poro_1VolFourierNode)

    CALL VAGFrac_SplitCellVolume( &
      NbCellLocal_Ncpus(commRank+1), &
      CellRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d" .AND. IdNodeLocal(:)%Frac == "n", &
      omegaFourierCell, &
      NodebyCellLocal, &
      CellThermalSourceVol, &
      NodeThermalSourceVol)

    CALL VAGFrac_SplitCellVolume( &
      NbFracLocal_Ncpus(commRank+1), &
      FracRocktypeLocal(2,:), &
      NbNodeLocal_Ncpus(commRank+1), &
      NodeRocktypeLocal(2,:), &
      IdNodeLocal(:)%T /= "d", &
      omegaFourierFrac, &
      NodebyFractureLocal, &
      FracThermalSourceVol, &
      NodeThermalSourceVol)

    ! check if vol is positive
    do k=1, NbCellLocal_Ncpus(commRank+1)
      if(PoroVolFourierCell(k)<eps) then
        print*, "vol fourier cell <= 0: ", k, PoroVolFourierCell(k)

        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

    do k=1, NbFracLocal_Ncpus(commRank+1)
      if(PoroVolFourierFrac(k)<eps) then

        print*, "vol fourier frac <= 0: ", k, PoroVolFourierFrac(k)
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

    ! check if vol is positive
    do k=1, NbCellLocal_Ncpus(commRank+1)
      if(Poro_1volFourierCell(k)<eps) then
        print*, "vol poro 1-phi fourier cell <= 0: ", k, Poro_1volFourierCell(k)

        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

    do k=1, NbFracLocal_Ncpus(commRank+1)
      if(Poro_1volFourierFrac(k)<eps) then

        print*, "vol poro 1-phi fourier frac <= 0: ", k, Poro_1volFourierFrac(k)
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
    end do

  end subroutine VAGFrac_VolsFourier

#endif

  ! Free
  subroutine VAGFrac_Free

    integer :: k

    ! darcy trans
    if( allocated(TkLocal_Darcy)) then
      do k=1, NbCellLocal_Ncpus(commRank+1)
        deallocate(TkLocal_Darcy(k)%pt)
      end do

      deallocate(TkLocal_Darcy)
    end if

    if( allocated(TkFracLocal_Darcy)) then
      do k=1, NbFracLocal_Ncpus(commRank+1)
        deallocate(TkFracLocal_Darcy(k)%pt)
      end do

      deallocate(TkFracLocal_Darcy)
    end if

#ifdef _THERMIQUE_

    ! fourier trans
    if( allocated(TkLocal_Fourier)) then
      do k=1, NbCellLocal_Ncpus(commRank+1)
        deallocate(TkLocal_Fourier(k)%pt)
      end do

      deallocate(TkLocal_Fourier)
    end if

    if( allocated(TkFracLocal_Fourier)) then
      do k=1, NbFracLocal_Ncpus(commRank+1)
        deallocate(TkFracLocal_Fourier(k)%pt)
      end do

      deallocate(TkFracLocal_Fourier)
    end if
#endif

    ! darcy porous vol
    deallocate(PoroVolDarcyCell)
    deallocate(PoroVolDarcyFrac)
    deallocate(PoroVolDarcyNode)

    ! darcy vol
    deallocate(VolDarcyCell)
    deallocate(VolDarcyFrac)
    deallocate(VolDarcyNode)

    ! fourier porous vol
#ifdef _THERMIQUE_

    deallocate(PoroVolFourierCell)
    deallocate(PoroVolFourierFrac)
    deallocate(PoroVolFourierNode)
    deallocate(Poro_1volFourierCell)
    deallocate(Poro_1volFourierFrac)
    deallocate(Poro_1volFourierNode)
    deallocate(CellThermalSourceVol)
    deallocate(FracThermalSourceVol)
    deallocate(NodeThermalSourceVol)
#endif

  end subroutine VAGFrac_Free


  !! Inverse of A
  subroutine VAGFrac_InvA(A, Ainv)

    double precision, intent(in) :: A(3,3)
    double precision, intent(out) :: Ainv(3,3)
    double precision :: coA(3,3)

    double precision :: det

    det = A(1,1)*( A(3,3)*A(2,2)-A(3,2)*A(2,3) )
    det = det - A(1,2)*( A(3,3)*A(2,1)-A(3,1)*A(2,3) ) 
    det = det + A(1,3)*( A(3,2)*A(2,1)-A(3,1)*A(2,2) ) 

    coA(1,1) = + ( A(3,3)*A(2,2)-A(3,2)*A(2,3) )/det 
    coA(1,2) = - ( A(3,3)*A(2,1)-A(3,1)*A(2,3) )/det 
    coA(1,3) = + ( A(3,2)*A(2,1)-A(3,1)*A(2,2) )/det  

    coA(2,1) = - ( A(3,3)*A(1,2)-A(3,2)*A(1,3) )/det 
    coA(2,2) = + ( A(3,3)*A(1,1)-A(3,1)*A(1,3) )/det 
    coA(2,3) = - ( A(3,2)*A(1,1)-A(3,1)*A(1,2) )/det  

    coA(3,1) = + ( A(2,3)*A(1,2)-A(2,2)*A(1,3) )/det 
    coA(3,2) = - ( A(2,3)*A(1,1)-A(2,1)*A(1,3) )/det 
    coA(3,3) = + ( A(2,2)*A(1,1)-A(2,1)*A(1,2) )/det  


    Ainv(1,1) = coA(1,1)
    Ainv(2,2) = coA(2,2)
    Ainv(3,3) = coA(3,3)

    Ainv(1,2) = coA(2,1)
    Ainv(2,1) = coA(1,2)

    Ainv(1,3) = coA(3,1)
    Ainv(3,1) = coA(1,3)

    Ainv(3,2) = coA(2,3)
    Ainv(2,3) = coA(3,2)

  end subroutine VAGFrac_InvA


  ! vol tetra
  subroutine VAGFrac_VolTetra(x1,x2,x3,x4,vol)

    ! volume du tetra defini par ses 4 points
    double precision, dimension(3), intent(in) :: &
        x1, x2, x3, x4

    double precision, intent(out) :: vol

    double precision, dimension(3) :: v
    double precision :: s, surf

    ! normale a 123 orientee vers 4 et surface 123
    call VAGFrac_VecNormalT(x1,x2,x3,x4,v,surf)

    ! hauteur 
    s = v(1)*(x1(1)-x4(1))+v(2)*(x1(2)-x4(2))+v(3)*(x1(3)-x4(3))
    s = dabs(s)

    vol = s*surf/3.d0

  end subroutine VAGFrac_VolTetra


  ! Vec normal 
  subroutine VAGFrac_VecNormalT(x1,x2,x3,x,v,surf)

    ! calcul du vecteur normal unitaire d'un triangle defini par les 
    ! coordonnees de ses trois points sortant par rapport 
    ! au point x,y,z > vx,vy,vz 
    ! + surface du triangle = surf 

    double precision, dimension(3), intent(in) :: x1 ,x2, x3, x
    double precision, dimension(3), intent(out) :: v
    double precision, intent(out) :: surf
    double precision, dimension(3) :: xt
    double precision :: s

    xt(:) = (x1(:)+x2(:)+x3(:))/3.d0

    v(1) = (x1(3)-x3(3))*(x1(2)-x2(2)) - (x1(2)-x3(2))*(x1(3)-x2(3))
    v(2) = (x1(1)-x3(1))*(x1(3)-x2(3)) - (x1(3)-x3(3))*(x1(1)-x2(1))
    v(3) = (x1(2)-x3(2))*(x1(1)-x2(1)) - (x1(1)-x3(1))*(x1(2)-x2(2))

    s = dsqrt(v(1)**2+v(2)**2+v(3)**2)

    surf = s/2.d0 
    v(:) = v(:)/s

    s = (xt(1)-x(1))*v(1) + (xt(2)-x(2))*v(2) + (xt(3)-x(3))*v(3)

    if (s .lt. 0.d0) then
      v(:) = - v(:)
    endif

  end subroutine VAGFrac_VecNormalT


  ! Idea: via num (local) des nodes
  subroutine VAGFrac_UnkFaceToUnkCell(k,i)

    ! k is cell
    ! i is num (local) of face    
    integer, intent(in) :: &
        k, i

    integer :: &
        in, idin, numin, &  ! num (face)
        idj, numj, &        ! num (local)
        j                   ! num (cell)

    integer :: nbNodeCell, nbNodeFace, nbFracCell

    nbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
    nbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)
    nbNodeFace = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)

    UnkFaceToUnkCell(:) = 0  ! size is nbNodeFace+nbFracFace

    ! Unk = (node,frac)
    ! one frac in a face frac
    do in = 1, nbNodeFace

      idin = NodebyFaceLocal%Pt(i) + in ! id of in
      numin = NodebyFaceLocal%Num(idin) ! num of in

      ! look for numin in NodebyCellLocal
      do j = 1, nbNodeCell

        idj = NodebyCellLocal%Pt(k) + j
        numj =  NodebyCellLocal%Num(idj)

        if(numj == numin) then ! idea is here, nums (local) are same
          UnkFaceToUnkCell(in) = j
          exit
        end if

      end do
    end do

    ! frac in Unk, one frac or zero frac in a face
    if( IdFaceLocal(i) == -2) then ! if i is frac
      do j = 1, nbFracCell

        idj = FracbyCellLocal%Pt(k) + j
        numj = FracbyCellLocal%Num(idj)

        if (numj == i) then
          UnkFaceToUnkCell(nbNodeFace+1) = j + nbNodeCell
          exit
        end if
      end do
    else ! if i is not frac
      UnkFaceToUnkCell(nbNodeFace+1) = 0
    end if

  end subroutine VAGFrac_UnkFaceToUnkCell


  ! *** Tests subroutines *** !

  ! Test Gradient: GkT
  subroutine VAGFrac_TestGkT(i, xk, xf, GkT)

    integer, intent(in) :: i
    double precision, intent(in) :: xk(3), xf(3) ! center of cell, face
    double precision, dimension(:,:), intent(in) :: GkT

    double precision :: &
        uk, uj, &
        sx, sy, sz, &
        xj(3)

    integer :: j, nj, nb

    ! test du gradient GkT 
    uk = 3.d0*xk(1) + 2.d0*xk(2) + xk(3) + 1.d0

    sx = 0.d0
    sy = 0.d0
    sz = 0.d0

    ! node
    do j=NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1) 

      nj = NodebyFaceLocal%Num(j)

      xj(:) = XNodeLocal(:,nj)
      ! print*, xj

      uj = 3.d0*xj(1) + 2.d0*xj(2) + xj(3) + 1.d0

      sx = sx + GkT(1, j-NodebyFaceLocal%Pt(i)) * (uj-uk)
      sy = sy + GkT(2, j-NodebyFaceLocal%Pt(i)) * (uj-uk)
      sz = sz + GkT(3, j-NodebyFaceLocal%Pt(i)) * (uj-uk)
    enddo

    ! frac
    xj(:) = xf(:)
    uj = 3.d0*xj(1) + 2.d0*xj(2) + xj(3) + 1.d0

    nb = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i) ! nb node of face i
    sx = sx + GkT(1, nb+1) * (uj-uk)
    sy = sy + GkT(2, nb+1) * (uj-uk)
    sz = sz + GkT(3, nb+1) * (uj-uk)

    if (dabs(sx+3.0)+ dabs(sy+2.0)+ dabs(sz+1.0) > 1.d-8) then
      print*, "Gradient error: "
      print*, "   Center of Cell: ", xk
      print*, "   Center of Face: ", xf
    end if

  end subroutine VAGFrac_TestGkT


  ! Test Gradient tangentiel: GfT
  subroutine VAGFrac_TestGfT(i, x1, x2, xf, v, GfT)

    integer, intent(in) :: i
    double precision, dimension(3), intent(in) :: &
        x1, x2, xf, v

    double precision, dimension(:,:) :: &
        GfT

    double precision :: &
        gx, gy, gz, &
        u1, u2, uf, ui, ps, &
        erx, ery, erz, &
        xi(3)

    integer :: is

    gx = 0.d0
    gy = 0.d0
    gz = 0.d0

    u1 = 3.d0*x1(1) + 2.d0*x1(2) + x1(3) - 1.d0
    u2 = 3.d0*x2(1) + 2.d0*x2(2) + x2(3) - 1.d0
    uf = 3.d0*xf(1) + 2.d0*xf(2) + xf(3) - 1.d0

    do is=1, NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)

      xi(:) = XNodeLocal(:, NodebyFaceLocal%Num( NodebyFaceLocal%Pt(i)+is) )
      ui = 3.d0*xi(1) + 2.d0*xi(2) + xi(3) - 1.d0

      gx = gx + GfT(1,is)*(uf-ui)
      gy = gy + GfT(2,is)*(uf-ui)
      gz = gz + GfT(3,is)*(uf-ui)
    enddo

    ps = 3.d0*v(1)+2.d0*v(2)+v(3) 

    erx = gx - (3.d0 - ps*v(1))
    ery = gy - (2.d0 - ps*v(2))
    erz = gz - (1.d0 - ps*v(3))

    if (dabs(erx)+ dabs(ery)+ dabs(erz) > 1.d-10) then
      print*, "Gradient Frac error: "
      print*, "   Center of Frac: ", xf
      stop
    end if

  end subroutine VAGFrac_TestGfT

end module VAGFrac
