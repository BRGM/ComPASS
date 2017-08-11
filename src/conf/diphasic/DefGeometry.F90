! *** Part of module GlobalMesh.F90 *** !


! IdFace is read from mesh file
! Set fracture face according to IdFace or manually
! Reset IdFace if necessary
subroutine GlobalMesh_SetFrac

  integer :: i, j
  double precision :: xi(3)

  NbFrac = 0

  do i=1, NbFace

     if( IdFace(i) == -2 ) then
        NbFrac = NbFrac + 1
     else
        IdFace(i) = 0
     end if
  end do

  ! IdFace(:) = 0
  
end subroutine GlobalMesh_SetFrac



SUBROUTINE GlobalMesh_SetFaceFlags
  INTEGER :: k, m
  DOUBLE PRECISION :: xk(3)

  NbFrac = 0
  IdFace = 0

  DO k=1,NbFace
    ! center of face
    xk(:) = 0.d0
    DO m = NodebyFace%Pt(k)+1, NodebyFace%Pt(k+1)
      xk(:) = xk(:) + XNode(:, NodebyFace%Num(m))
    ENDDO
    xk(:) = xk(:)/dble(NodebyFace%Pt(k+1) - NodebyFace%Pt(k))

    IF(xk(3) <= 1.d0)THEN
      FaceFlags(k) = 2
    ELSE
      FaceFlags(k) = 1
    ENDIF
  ENDDO
END SUBROUTINE GlobalMesh_SetFaceFlags


SUBROUTINE GlobalMesh_SetCellFlags
  INTEGER :: k, m
  DOUBLE PRECISION :: xk(3)

  DO k=1,NbCell
    xk(:) = 0.d0
    DO m = NodebyCell%Pt(k)+1, NodebyCell%Pt(k+1)
      xk(:) = xk(:) + XNode(:, NodebyCell%Num(m))
    ENDDO
    xk(:) = xk(:)/dble(NodebyCell%Pt(k+1) - NodebyCell%Pt(k))

    IF(xk(3) <= 1.d0)THEN
      CellFlags(k) = 2
    ELSE
      CellFlags(k) = 1
    ENDIF
  ENDDO
END SUBROUTINE GlobalMesh_SetCellFlags



! Set dirichlet boundary
!    cartesian mesh:   manually
!    other mesh: manually or
!                according to IdNodeFromfile could be served
! The node index read from mesh file is: IdNodeFromfile (integer)
subroutine GlobalMesh_SetDirBC

  integer :: i, j

  double precision :: xi, yi, zi

  ! IdNode
  if(.not. allocated(Idnode)) then
     allocate(IdNode(NbNode))
  end if

  NbDirNodeP = 0
#ifdef _THERMIQUE_
  NbDirNodeT = 0
#endif

  do i=1, NbNode

     xi = XNode(1,i)
     yi = XNode(2,i)
     zi = XNode(3,i)

     if( abs(zi-Mesh_zmin)<eps &
       .or. abs(zi-Mesh_zmax)<eps ) then 

        NbDirNodeP = NbDirNodeP + 1
        IdNode(i)%P = "d"

#ifdef _THERMIQUE_
        NbDirNodeT = NbDirNodeT + 1
        IdNode(i)%T = "d"
#endif
     else
        IdNode(i)%P = "i"
#ifdef _THERMIQUE_
        IdNode(i)%T = "i"
#endif
     end if

  end do
  
end subroutine GlobalMesh_SetDirBC


!> \brief User set well for car mesh
subroutine GlobalMesh_SetWellCar(nx,ny,nz)

  integer :: i, j, iwell, jwell
  integer, intent(in) :: nx, ny, nz

  ! one injection well
  NbWellInj = 0
  allocate(NbEdgebyWellInj(NbWellInj))
  allocate(NumNodebyEdgebyWellInj(2,nz,NbWellInj))
  
   NbEdgebyWellInj(:) = nz
   NumNodebyEdgebyWellInj(:,:,:) = -1    

  iwell = 2*nx/3
  jwell = ny/2

  do i=1,NbWellInj
      do j=1,NbEdgebyWellInj(i)
         NumNodebyEdgebyWellInj(1,j,i) = j*(nx+1)*(ny+1) + jwell * (nx+1) + iwell + 1
         NumNodebyEdgebyWellInj(2,j,i) = (j-1)*(nx+1)*(ny+1) + jwell * (nx+1) + iwell + 1
      end do
   end do

  ! 1 production well
  NbWellProd = 0
  allocate(NbEdgebyWellProd(NbWellProd))
  allocate(NumNodebyEdgebyWellProd(2,nz,NbWellProd))

  NbEdgebyWellProd(:) = nz
  NumNodebyEdgebyWellProd(:,:,:) = -1    

  iwell = nx/3
  jwell = ny/2
  
  do i=1,NbWellProd
     do j=1,NbEdgebyWellProd(i)
        NumNodebyEdgebyWellProd(1,j,i) = j*(nx+1)*(ny+1) + jwell * (nx+1) + iwell + 1
        NumNodebyEdgebyWellProd(2,j,i) = (j-1)*(nx+1)*(ny+1) + jwell * (nx+1) + iwell + 1
     end do
  end do
  
end subroutine GlobalMesh_SetWellCar
