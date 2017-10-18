!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

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

!     xi = XNode(1,i)
!     yi = XNode(2,i)
!     zi = XNode(3,i)
!
!     if(  abs( abs(xi)-Mesh_xmin)<eps .or. &
!          abs( abs(yi)-Mesh_ymin)<eps .or. &
!          abs( abs(xi)-Mesh_xmax)<eps .or. &
!          abs( abs(yi)-Mesh_ymax)<eps) then
!
!        NbDirNodeP = NbDirNodeP + 1
!        IdNode(i)%P = "d"
!
!#ifdef _THERMIQUE_
!        NbDirNodeT = NbDirNodeT + 1
!        IdNode(i)%T = "d"
!#endif
!     else
        IdNode(i)%P = "i"
#ifdef _THERMIQUE_
        IdNode(i)%T = "i"
#endif
     !end if

  end do

end subroutine GlobalMesh_SetDirBC


!> \brief User set well for car mesh
subroutine GlobalMesh_SetWellCar(nx,ny,nz)

  integer :: i, j, iwell, jwell
  integer, intent(in) :: nx, ny, nz

  ! one injection well
  NbWellInj = 1
  allocate(NbEdgebyWellInj(NbWellInj))
  allocate(NumNodebyEdgebyWellInj(2,nz,NbWellInj))

   NbEdgebyWellInj(1) = nz
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
  NbWellProd = 1
  allocate(NbEdgebyWellProd(NbWellProd))
  allocate(NumNodebyEdgebyWellProd(2,nz,NbWellProd))

  NbEdgebyWellProd(1) = nz
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

!> \brief User set porosity
!!
!! PorositeCell contains the porosity of each cell;
!! PorositeFace is restricted to PorositeFrac
!! then distributed to PorositeFracLocal
subroutine GlobalMesh_SetPorosite

  allocate(PorositeCell(NbCell))
  allocate(PorositeFace(NbFace))

  PorositeCell(:) = 1.d-1
  PorositeFace(:) = 4.d-1

  ! TODO: read from file if needed

end subroutine GlobalMesh_SetPorosite
