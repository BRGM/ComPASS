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
   IdFace = 0

   do i = 1, NbFace

      if (IdFace(i) == -2) then
         NbFrac = NbFrac + 1
      else
         IdFace(i) = 0
      end if
   end do

end subroutine GlobalMesh_SetFrac

SUBROUTINE GlobalMesh_SetFaceFlags
   INTEGER :: k, m
   DOUBLE PRECISION :: xk(3)

   DO k = 1, NbFace
      ! center of face
      xk(:) = 0.d0
      DO m = NodebyFace%Pt(k) + 1, NodebyFace%Pt(k + 1)
         xk(:) = xk(:) + XNode(:, NodebyFace%Num(m))
      ENDDO
      xk(:) = xk(:)/dble(NodebyFace%Pt(k + 1) - NodebyFace%Pt(k))

      IF (ABS(xk(3) - Mesh_zmin) < eps) THEN
         FaceFlags(k) = 1
      ELSEIF (ABS(xk(3) - Mesh_zmax) < eps) THEN
         FaceFlags(k) = 2
      ELSE
         FaceFlags(k) = 0
      ENDIF
   ENDDO
END SUBROUTINE GlobalMesh_SetFaceFlags

SUBROUTINE GlobalMesh_SetFracRocktype

   FracRocktype(1, :) = 1

#ifdef _THERMIQUE_
   FracRocktype(2, :) = 1
#endif
END SUBROUTINE GlobalMesh_SetFracRocktype

SUBROUTINE GlobalMesh_SetCellFlags
   INTEGER :: k, m, i
   DOUBLE PRECISION :: xk(3)
   INTEGER :: NbThermalSource
   DOUBLE PRECISION :: ThermalSourceRadius(3)
   DOUBLE PRECISION :: ThermalSourceOffset(3)
   DOUBLE PRECISION :: ThermalSourceOffset0(3)
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: ThermalSource

   NbThermalSource = 10 ! nb of points
   ThermalSourceRadius = (/.5d0, 1.d0, 10.d0/) ! delta_x, delta_y, delta_z
   ThermalSourceOffset = (/40.d0, 0.d0, 0.d0/)
   ThermalSourceOffset0 = (/20.5d0, 0.5d0, 15.d0/)

   ALLOCATE (ThermalSource(3, NbThermalSource))

   DO i = 1, NbThermalSource
      ThermalSource(:, i) = ThermalSourceOffset0 + (i - 1)*ThermalSourceOffset
   ENDDO

   DO k = 1, NbCell
      xk(:) = 0.d0
      DO m = NodebyCell%Pt(k) + 1, NodebyCell%Pt(k + 1)
         xk(:) = xk(:) + XNode(:, NodebyCell%Num(m))
      ENDDO
      xk(:) = xk(:)/dble(NodebyCell%Pt(k + 1) - NodebyCell%Pt(k))

!    IF(MOD(INT(xk(1)),2) == MOD(INT(xk(3)),2))THEN
      IF (xk(3) <= 1.5d0) THEN
         CellFlags(k) = 2
      ELSE
         CellFlags(k) = 1
      ENDIF

      DO i = 1, NbThermalSource
         IF (ALL(ABS(xk - ThermalSource(:, i)) < ThermalSourceRadius)) THEN
            CellFlags(k) = CellFlags(k) + 2
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE (ThermalSource)
END SUBROUTINE GlobalMesh_SetCellFlags

SUBROUTINE GlobalMesh_SetCellRocktype

   CellRocktype(1, :) = MOD(CellFlags - 1, 2) + 1

#ifdef _THERMIQUE_
   CellRocktype(2, :) = MOD(CellFlags - 1, 2) + 1
#endif
!  CellTRocktype = CellFlags
END SUBROUTINE GlobalMesh_SetCellRocktype

#ifdef _THERMIQUE_
SUBROUTINE GlobalMesh_SetCellThermalSourceType

   CellThermalSourceType = MERGE(1, 0, CellFlags > 2)
END SUBROUTINE GlobalMesh_SetCellThermalSourceType

SUBROUTINE GlobalMesh_SetFracThermalSourceType

   FracThermalSourceType = 0
END SUBROUTINE GlobalMesh_SetFracThermalSourceType
#endif

SUBROUTINE GlobalMesh_SetNodeFlags
   INTEGER :: i, k, kpt
   INTEGER :: f
   DOUBLE PRECISION :: xk(3)

   DO i = 1, NbNode

      kpt = FacebyNode%Pt(i) + 1
      k = FacebyNode%Num(kpt)
      f = FaceFlags(k)
      DO kpt = FacebyNode%Pt(i) + 2, FacebyNode%Pt(i + 1)
         k = FacebyNode%Num(kpt)

         f = MAX(f, FaceFlags(k))
      ENDDO
      NodeFlags(i) = f
   ENDDO
END SUBROUTINE GlobalMesh_SetNodeFlags

! Set dirichlet boundary
! cartesian mesh:manually
! other mesh:manually or
! according to IdNodeFromfile could be served
! The node index read from mesh file is:IdNodeFromfile(integer)
subroutine GlobalMesh_SetDirBC

   integer :: i, j

   double precision :: xi, yi, zi

   ! IdNode
   if (.not. allocated(Idnode)) then
      allocate (IdNode(NbNode))
   end if

   NbDirNodeP = 0
#ifdef _THERMIQUE_
   NbDirNodeT = 0
#endif

   do i = 1, NbNode

      xi = XNode(1, i)
      yi = XNode(2, i)
      zi = XNode(3, i)

      IF (NodeFlags(i) /= 0) THEN
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
subroutine GlobalMesh_SetWellCar(nx, ny, nz)

   integer :: i, j, iwell, jwell
   integer, intent(in) :: nx, ny, nz

   ! one injection well
   NbWellInj = 0
   allocate (NbEdgebyWellInj(NbWellInj))
   allocate (NumNodebyEdgebyWellInj(2, nz, NbWellInj))

   NbEdgebyWellInj(:) = nz
   NumNodebyEdgebyWellInj(:, :, :) = -1

   iwell = 2*nx/3
   jwell = ny/2

   do i = 1, NbWellInj
      do j = 1, NbEdgebyWellInj(i)
         NumNodebyEdgebyWellInj(1, j, i) = j*(nx + 1)*(ny + 1) + jwell*(nx + 1) + iwell + 1
         NumNodebyEdgebyWellInj(2, j, i) = (j - 1)*(nx + 1)*(ny + 1) + jwell*(nx + 1) + iwell + 1
      end do
   end do

   ! 1 production well
   NbWellProd = 0
   allocate (NbEdgebyWellProd(NbWellProd))
   allocate (NumNodebyEdgebyWellProd(2, nz, NbWellProd))

   NbEdgebyWellProd(:) = nz
   NumNodebyEdgebyWellProd(:, :, :) = -1

   iwell = nx/3
   jwell = ny/2

   do i = 1, NbWellProd
      do j = 1, NbEdgebyWellProd(i)
         NumNodebyEdgebyWellProd(1, j, i) = j*(nx + 1)*(ny + 1) + jwell*(nx + 1) + iwell + 1
         NumNodebyEdgebyWellProd(2, j, i) = (j - 1)*(nx + 1)*(ny + 1) + jwell*(nx + 1) + iwell + 1
      end do
   end do

end subroutine GlobalMesh_SetWellCar
