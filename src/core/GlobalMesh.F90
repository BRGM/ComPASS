!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

!> \brief Main subroutine GlobalMesh_make.   <br>
!! Contains the connectivity of the global mesh.  <br>
!! 1. Read mesh   <br>
!!      nbCell, nbFace, nbNode,   <br>
!!      NodebyFace, FacebyCell, XNode,   <br>
!!    and IdFace, Dir Part   <br>
!! 2.a. calcul CellbyCell for metis   <br>
!!    use NodebyFace, FacebyCell     <br>
!!    to -> NodebyCell -> CellbyNode -> CellbyCell   <br>
!! 2.b. calcul CellbyFace using FacebyCell    <br>
!! 3. IdFace, Mesh Part

module GlobalMesh
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   ! This for array that are interfaced with python/C++
   use iso_c_binding, only: c_int, c_double, c_int8_t
   use mpi, only: MPI_Abort
   use CommonMPI
   use CommonType
   use DefModel
   use DefWell
   use DefMSWell

#else
   ! This for array that are interfaced with python/C++
   use iso_c_binding, only: c_int, c_double, c_int8_t

   use mpi, only: MPI_Abort

   use CommonMPI, only: compass_comm_world
   use CommonType, only: Type_IdNode, CSR, &
                         CommonType_deallocCSR
   use DefModel, only: indthermique

   use DefWell, only: &
      TYPE_CSRDataNodeWell, NodeDatabyWellInj, NodeDatabyWellProd, &
      DefWell_deallocCSRDataWell, &
      GlobalWellGeometries, DefWell_clear_global_well_geometries

   use DefMSWell, only: &
      NodeDatabyMSWell
   ! use SchemeParameters ! FIXME: to be removed: for eps

#endif
   implicit none

#ifdef _DISPMODULE_
   use dispmodule
#endif

   ! Mesh
   integer :: &
      NbCell, &  !< Total number of cells
      NbNode, &  !< Total number of nodes
      NbFace, &  !< Total number of faces
      NbFrac, &  !< Total number of fracture faces
      NbDirNodeP !< Total number of Dirichlet nodes Darcy

   type(GlobalWellGeometries) :: GeomInj, GeomProd, GeomMSWell

#ifdef _THERMIQUE_
   integer :: &
      NbDirNodeT    !< Total number of Dirichlet nodes Fourier
#endif

   double precision, protected :: &
      Mesh_xmax, & !< Global size of the mesh xmax, known by proc master only
      Mesh_xmin, & !< Global size of the mesh xmin, known by proc master only
      Mesh_ymax, & !< Global size of the mesh ymax, known by proc master only
      Mesh_ymin, & !< Global size of the mesh ymin, known by proc master only
      Mesh_zmax, & !< Global size of the mesh zmax, known by proc master only
      Mesh_zmin    !< Global size of the mesh zmin, known by proc master only

   ! kind=c_double is because array is interfaced with python/C++ (cf. target attirbute)
   real(kind=c_double), allocatable, dimension(:, :), protected, target :: &
      XNode      !< Global coordinates of nodes

   integer(c_int), allocatable, dimension(:), target :: &
      NodeFlags, &
      CellFlags, &
      FaceFlags

   integer(c_int8_t), allocatable, dimension(:), target :: &
      CellTypes, &
      FaceTypes

   integer(c_int), allocatable, dimension(:, :), target :: &
      NodeRocktype, &
      CellRocktype, &
      FracRocktype

   type(CSR), protected :: &
      FacebyCell, & !< CSR list of Faces surrounding each Cell
      NodebyFace    !< CSR list of Nodes surrounding each Face

   ! Connectivities
   type(CSR), protected :: &
      NodebyCell, &  !< CSR list of Nodes surrounding each Cell
      FacebyNode, & !< CSR list of Faces surrounding each Node
      CellbyNode, & !< CSR list of Cells surrounding each Node
      CellbyCell, & !< CSR list of Cells surrounding each Cell
      CellbyFace    !< CSR list of Cells surrounding each Face

   ! Connectivities Well
   type(CSR), protected :: &
      NodebyWellInj, &   !< CSR list of Nodes of each injection Well
      NodebyWellProd, &  !< CSR list of Nodes of each production Well
      NodebyMSWell       !< CSR list of Nodes of each MSWell

   integer, dimension(:), allocatable :: &
      MSWellNodebyNode     !<  mswell global node numbering

   !type(FractureInfoCOC), protected :: &
   !type(COC), protected :: &
   type(CSR), protected :: &
      FracbyNode

   !integer(c_int), allocatable, dimension(:), target :: &
   !  FacetoFracGlobal, FractoFaceGlobal

   ! IdFace is filled in GlobalMesh_ReadMesh for Dir
   !                     GlobalMesh_Frac for Frac: IdFace(fracface)=-2
   !> \todo FIXME: protected attribute has been removed
   integer(c_int), allocatable, dimension(:), target :: &
      IdFace !< Identifier of each Face (1 to 6 to identify boundary faces and -2 for fracture faces)

   ! IdCell
   integer, allocatable, dimension(:), protected :: &
      IdCell !< Integer Identifier of each Cell, useful if heterogeneous domain (different rocktype)

   ! IdNode
   ! FIXME: protected attribute has been removed
   type(Type_IdNode), allocatable, dimension(:), target :: &
      IdNode !< Characters Identifiers of each node (related to node own or ghost; boundaries and fractures)

   ! Porosite
   ! FIXME: protected has been removed to access arrays from C
   real(c_double), allocatable, dimension(:), target :: &
      PorositeCell, & !< Porosity of each Cell
      PorositeFrac    !< Porosity of each fracture face

   ! Permeability
   !> \todo FIXME: protected has been removed to access array from C
   real(c_double), allocatable, dimension(:, :, :), target :: &
      PermCell !< Permeability tensor for each cell
   ! FIXME: protected has been removed to access array from C
   real(c_double), allocatable, dimension(:), target :: &
      PermFrac !< Permeability scalr value for each fracture face

#ifdef _THERMIQUE_
   ! Thermal conductivity
   real(c_double), allocatable, dimension(:, :, :), target :: &
      CondThermalCell !< Thermal conductivity tensor for each cell
   real(c_double), allocatable, dimension(:), target :: &
      CondThermalFrac !< Thermal conductivity scalar value for each fracture face

   ! Thermal source
   !integer(c_int), allocatable, dimension(:), target :: CellThermalSourceType
   !integer(c_int), allocatable, dimension(:), target :: FracThermalSourceType

   real(c_double), allocatable, dimension(:), target :: CellThermalSource
   real(c_double), allocatable, dimension(:), target :: FracThermalSource
#endif

   ! Used to ouput well information
   integer, protected :: fdGm
   ! integer, protected :: fdGm_unit

   ! main subroutine in this module
   public :: &
      !GlobalMesh_Make_read_file, &
      !GlobalMesh_Make_post_read, &
      GlobalMesh_free

   !> \todo FIXME: All routines are made public here
   !private :: &
   public :: &
      GlobalMesh_MeshBoundingBox, & ! computes mesh bounding box (Mesh_xmin, Mesh_xmax...)
      !GlobalMesh_ReadMeshCar,          & ! generate cartesian mesh
      !GlobalMesh_ReadMeshFromFile,     & ! read mesh from file
      GlobalMesh_FaceByNodeGlobal, & ! make FacebyNode
      GlobalMesh_CellByNodeGlobal, & ! make CellbyNode
      GlobalMesh_CellByCellGlobal, & ! make CellbyCell
      GlobalMesh_CellbyFaceGlobal, & ! make CellbyFace
      GlobalMesh_NodeOfFrac, & ! IdNode()%Frac
      GlobalMesh_FracbyNode, & ! Make FracbyNode for Well Index
      GlobalMesh_WellConnectivity, & ! NodebyWell and NodeDatabyWell
      GlobalMesh_SetFrac, &
      GlobalMesh_create_mesh, &
      GlobalMesh_allocate_petrophysics, &
      GlobalMesh_set_all_rocktypes

contains

   !> \brief Determine Bounding box of the mesh
   !! Fill Mesh_xmin, Mesh_xmax and similar
   subroutine GlobalMesh_MeshBoundingBox() &
      bind(C, name="GlobalMesh_mesh_bounding_box")

      integer :: i

      Mesh_xmax = XNode(1, 1)
      Mesh_xmin = XNode(1, 1)
      Mesh_ymax = XNode(2, 1)
      Mesh_ymin = XNode(2, 1)
      Mesh_zmax = XNode(3, 1)
      Mesh_zmin = XNode(3, 1)

      do i = 2, NbNode

         if (XNode(1, i) < Mesh_xmin) then
            Mesh_xmin = XNode(1, i)
         else if (XNode(1, i) > Mesh_xmax) then
            Mesh_xmax = XNode(1, i)
         end if

         if (XNode(2, i) < Mesh_ymin) then
            Mesh_ymin = XNode(2, i)
         else if (XNode(2, i) > Mesh_ymax) then
            Mesh_ymax = XNode(2, i)
         end if

         if (XNode(3, i) < Mesh_zmin) then
            Mesh_zmin = XNode(3, i)
         else if (XNode(3, i) > Mesh_zmax) then
            Mesh_zmax = XNode(3, i)
         end if

      end do

      print *, "Bounding box:"
      print *, Mesh_xmin, "< X <", Mesh_xmax
      print *, Mesh_ymin, "< Y <", Mesh_ymax
      print *, Mesh_zmin, "< Z <", Mesh_zmax

   end subroutine GlobalMesh_MeshBoundingBox

   !> \brief Set number of frac faces ?
   subroutine GlobalMesh_SetFrac() &
      bind(C, name="GlobalMesh_set_frac")

      integer :: i

      NbFrac = 0
      do i = 1, NbFace
         if (IdFace(i) == -2) then
            NbFrac = NbFrac + 1
         else
            IdFace(i) = 0
         end if
      end do

   end subroutine GlobalMesh_SetFrac

   !> \brief Compute all global connectivities
   !! FaceByNode, CellByNode, CellbyCell, CellbyFace
   subroutine GlobalMesh_Compute_all_connectivies() &
      bind(C, name="GlobalMesh_compute_all_connectivies")

      call GlobalMesh_FaceByNodeGlobal
      call GlobalMesh_CellByNodeGlobal
      call GlobalMesh_CellbyCellGlobal
      call GlobalMesh_CellbyFaceGlobal

   end subroutine GlobalMesh_Compute_all_connectivies

!> \brief Allocate the rocktype vectors for node, frac, cell
   subroutine GlobalMesh_allocate_rocktype()

      ALLOCATE (NodeRocktype(IndThermique + 1, Nbnode))
      NodeRocktype = 1
      ALLOCATE (CellRocktype(IndThermique + 1, NbCell))
      CellRocktype = 1
      ALLOCATE (FracRocktype(IndThermique + 1, NbFace))
      FracRockType = 1

   end subroutine GlobalMesh_allocate_rocktype

!> \brief Deallocate the rocktype vectors for node, frac, cell
   subroutine GlobalMesh_deallocate_rocktypes()

      deallocate (NodeRocktype)
      deallocate (FracRocktype)
      deallocate (CellRocktype)

   end subroutine GlobalMesh_deallocate_rocktypes

   subroutine GlobalMesh_allocate_petrophysics() &
      bind(C, name="GlobalMesh_allocate_petrophysics")

      allocate (PorositeCell(NbCell))
      PorositeCell(:) = 0.d0
      allocate (PorositeFrac(NbFace))
      PorositeFrac(:) = 0.d0

      allocate (PermCell(3, 3, NbCell))
      PermCell(:, :, :) = 0.d0

      allocate (PermFrac(NbFace))
      PermFrac(:) = 0.d0

#ifdef _THERMIQUE_

      allocate (CondThermalCell(3, 3, NbCell))
      CondThermalCell(:, :, :) = 0.d0

      allocate (CondThermalFrac(NbFace))
      CondThermalFrac(:) = 0.d0

      allocate (CellThermalSource(NbCell))
      CellThermalSource = 0.d0
      allocate (FracThermalSource(NbFace))
      FracThermalSource = 0.d0

#endif

!call CommonType_deallocCSR(FractureIdbyNode)

   end subroutine GlobalMesh_allocate_petrophysics

!> \brief Set rocktype for a type of cv
!! (node or frac or cell)
   SUBROUTINE GlobalMesh_SetRocktype( &
      NbNode, &
      NbElem, &
      IsRocktypeNode, &
      ElemRocktype, &
      ElemPermeability, &
      ElembyNode, &
      Rocktype)

      INTEGER, INTENT(IN) :: NbNode
      INTEGER, INTENT(IN) :: NbElem
      LOGICAL, INTENT(IN) :: IsRocktypeNode(NbNode)
      INTEGER, INTENT(IN) :: ElemRocktype(NbElem)
      DOUBLE PRECISION, INTENT(IN) :: ElemPermeability(NbElem)
      TYPE(CSR), INTENT(IN) :: ElembyNode

      INTEGER, INTENT(INOUT) :: Rocktype(NbNode)

      INTEGER :: i
      INTEGER :: kpt, k
      INTEGER :: rt
      DOUBLE PRECISION :: v, vk

      DO i = 1, NbNode
         IF (IsRocktypeNode(i)) THEN
            kpt = ElembyNode%Pt(i) + 1
            k = ElembyNode%Num(kpt)

            rt = ElemRocktype(k)
            v = ElemPermeability(k)
            DO kpt = ElembyNode%Pt(i) + 2, ElembyNode%Pt(i + 1)
               k = ElembyNode%Num(kpt)

               vk = ElemPermeability(k)
               IF (rt /= ElemRocktype(k) .AND. vk > v) THEN
                  rt = ElemRocktype(k)
                  v = vk
               ENDIF
            ENDDO

            Rocktype(i) = rt
         ENDIF
      ENDDO
   END SUBROUTINE GlobalMesh_SetRocktype

!> \brief Set all rocktypes
!! for node, cell, frac
!! using the permeability and the thermal conductivity
   subroutine GlobalMesh_set_all_rocktypes() &
      bind(C, name="GlobalMesh_set_all_rocktypes")

      CALL GlobalMesh_SetRocktype( &
         NbNode, NbCell, &
         IdNode%Frac /= "y", &
         CellRocktype(1, :), &
         MAXVAL(RESHAPE(PermCell, (/9, NbCell/)), 1), &
         CellbyNode, &
         NodeRocktype(1, :))

      CALL GlobalMesh_SetRocktype( &
         NbNode, NbFace, &
         IdNode%Frac == "y", &
         FracRocktype(1, :), &
         PermFrac, &
         !FractureIdbyNode, &
         FracbyNode, &
         NodeRocktype(1, :))

#ifdef _THERMIQUE_

      CALL GlobalMesh_SetRocktype( &
         NbNode, NbCell, &
         IdNode%Frac /= "y", &
         CellRocktype(2, :), &
         MAXVAL(RESHAPE(CondThermalCell, (/9, NbCell/)), 1), &
         CellbyNode, &
         NodeRocktype(2, :))

      CALL GlobalMesh_SetRocktype( &
         NbNode, &
         NbFace, &
         IdNode%Frac == "y", &
         FracRocktype(2, :), &
         CondThermalFrac, &
         !FractureIdbyNode, &
         FracbyNode, &
         NodeRocktype(2, :))

#endif

   end subroutine GlobalMesh_set_all_rocktypes

   subroutine GlobalMesh_allocate_flags()

      call GlobalMesh_deallocate_flags

      allocate (NodeFlags(NbNode))
      allocate (CellFlags(NbCell))
      allocate (FaceFlags(NbFace))
      allocate (CellTypes(NbCell))
      allocate (FaceTypes(NbFace))

      ! Init all flags
      NodeFlags = 0
      CellFlags = 0
      FaceFlags = 0

   end subroutine GlobalMesh_allocate_flags

   subroutine GlobalMesh_deallocate_flags()

      if (allocated(NodeFlags)) deallocate (NodeFlags)
      if (allocated(CellFlags)) deallocate (CellFlags)
      if (allocated(FaceFlags)) deallocate (FaceFlags)
      if (allocated(CellTypes)) deallocate (CellTypes)
      if (allocated(FaceTypes)) deallocate (FaceTypes)

   end subroutine GlobalMesh_deallocate_flags

   !> \brief Fill global connectivity CellbyCell with the numero of
   !! the cells surrounding one cell.
   !!
   !! The vector is used to build the
   !!   Partition mesh with Metis.    <br>
   !! Output:
   !!   CellbyCellGlobal   <br>
   !! Input:
   !!   CellbyNode, NodebyCell
   subroutine GlobalMesh_CellByCellGlobal

      integer :: i, j, k
      integer, allocatable, dimension(:) ::  colorCell
      integer :: nbtempNode, nbtempCell, counterNumCellbyCell, &
                 beginNode, loadCell, loadNode, beginCell

      counterNumCellbyCell = 0
      allocate (colorCell(NbCell))
      colorCell(:) = 0

      CellbyCell%Nb = NbCell
      allocate (CellbyCell%Pt(NbCell + 1))

      ! 1st step - counting
      do i = 1, NbCell
         nbtempNode = NodebyCell%Pt(i + 1) - NodebyCell%Pt(i)
         beginNode = NodebyCell%Pt(i) + 1
         ! Loop over the Nodes in Cell i
         do j = 1, nbtempNode
            loadNode = NodebyCell%Num(beginNode + j - 1)
            nbtempCell = CellbyNode%Pt(loadNode + 1) - CellbyNode%Pt(loadNode)
            beginCell = CellbyNode%Pt(loadNode) + 1
            ! Loop over the Cells surrounding node j
            do k = 1, nbtempCell
               loadCell = CellbyNode%Num(beginCell + k - 1)
               if ((colorCell(loadCell) == i) .or. (loadCell == i)) then
                  ! Cell k has already been counted, nothing is done
               else
                  ! Cell k is counting
                  colorCell(loadCell) = i
                  counterNumCellbyCell = counterNumCellbyCell + 1
               endif
            enddo
         enddo
      enddo

      allocate (CellbyCell%Num(counterNumCellbyCell))

      ! 2nd step - filling
      counterNumCellbyCell = 0
      colorCell(:) = 0

      CellbyCell%Pt(1) = 0
      do i = 1, NbCell
         nbtempNode = NodebyCell%Pt(i + 1) - NodebyCell%Pt(i)
         beginNode = NodebyCell%Pt(i) + 1
         ! Loop over the Nodes in Cell i
         do j = 1, nbtempNode
            loadNode = NodebyCell%Num(beginNode + j - 1)
            nbtempCell = CellbyNode%Pt(loadNode + 1) - CellbyNode%Pt(loadNode)
            beginCell = CellbyNode%Pt(loadNode) + 1
            ! Loop over the Cells surrounding node j
            do k = 1, nbtempCell
               loadCell = CellbyNode%Num(beginCell + k - 1)
               if ((colorCell(loadCell) == i) .or. (loadCell == i)) then
                  ! Cell k has already been stored, nothing is done
               else
                  ! Cell k is stored
                  colorCell(loadCell) = i
                  counterNumCellbyCell = counterNumCellbyCell + 1
                  CellbyCell%Num(counterNumCellbyCell) = loadCell
               endif
            enddo
         enddo
         CellbyCell%Pt(i + 1) = counterNumCellbyCell
      enddo

      deallocate (colorCell)

   end subroutine GlobalMesh_CellByCellGlobal

   !> \brief This subroutine is not called.
   !! Build NodebyCell using NodebyFace and FacebyCell.
   !! Output:
   !!   NodebyCell   <br>
   !! Use:
   !!   NodebyFace, FacebyCell
   subroutine GlobalMesh_NodeByCellGlobal

      integer :: i, j, k
      integer :: counterNumNodebyCell = 0
      integer, allocatable, dimension(:) :: colorNodes
      integer :: beginFace, nbFacetempCell
      integer :: beginNode, nbNodetempFace, faceLoad

      allocate (colorNodes(NbNode))
      colorNodes(:) = 0

      ! 1st step - counting
      ! Loop over every Cell
      do i = 1, NbCell
         nbFacetempCell = FacebyCell%Pt(i + 1) - FacebyCell%Pt(i)
         beginFace = FacebyCell%Pt(i) + 1

         ! Loop over the Faces of Cell i
         do j = 1, nbFacetempCell
            faceLoad = FacebyCell%Num(beginFace + j - 1)
            nbNodetempFace = NodebyFace%Pt(faceLoad + 1) - NodebyFace%Pt(faceLoad)
            beginNode = NodebyFace%Pt(faceLoad) + 1

            ! Loop over the Nodes of Face j
            do k = 1, nbNodetempFace
               if (colorNodes(NodebyFace%Num(beginNode + k - 1)) == i) then
                  ! Node k has already been counted, nothing is done
               else
                  ! Node k is counted
                  colorNodes(NodebyFace%Num(beginNode + k - 1)) = i
                  counterNumNodebyCell = counterNumNodebyCell + 1
               endif
            enddo
         enddo
      enddo

      NodebyCell%Nb = NbCell
      allocate (NodebyCell%Pt(NbCell + 1), NodebyCell%Num(counterNumNodebyCell))
      NodebyCell%Pt(1) = 0

      ! 2nd step - filling
      colorNodes(:) = 0
      counterNumNodebyCell = 0

      ! Loop over every Cell
      do i = 1, NbCell
         nbFacetempCell = FacebyCell%Pt(i + 1) - FacebyCell%Pt(i)
         beginFace = FacebyCell%Pt(i) + 1

         ! Loop over the Faces of Cell i
         do j = 1, nbFacetempCell
            faceLoad = FacebyCell%Num(beginFace + j - 1)
            nbNodetempFace = NodebyFace%Pt(faceLoad + 1) - NodebyFace%Pt(faceLoad)
            beginNode = NodebyFace%Pt(faceLoad) + 1

            ! Loop over the Nodes of Face j
            do k = 1, nbNodetempFace
               if (colorNodes(NodebyFace%Num(beginNode + k - 1)) == i) then
                  ! Node k has already been stored, nothing is done
               else
                  ! Node k is stored
                  counterNumNodebyCell = counterNumNodebyCell + 1
                  colorNodes(NodebyFace%Num(beginNode + k - 1)) = i
                  NodebyCell%Num(counterNumNodebyCell) = NodebyFace%Num(beginNode + k - 1)
               endif
            enddo
         enddo
         NodebyCell%Pt(i + 1) = counterNumNodebyCell
      enddo

      deallocate (colorNodes)

   end subroutine GlobalMesh_NodeByCellGlobal

   !> \brief Make CellbyNode using NodebyCell.
   !! Output:
   !!  CellbyNode   <br>
   !! Use:
   !!   NodebyCell
   subroutine GlobalMesh_CellByNodeGlobal

      integer :: i, j
      integer, allocatable, dimension(:)  ::  nbCellbyNode
      integer :: counterNumCellbyNode, nbtempCell, beginNode, loadNode

      allocate (nbCellbyNode(NbNode))
      nbCellbyNode(:) = 0
      counterNumCellbyNode = 0

      CellbyNode%Nb = NbNode
      allocate (CellbyNode%Pt(NbNode + 1))

      ! 1st step - counting
      do i = 1, NbCell
         nbtempCell = NodebyCell%Pt(i + 1) - NodebyCell%Pt(i)
         beginNode = NodebyCell%Pt(i) + 1
         ! Loop over the nodes of cell i
         do j = 1, nbtempCell
            loadNode = NodebyCell%Num(beginNode + j - 1)
            ! Number of cells surrounding node j
            nbCellbyNode(loadNode) = nbCellbyNode(loadNode) + 1
         enddo
      enddo

      CellbyNode%Pt(:) = 0
      do i = 1, NbNode
         counterNumCellbyNode = counterNumCellbyNode + nbCellbyNode(i)
         CellbyNode%Pt(i + 1) = CellbyNode%Pt(i) + nbCellbyNode(i)
      enddo

      allocate (CellbyNode%Num(counterNumCellbyNode))

      ! 2nd step - filling
      nbCellbyNode(:) = 0
      do i = 1, NbCell
         nbtempCell = NodebyCell%Pt(i + 1) - NodebyCell%Pt(i)
         beginNode = NodebyCell%Pt(i) + 1

         ! Loop over the nodes of cell i
         do j = 1, nbtempCell
            loadNode = NodebyCell%Num(beginNode + j - 1)
            nbCellbyNode(loadNode) = nbCellbyNode(loadNode) + 1
            ! Number of cells surrounding node j
            CellbyNode%Num(nbCellbyNode(loadNode) + CellbyNode%Pt(loadNode)) = i
         enddo

      enddo

      deallocate (nbCellbyNode)

   end subroutine GlobalMesh_CellByNodeGlobal

   !> \brief Make CellbyFace using FacebyCell.
   !! Output:
   !!  CellbyFace   <br>
   !! Use:
   !!  FacebyCell
   subroutine GlobalMesh_CellByFaceGlobal

      integer :: i, k, nuf
      integer, allocatable, dimension(:) ::  cptMaille

      allocate (cptMaille(NbFace))
      cptMaille(:) = 0

      CellbyFace%Nb = NbFace
      allocate (CellbyFace%Pt(NbFace + 1))

      ! 1st step - counting
      do k = 1, NbCell
         do i = FacebyCell%Pt(k) + 1, FacebyCell%Pt(k + 1)
            nuf = FacebyCell%Num(i)
            cptMaille(nuf) = cptMaille(nuf) + 1
         enddo
      enddo

      CellbyFace%Pt(1) = 0
      do i = 1, NbFace
         CellbyFace%Pt(i + 1) = CellbyFace%Pt(i) + cptMaille(i)
      enddo

      allocate (CellbyFace%Num(CellbyFace%Pt(NbFace + 1)))

      cptMaille(:) = 0
      ! 2nd step - filling
      do k = 1, NbCell
         do i = FacebyCell%Pt(k) + 1, FacebyCell%Pt(k + 1)
            nuf = FacebyCell%Num(i)
            cptMaille(nuf) = cptMaille(nuf) + 1
            CellbyFace%Num(CellbyFace%Pt(nuf) + cptMaille(nuf)) = k
         enddo
      enddo

      deallocate (cptMaille)

   end subroutine GlobalMesh_CellByFaceGlobal

   !> \brief Deallocate vectors of GlobalMesh.
   subroutine GlobalMesh_free

      deallocate (XNode)

      call GlobalMesh_deallocate_flags
      call GlobalMesh_deallocate_rocktypes

      call CommonType_deallocCSR(FacebyCell)
      call CommonType_deallocCSR(NodebyFace)
      call CommonType_deallocCSR(NodebyCell)
      call CommonType_deallocCSR(FacebyNode)
      call CommonType_deallocCSR(CellbyNode)
      !call CommonType_deallocFracInfoCOC(FracbyNode)
      !call CommonType_deallocCOC(FracbyNode)
      call CommonType_deallocCSR(FracbyNode)
      call CommonType_deallocCSR(CellbyCell)
      call CommonType_deallocCSR(CellbyFace)
      call CommonType_deallocCSR(NodebyWellInj)
      call CommonType_deallocCSR(NodebyWellProd)
      call CommonType_deallocCSR(NodebyMSWell)
      call DefWell_deallocCSRDataWell(NodeDatabyWellInj)
      call DefWell_deallocCSRDataWell(NodeDatabyWellProd)
      call DefWell_deallocCSRDataWell(NodeDatabyMSWell)

      deallocate (IdCell)
      deallocate (IdFace)
      deallocate (IdNode)
      call DefWell_clear_global_well_geometries(GeomInj)
      call DefWell_clear_global_well_geometries(GeomProd)
      call DefWell_clear_global_well_geometries(GeomMSWell)
      deallocate (PermCell)
      deallocate (PermFrac)
      deallocate (MSWellNodebyNode)

#ifdef _THERMIQUE_
      deallocate (CondThermalCell)
      deallocate (CondThermalFrac)
      deallocate (CellThermalSource)
      deallocate (FracThermalSource)
#endif

   end subroutine GlobalMesh_free

   !> \brief Fill IdNode()%Frac using IdFace.
   !!
   !! IdNode()\%Frac contains a character:                 <br>
   !! 'y': yes if node is in a fracture                    <br>
   !! 'n': no if node is not in a fracture
   subroutine GlobalMesh_NodeOfFrac() &
      bind(C, name="GlobalMesh_node_of_frac")

      integer :: k, i, numi
      integer :: Ierr, errcode ! used for MPI_Abort
      if (allocated(IdNode)) then
         write (*, *) 'IdNode should not be allocated'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      else
         allocate (IdNode(NbNode))
      end if

      do k = 1, NbNode
         IdNode(k)%Proc = "o" ! not used
         IdNode(k)%Frac = "n"
      end do

      do k = 1, NbFace
         if (IdFace(k) == -2) then
            do i = NodebyFace%Pt(k) + 1, NodebyFace%Pt(k + 1)
               numi = NodebyFace%Num(i)
               IdNode(numi)%Frac = "y" ! node i is in frac k
            end do
         end if
      end do
   end subroutine GlobalMesh_NodeOfFrac

   !> \brief Make FacebyNode using NodebyFace.   <br>
   !! Output:
   !!  FacebyNode   <br>
   !! Use:
   !!   NodebyFace
   subroutine GlobalMesh_FaceByNodeGlobal

      integer :: i, j
      integer, allocatable, dimension(:)  ::  nbFacebyNode
      integer :: counterNumFacebyNode, nbtempFace, beginNode, loadNode

      allocate (nbFacebyNode(NbNode))
      nbFacebyNode(:) = 0
      counterNumFacebyNode = 0

      FacebyNode%Nb = NbNode
      allocate (FacebyNode%Pt(NbNode + 1))

      ! 1st step - counting
      do i = 1, NbFace
         nbtempFace = NodebyFace%Pt(i + 1) - NodebyFace%Pt(i)
         beginNode = NodebyFace%Pt(i) + 1
         ! Loop over the nodes of cell i
         do j = 1, nbtempFace
            loadNode = NodebyFace%Num(beginNode + j - 1)
            ! Number of cells surrounding node j
            nbFacebyNode(loadNode) = nbFacebyNode(loadNode) + 1
         enddo
      enddo

      FacebyNode%Pt(:) = 0
      do i = 1, NbNode
         counterNumFacebyNode = counterNumFacebyNode + nbFacebyNode(i)
         FacebyNode%Pt(i + 1) = FacebyNode%Pt(i) + nbFacebyNode(i)
      enddo

      allocate (FacebyNode%Num(counterNumFacebyNode))

      ! 2nd step - filling
      nbFacebyNode(:) = 0
      do i = 1, NbFace
         nbtempFace = NodebyFace%Pt(i + 1) - NodebyFace%Pt(i)
         beginNode = NodebyFace%Pt(i) + 1

         ! Loop over the nodes of cell i
         do j = 1, nbtempFace
            loadNode = NodebyFace%Num(beginNode + j - 1)
            nbFacebyNode(loadNode) = nbFacebyNode(loadNode) + 1
            ! Number of cells surrounding node j
            FacebyNode%Num(nbFacebyNode(loadNode) + FacebyNode%Pt(loadNode)) = i
         enddo

      enddo

      deallocate (nbFacebyNode)

   end subroutine GlobalMesh_FaceByNodeGlobal

   !> \brief Make CSR FracbyNode with the number of fracture face by Node.
   !!
   !! One line corresponds to one Node, if there is no fracture in node i: \%Pt(i+1) = \%Pt   <br>
   !! Output:
   !!  FracbyNode   <br>
   !! Use:
   !!  NodebyFace, IdFace
   subroutine GlobalMesh_FracbyNode() &
      bind(C, name="GlobalMesh_frac_by_node")

      integer :: i, n, num_face, num_frac, num_node, npt
      integer, allocatable, dimension(:) :: comptNode

      allocate (comptNode(NbNode))
      comptNode(:) = 0

      ! counting
      do num_face = 1, NbFace
         if (IdFace(num_face) == -2) then ! fracface
            do n = NodebyFace%Pt(num_face) + 1, NodebyFace%Pt(num_face + 1)
               num_node = NodebyFace%Num(n)
               comptNode(num_node) = comptNode(num_node) + 1
            enddo
         endif
      enddo

      ! Filling
      FracbyNode%Nb = NbNode
      allocate (FracbyNode%Pt(NbNode + 1))
      FracbyNode%Pt(:) = 0
      do i = 1, NbNode
         FracbyNode%Pt(i + 1) = FracbyNode%Pt(i) + comptNode(i)
      enddo

      allocate (FracbyNode%Num(FracbyNode%Pt(NbNode + 1)))
      comptNode(:) = 0
      num_frac = 0
      do num_face = 1, NbFace
         if (IdFace(num_face) == -2) then ! fracface
            num_frac = num_frac + 1
            do n = NodebyFace%Pt(num_face) + 1, NodebyFace%Pt(num_face + 1)
               num_node = NodebyFace%Num(n)
               comptNode(num_node) = comptNode(num_node) + 1
               npt = FracbyNode%Pt(num_node) + comptNode(num_node)
               FracbyNode%Num(npt) = num_face
            enddo
         endif
      enddo

      deallocate (comptNode)

   end subroutine GlobalMesh_FracbyNode

   !> \brief Build NodebyWell and NodeDatabyWell after the lecture of Filemesh.
   !!
   !! The order of the storage of the node of each well is important :               <br>
   !! the parent must be stored after the son(s).
   subroutine BuildWellConnectivity(geometries, NodebyWell, NodeDatabyWell)

      type(GlobalWellGeometries), target, intent(in) :: geometries

      type(CSR), intent(out) :: NodebyWell
      type(TYPE_CSRDataNodeWell), intent(out) :: NodeDatabyWell

      integer :: i, j, ival, nb_wells, max_nb_edges_per_well, total_nb_edges
      integer, allocatable, dimension(:) :: OrderedNodes
      integer, allocatable, dimension(:) :: Parents
      integer(c_int), dimension(:), pointer :: NbEdgebyWell
      integer(c_int), dimension(:, :, :), pointer :: NumNodebyEdgebyWell

      NbEdgebyWell => geometries%NbEdgebyWell
      NumNodebyEdgebyWell => geometries%NumNodebyEdgebyWell

      nb_wells = geometries%Nb
      total_nb_edges = sum(NbEdgebyWell)

      allocate (NodebyWell%Pt(nb_wells + 1))
      allocate (NodebyWell%Num(total_nb_edges + nb_wells))
      allocate (NodeDatabyWell%Pt(nb_wells + 1))
      allocate (NodeDatabyWell%Num(total_nb_edges + nb_wells))
      allocate (NodeDatabyWell%Val(total_nb_edges + nb_wells))

      NodebyWell%Nb = nb_wells
      NodeDatabyWell%Nb = nb_wells
      NodebyWell%Num = -1
      NodeDatabyWell%Num = -1
      NodeDatabyWell%Val(:)%PtParent = -1
      NodeDatabyWell%Val(:)%RelParent = -1

      ! FIXME: Use a CSR initialization routine
      NodebyWell%Pt(1) = 0
      do i = 1, nb_wells
         NodebyWell%Pt(i + 1) = NodebyWell%Pt(i) + NbEdgebyWell(i) + 1
      end do
      NodeDatabyWell%Pt(:) = NodebyWell%Pt(:)

      ! temporary arrays
      max_nb_edges_per_well = size(NumNodebyEdgebyWell, 2)
      allocate (OrderedNodes(max_nb_edges_per_well + 1))
      allocate (Parents(max_nb_edges_per_well + 1))

      ! sort the nodes, following the well path
      do i = 1, nb_wells
         call SortWellNodes(NumNodebyEdgebyWell(1:2, 1:NbEdgebyWell(i), i), OrderedNodes, Parents)
         do j = 1, NbEdgebyWell(i) + 1
            ival = NodebyWell%Pt(i) + j
            NodebyWell%Num(ival) = OrderedNodes(j)
            NodeDatabyWell%Val(ival)%Parent = Parents(j)
         end do
      end do
      NodeDatabyWell%Num(:) = NodebyWell%Num(:)

      deallocate (OrderedNodes)
      deallocate (Parents)

   end subroutine BuildWellConnectivity

   !> \brief Sort the Nodes of each Well to store parent after each son(s)
   subroutine SortWellNodes(Edges, OrderedNodes, Parents)

      integer, dimension(:, :), intent(in) :: Edges
      integer, dimension(:), intent(out) :: OrderedNodes
      integer, dimension(:), intent(out) :: Parents
      integer, allocatable, dimension(:) :: LocalIdx
      integer, allocatable, dimension(:) :: NbEdgesbyNode
      integer :: ne, ig, pg, il, pl, i, j, io
      logical :: found
      integer :: Ierr, errcode ! used for MPI_Abort

      OrderedNodes = 0; Parents = 0

      ne = size(Edges, 2)
      allocate (LocalIdx(NbNode))
      allocate (NbEdgesbyNode(ne + 1))

      ! LocalIdx = ne+1 ! no parent detection
      ! -1: initialize
      ! 0:  head
      LocalIdx = -1
      do i = 1, ne
         LocalIdx(Edges(2, i)) = i ! indexing using son number
      end do

      ! call disp('Edges(1,...) ', Edges(1, 1:size(Edges, 2)), unit=fdGm_unit)
      ! call disp('Edges(2,...) ', Edges(2, 1:size(Edges, 2)), unit=fdGm_unit)
      ! call disp('taking a global index, what is its local index in the sons list ?: LocalIdx ', LocalIdx, unit=fdGm_unit)

      ! since we are indexing by the sons, we have to detect the heads i.e. the edges with no parents
      do i = 1, ne
         found = .false.
         do j = 1, ne
            if (Edges(1, i) == Edges(2, j)) then
               found = .true.
            end if
         end do
         if (.not. found) then
            LocalIdx(Edges(1, i)) = 0
         end if
      end do

      ! counting number of edges per node
      NbEdgesbyNode = 0
      do i = 1, ne
         ig = Edges(2, i); pg = Edges(1, i)
         il = LocalIdx(ig); pl = LocalIdx(pg)
         ! write(fdGm,*) '---->', ig, pg, il, pl
         NbEdgesbyNode(il) = NbEdgesbyNode(il) + 1
         if (pl > 0) then
            NbEdgesbyNode(pl) = NbEdgesbyNode(pl) + 1
         end if
         !write(fdGm,'(4(a3,i3))') 'ig', ig, 'il', il, 'pg', pg, 'pl', pl
         !write(fdGm,'(a22,i3,a8,i3,a1,i3,i3)') 'NbEdgesbyNode(current=', ig, ',parent=',pg ,')' , NbEdgesbyNode(il)
         !write(fdGm,*)
      end do
      ! call disp('NbEdgesbyNode ', NbEdgesbyNode, unit=fdGm_unit) ! sons connected to how many nodes ?

      io = 1
      do i = 1, ne
         ig = Edges(2, i); il = LocalIdx(ig)
         !write(fdGm,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
         !write(fdGm,*) '> ig', ig, 'NbEdgesbyNode', NbEdgesbyNode(il)
         ! call disp('OrderedNodes ', OrderedNodes, 'i3', unit=fdGm_unit)
         ! call disp('Parents      ', Parents, 'i3', unit=fdGm_unit)
         if (Edges(1, i) == Edges(2, i)) then
            !write(*,*) 'single/orphan node', Edges(1, i), 'not permitted, check your mesh !'
            call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
         end if
         if (NbEdgesbyNode(il) == 1) then
            ! is a queue or head
            !write(fdGm,*) '--> current queue', ig, 'io', io
            do ! infinite loop, we should detect pathologic cases like loops
               pg = Edges(1, il); pl = LocalIdx(pg)
               ! if(pl > 0) then
               !    write(fdGm,*) ' parent of', ig, ' is ', pg, 'connected to ', NbEdgesbyNode(pl), ' edges'
               ! else
               !    write(fdGm,*) ' parent of', ig, ' is ', pg
               ! end if
               if (pl == 0) then
                  ! head
                  ! write(fdGm,*) '(h) ig ', ig, ' pg ', pg
                  ! store this point
                  OrderedNodes(io) = Edges(2, il)
                  Parents(io) = Edges(1, il)
                  io = io + 1
                  ! store a point with a (-1) parent ...
                  OrderedNodes(io) = Edges(1, il)
                  Parents(io) = -1  ! No parents
                  io = io + 1
                  ! write(*,*) '++ (p,s)', OrderedNodes(io), Parents(io)
                  ! call disp('OrderedNodes ', OrderedNodes, 'i3', unit=fdGm_unit)
                  ! call disp('Parents      ', Parents, 'i3', unit=fdGm_unit)
                  exit
               else if (NbEdgesbyNode(pl) == 2) then
                  ! node between 2 edges, OK, continue to go up
                  OrderedNodes(io) = Edges(2, il)
                  Parents(io) = Edges(1, il)
                  ! write(fdGm,*) '(++,continue) new current node', ig
                  ! write(*,*) '++ (p,s)', OrderedNodes(io), Parents(io)
                  ! call disp('OrderedNodes ', OrderedNodes, 'i3', unit=fdGm_unit)
                  ! call disp('Parents      ', Parents, 'i3', unit=fdGm_unit)
                  io = io + 1
                  ig = pg; il = LocalIdx(ig) ! new ig, new il for next iteration
                  cycle
               else if (NbEdgesbyNode(pl) > 2) then
                  ! we have more than 2 edges connected to this node, we can quit the loop
                  NbEdgesbyNode(pl) = NbEdgesbyNode(pl) - 1
                  OrderedNodes(io) = Edges(2, il)
                  Parents(io) = Edges(1, il)
                  io = io + 1
                  ! write(*,*) '++ (p,s)', OrderedNodes(io), Parents(io)
                  ! call disp('OrderedNodes ', OrderedNodes, 'i3', unit=fdGm_unit)
                  ! call disp('Parents      ', Parents, 'i3', unit=fdGm_unit)
                  ! write(fdGm,*) '(--,break) decrement NbEdgesbyNode for parent', pg, ' and quit loop for node', ig
                  exit
               else if (NbEdgesbyNode(pl) == 0) then
                  !write(fdGm,*) '(break) orphan node ', ig, ' is not connected to a well'
                  exit
               end if
            end do
         end if
      end do

      ! some tests
      if (OrderedNodes(1) == 0) then
         ! no node has been stored for this well, error
         write (*, *) 'no node in OrderedNodes (maybe no queue detected, or cyclic well)'
         write (*, *) 'are you sure about the format of your Filemesh?'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      endif
      do i = 2, io - 1

      end do

      ! TODO, exit if no edges in %Val%Parent and/or %Num are zeross

      deallocate (LocalIdx)
      deallocate (NbEdgesbyNode)

   end subroutine SortWellNodes

   !> \brief Build global connectivity of injection and production wells.
   subroutine GlobalMesh_WellConnectivity() &
      bind(C, name="build_well_connectivity")

      call BuildWellConnectivity(GeomInj, NodebyWellInj, NodeDatabyWellInj)
      call BuildWellConnectivity(GeomProd, NodebyWellProd, NodeDatabyWellProd)
      call BuildWellConnectivity(GeomMSWell, NodebyMSWell, NodeDatabyMSWell)
      call BuildWellNodesNumbering(NodebyMSWell, MSWellNodebyNode)

   end subroutine GlobalMesh_WellConnectivity

   !> \brief Build MSWellNodebyMSWell
   !!  uses well node  numbering
   subroutine BuildWellNodesNumbering(NodebyWell, WellNodebyNode)

      type(CSR), intent(in) :: NodebyWell
      integer, dimension(:), allocatable :: WellNodebyNode
      integer :: i, s, well_node_idx, nb_wells

      nb_wells = NodebyWell%Nb
      allocate (WellNodebyNode(NbNode))
      WellNodebyNode(:) = -1

      well_node_idx = 1
      do i = 1, nb_wells
         do s = NodebyWell%Pt(i) + 1, NodebyWell%Pt(i + 1)
            WellNodebyNode(NodebyWell%Num(well_node_idx)) = well_node_idx
            well_node_idx = well_node_idx + 1
         end do
      end do

   end subroutine BuildWellNodesNumbering

   subroutine fill_CSR(ptr, indices, csrdata, c_indexing)

      integer(c_int), dimension(:), intent(in) :: ptr
      integer(c_int), dimension(:), intent(in) :: indices
      type(CSR), intent(inout) :: csrdata
      logical, optional, intent(in) :: c_indexing
      integer :: n

      n = size(ptr) - 1
      csrdata%Nb = n
      if (allocated(csrdata%Pt)) deallocate (csrdata%Pt)
      allocate (csrdata%Pt(n + 1))
      csrdata%Pt = ptr
      if (allocated(csrdata%Num)) deallocate (csrdata%Num)
      allocate (csrdata%Num(ptr(n + 1)))
      if (present(c_indexing) .and. c_indexing) then
         csrdata%Num = indices + 1
      else
         csrdata%Num = indices
      end if
      if (allocated(csrdata%Val)) deallocate (csrdata%Val)
      ! CHECKME: val is NOT reallocated... should use COC strtucture here

   end subroutine fill_CSR

   !> \todo CHECKME/IMPROVE: data is copied... should work with C structures ?!
   subroutine GlobalMesh_create_mesh(nodes, &
                                     cell_faces_ptr, cell_faces_val, &
                                     cell_nodes_ptr, cell_nodes_val, &
                                     face_nodes_ptr, face_nodes_val, &
                                     cell_id, face_id, &
                                     c_indexing)

      real(c_double), dimension(:, :), intent(in) :: nodes
      integer(c_int), dimension(:), intent(in) :: cell_faces_ptr
      integer(c_int), dimension(:), intent(in) :: cell_faces_val
      integer(c_int), dimension(:), intent(in) :: cell_nodes_ptr
      integer(c_int), dimension(:), intent(in) :: cell_nodes_val
      integer(c_int), dimension(:), intent(in) :: face_nodes_ptr
      integer(c_int), dimension(:), intent(in) :: face_nodes_val
      integer(c_int), dimension(:), intent(in) :: cell_id
      integer(c_int), dimension(:), intent(in) :: face_id
      logical, optional, intent(in) :: c_indexing
      logical :: check, use_c_indexing
      integer :: Ierr, errcode ! used for MPI_Abort

      NbNode = size(nodes, 2)
      NbCell = size(cell_id)
      NbFace = size(face_id)

      check = .True.
      if (size(nodes, 1) /= 3) check = .False.
      if (size(cell_faces_ptr) /= NbCell + 1) check = .False.
      if (size(cell_faces_val) /= cell_faces_ptr(NbCell + 1)) check = .False.
      if (size(cell_nodes_ptr) /= NbCell + 1) check = .False.
      if (size(cell_nodes_val) /= cell_nodes_ptr(NbCell + 1)) check = .False.
      if (size(face_nodes_ptr) /= NbFace + 1) check = .False.
      if (size(face_nodes_val) /= face_nodes_ptr(NbFace + 1)) check = .False.
      if (.not. check) then
         write (*, *) 'Unconsistent mesh data!'
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      if (allocated(XNode)) deallocate (XNode)
      allocate (XNode(3, NbNode))
      XNode = nodes

      use_c_indexing = present(c_indexing) .and. c_indexing
      call fill_csr(cell_faces_ptr, cell_faces_val, FacebyCell, use_c_indexing)
      call fill_csr(cell_nodes_ptr, cell_nodes_val, NodebyCell, use_c_indexing)
      call fill_csr(face_nodes_ptr, face_nodes_val, NodebyFace, use_c_indexing)

      if (allocated(IdCell)) deallocate (IdCell)
      allocate (IdCell(NbCell))
      IdCell = cell_id

      if (allocated(IdFace)) deallocate (IdFace)
      allocate (IdFace(NbFace))
      IdFace = face_id

      call DefWell_clear_global_well_geometries(GeomInj)
      call DefWell_clear_global_well_geometries(GeomProd)
      call DefWell_clear_global_well_geometries(GeomMSWell)

      call GlobalMesh_allocate_flags

   end subroutine GlobalMesh_create_mesh

end module GlobalMesh
