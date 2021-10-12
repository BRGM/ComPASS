!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

!> \brief Main subroutine LocalMesh_Make.  <br>
!! Contains the local connectivity (after partition) of all CPUs.
!!
!! Set the local connectivity once the partition has been done.
!! Determine if the objects (Cells, Faces, Nodes, Wells) are own or ghost.
!! The main subroutine is LocalMesh_Make, called only by proc master,
!! then the informations are send in subroutine MeshSchema_make.
module LocalMesh

#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use CommonType
   use CommonMPI
   use GlobalMesh
   use DefModel
   use DefWell
#else
   use iso_c_binding, only: c_int
   use CommonType, only: &
      CSR, FamilyDOFId, FamilyDOFIdCOC, Array1IdNode, &
      ARRAY1dble, ARRAY3dble, ARRAY2Int, ARRAY1Int8, ARRAY1Int, ARRAY2dble, &
      CommonType_deallocCSR, copy_sparsity_pattern
   use CommonMPI, only: Ncpus, CommonMPI_abort

   use GlobalMesh, only: &
      Mesh_xmax, Mesh_xmin, Mesh_ymin, Mesh_ymax, Mesh_zmax, Mesh_zmin, &
      CellFlags, FaceFlags, &
      CellbyCell, FacebyCell, NodebyCell, &
      CellbyFace, NodebyFace, &
      CellbyNode, &
      IdNode, IdFace, IdCell, &
      NodebyWellInj, NodebyWellProd, &
      NbWellInj, NbFace, NbCell, NbNode, NbWellProd, &
      CellThermalSource, FracThermalSource, CondThermalFrac, &
      PermFrac, PorositeCell, PorositeFrac, NodeFlags, &
      CellTypes, FaceTypes, &
      NodeRocktype, CellRocktype, FracRocktype, &
      XNode, PermCell, CondThermalCell

   ! use PartitionMesh, only: &
   ! ProcbyCell

   use DefModel, only: &
      IndThermique

   use DefWell, only: &
      TYPE_CSRDataNodeWell, WellData_type, &
      DataWellInj, DataWellProd, NodeDatabyWellProd, NodeDatabyWellInj
#endif

   ! 1. FacebyCellRes_Ncpus, NodebyFaceRes_Ncpus,
   !    NodebyCellRes_Ncpus, CellbyFaceRes_Ncpus
   !    CellbyProc, FacebyProc, NodebyProc

   ! 2. FacebyCellLocal_Ncpus, NodebyFaceLocal_Ncpus,
   !    NodebyCellLocal_Ncpus, CellbyFaceLocal_Ncpus

   ! 3. CellbyNodeOwn_Ncpus, NodebyNodeOwn_Ncpus, FracbyNodeOwn_Ncpus
   !    CellbyFracOwn_Ncpus, FacebyFracOwn_Ncpus, NodebyFracOwn_Ncpus
   !    FracbyCellLocal_Ncpus

   implicit none

   ! Outputs

   !! Local Mesh info res/own, res=own+ghost
   !! Rq: we will define NbCellOwn_Ncpus in MeshSchema.F90,
   !!     we want all the vectors used in schema are defined in MeshSchema.F90,
   !!     thus, here we use a different name: "S" is added in NbCellOwn"S"_Ncpus
   integer, dimension(:), allocatable, public :: &
      NbCellResS_Ncpus, & !< Number of Cells res=own+ghost for all CPUs
      NbCellOwnS_Ncpus, & !< Number of Cells own for all CPUs
      NbFaceResS_Ncpus, & !< Number of Faces res=own+ghost for all CPUs
      NbFaceOwnS_Ncpus, & !< Number of Faces own for all CPUs
      NbNodeResS_Ncpus, & !< Number of Nodes res=own+ghost for all CPUs
      NbNodeOwnS_Ncpus, & !< Number of Nodes own for all CPUs
      NbFracResS_Ncpus, & !< Number of Fracture Faces res=own+ghost for all CPUs
      NbFracOwnS_Ncpus    !< Number of Fracture Faces own for all CPUs

   ! Well
   integer, dimension(:), allocatable, public :: &
      NbWellInjResS_Ncpus, & !< Number of Injection Wells res=own+ghost for all CPUs
      NbWellInjOwnS_Ncpus, & !< Number of Injection Wells own for all CPUs
      NbWellProdResS_Ncpus, & !< Number of Production Wells res=own+ghost for all CPUs
      NbWellProdOwnS_Ncpus    !< Number of Production Wells own for all CPUs

   ! Well connectivities in global number
   type(CSR), dimension(:), allocatable, public :: &
      NodebyWellInjRes_Ncpus, & !< Node (global number) of injection well res=own+ghost for all CPUs
      NodebyWellProdRes_Ncpus   !< Node (global number) of production well res=own+ghost for all CPUs

   ! Well connectivites in index (local)
   type(CSR), dimension(:), allocatable, public :: &
      NodebyWellInjLocal_Ncpus, & !< Node (local number) of injection well res=own+ghost for all CPUs
      NodebyWellProdLocal_Ncpus   !< Node (local number) of production well res=own+ghost for all CPUs

   Double precision, protected :: &
      meshSizeS_xmax, & !< Temporary vector to send global size of mesh xmax
      meshSizeS_xmin, & !< Temporary vector to send global size of mesh xmin
      meshSizeS_ymax, & !< Temporary vector to send global size of mesh ymax
      meshSizeS_ymin, & !< Temporary vector to send global size of mesh ymin
      meshSizeS_zmax, & !< Temporary vector to send global size of mesh zmax
      meshSizeS_zmin    !< Temporary vector to send global size of mesh zmin

   !! Local mesh and connectivity in num (global) for all procs
   type(CSR), dimension(:), allocatable, public :: &
      FacebyCellRes_Ncpus, & !< Connectivity FacebyCell (in num global) for local Cells for all procs
      NodebyFaceRes_Ncpus, & !< Connectivity NodebyFace (in num global) for local Faces for all procs
      NodebyCellRes_Ncpus, & !< Connectivity NodebyCell (in num global) for local Cells for all procs
      CellbyFaceRes_Ncpus    !< Connectivity CellbyFace (in num global) for local Faces for all procs

   !! Local mesh and connectivity in num (local) for all procs
   !! Using *Res_Ncpus to calculate the following vectors
   type(CSR), dimension(:), allocatable, public :: &
      FacebyCellLocal_Ncpus, &  !< Connectivity FacebyCell (in num local) for local Cells for all procs
      NodebyFaceLocal_Ncpus, &  !< Connectivity NodebyFace (in num local) for local Faces for all procs
      CellbyFaceLocal_Ncpus     !< Connectivity CellbyFace (in num local) for local Faces for all procs

   !! Local elements by proc in num (global) in order
   !!  | elements own | elements ghost (which are own for proc i) | elements ghost (which are own for proc j>i) |
   type(CSR), dimension(:), allocatable, public :: &
      CellbyProc, & !< Local Cells by proc in num global (CSR)
      FacebyProc, & !< Local Faces by proc in num global (CSR)
      NodebyProc, & !< Local Nodes by proc in num global (CSR)
      FracbyProc    !< Local Fracture faces by proc in num global (CSR)

   !! Local wells by proc in num (global)
   !!  | wells own | wells ghost (which are own for proc i) | wells ghost (which are own for proc j>i) | ...
   type(CSR), dimension(:), allocatable, public :: &
      WellInjbyProc, & !< Local injection Wells by proc in num global (CSR)
      WellProdbyProc   !< Local production Wells by proc in num global (CSR)

   !! DataWellRes_Ncpus(num_welllocal,num_proc)
   type(WellData_type), allocatable, dimension(:, :), public :: DataWellInjRes_Ncpus !< Data of injection well (Radius,...) for local well
   type(WellData_type), allocatable, dimension(:, :), public :: DataWellProdRes_Ncpus !< Data of production well (Radius,...) for local well

   !! NumNode/Frac/WellbyProc_Ncpus(ip): node/frac/well num in the proc that this object is own
   type(FamilyDOFIdCOC), dimension(:), allocatable, public :: &
      NumNodebyProc_Ncpus, &
      NumFracbyProc_Ncpus, &
      NumWellInjbyProc_Ncpus, &
      NumWellProdbyProc_Ncpus

   !! IdFaceRes_Ncpus
   !   = -2 if Frac
   !   >= 1 if Dir
   !   = 0 else
   type(Array1Int), dimension(:), allocatable, public :: &
      IdFaceRes_Ncpus

   !! IdNodeRes_Ncpus
   type(Array1IdNode), dimension(:), allocatable, public :: &
      IdNodeRes_Ncpus

   !! IdCell
   type(ARRAY1Int), dimension(:), allocatable, public :: &
      IdCellRes_Ncpus

   !! The following vectors are used to the strucutre of Jacobian
   type(CSR), dimension(:), allocatable, public :: &
      CellbyNodeOwn_Ncpus, &  ! from NodebyCellLocal_Ncpus
      NodebyNodeOwn_Ncpus, &  ! from NodebyCellLocal_Ncpus, CellbyNodeOwn_Ncpus
      FracbyNodeOwn_Ncpus, &  ! from NodebyFaceLocal_Ncpus
      !
      NodebyCellLocal_Ncpus, &  ! from NodebyCellRes_Ncpus
      FracbyCellLocal_Ncpus, &  ! from FacebyCellLocal_Ncpus
      !
      NodebyFracOwn_Ncpus, &  ! from NodebyFaceLocal_Ncpus
      ! WARNING these are not the fracture nodes but the set of nodes
      !         of the two cells on each side of the fracture
      CellbyFracOwn_Ncpus, &  ! from CellbyFaceLocal_Ncpus
      FracbyFracOwn_Ncpus     ! from FracbyCellLocal_Ncpus, CellbyFracOwn_Ncpus

   !! The following vectors are used to the strucutre of Jacobian
   type(CSR), dimension(:), allocatable, public :: &
      WellInjbyNodeOwn_Ncpus, &  ! numero (local) of well inj connected to this node own
      WellProdbyNodeOwn_Ncpus    ! numero (local) of well prod connected to this node own

   type(CSR), dimension(:), allocatable, public :: &
      FacebyNodeOwn_Ncpus ! only used to make FracbyNodeOwn_Ncpus

   !! Local frac face to face
   type(ARRAY1Int), dimension(:), allocatable, public :: &
      FracToFaceLocal_Ncpus, &
      FaceToFracLocal_Ncpus

   !! XNode Local
   type(ARRAY2dble), dimension(:), allocatable, target, public :: &
      XNodeRes_Ncpus

   type(ARRAY1Int), dimension(:), allocatable, target, public :: &
      NodeFlags_Ncpus, &
      CellFlags_Ncpus, &
      FaceFlags_Ncpus

   type(ARRAY1Int8), dimension(:), allocatable, target, public :: &
      CellTypes_Ncpus, &
      FaceTypes_Ncpus

   ! Rocktype
   type(ARRAY2Int), dimension(:), allocatable, target, public :: &
      NodeRocktype_Ncpus, &
      CellRocktype_Ncpus, &
      FracRocktype_Ncpus

   ! Porosite
   type(ARRAY1dble), dimension(:), allocatable, public :: &
      PorositeCell_Ncpus, &
      PorositeFrac_Ncpus

   ! Permeability
   type(ARRAY3dble), dimension(:), allocatable, public :: &
      PermCellLocal_Ncpus
   type(ARRAY1dble), dimension(:), allocatable, public :: &
      PermFracLocal_Ncpus

#ifdef _THERMIQUE_
   ! Thermal conductivity
   type(ARRAY3dble), dimension(:), allocatable, public :: &
      CondThermalCellLocal_Ncpus
   type(ARRAY1dble), dimension(:), allocatable, public :: &
      CondThermalFracLocal_Ncpus

   ! Thermal source
   TYPE(ARRAY1dble), DIMENSION(:), ALLOCATABLE, PUBLIC :: &
      CellThermalSource_Ncpus, &
      FracThermalSource_Ncpus
#endif

   ! Data of Node by Well (Parent, PtParent, WID, WIF)
   type(TYPE_CSRDataNodeWell), dimension(:), allocatable, public :: &
      NodeDatabyWellInjLocal_Ncpus, &
      NodeDatabyWellProdLocal_Ncpus

   ! tmp vectors used to transform from num (global) to num (local)
   integer, dimension(:), allocatable, private :: &
      localbyGlobalCell, &
      localbyGlobalFace, &
      localbyGlobalNode

   public :: &
      LocalMesh_Make, &
      LocalMesh_Free

   private :: &
      LocalMesh_CellbyProc, & ! make CellbyProc(ip1) NbCellResS/OwnS_Ncpus, ip1=ip+1, ip is number of proc
      LocalMesh_FacebyCellRes, & ! make FacebyCellRes_Ncpus(ip1)
      LocalMesh_FacebyProc, & ! make FacebyProc(ip1)
      LocalMesh_NodebyCellRes, & ! make NodebyCellRes_Ncpus(ip1)
      LocalMesh_NodebyProc, & ! make NodebyProc(ip1), NbNodeRes/OwnS_Ncpus(ip1)
      LocalMesh_NodebyFaceRes, & ! make NodebyFaceRes_Ncpus(ip1)
      LocalMesh_CellbyFaceRes, & ! make CellbyFaceRes_Ncpus(ip1)
      !
      LocalMesh_WellbyProc, & ! make WellInjbyProc(ip1), DataWellInjRes_Ncpus(ip1)   (and WellProd)
      LocalMesh_NodebyWellRes, & ! make NodebyWellInjRes_Ncpus(ip1)   (and WellProd)
      !
      LocalMesh_GlobalToLocal, & ! make vectors used to trasnform num (global) tp num(local)
      !
      LocalMesh_FacebyCellLocal, & ! make FacebyCellLocal_Ncpus(ip1)
      LocalMesh_NodebyCellLocal, & ! make NodebyCellLocal_Ncpus(ip1)
      LocalMesh_NodebyFaceLocal, & ! make NodebyFaceLocal_Ncpus(ip1)
      LocalMesh_CellbyFaceLocal, & ! make CellbyNodeLocal_Ncpus(ip1)
      !
      LocalMesh_NodebyWellLocal, & ! make NodebyWellInjLocal_Ncpus(ip1)   (and WellProd)
      !
      LocalMesh_IdCellRes, & ! make IdCellRes_Ncpus(ip1)
      LocalMesh_IdFaceRes, & ! make IdFaceRes_Ncpus(ip1)
      LocalMesh_IdNodeRes, & ! make IdNodeRes_Ncpus(ip1), NbNodeNotDirRes/OwnS_Ncpus(ip1)
      !
      LocalMesh_FracbyCellLocal, & ! make FracbyCellLocal_Ncpus(ip1)
      LocalMesh_FracbyProc, & ! make FracbyProc(ip1)
      !
      ! NodebyCellLocal-> CellbyNodeOwn
      !                 + NodebyCellLocal -> NodebyNodeOwn
      ! NodebyFaceLocal-> FacebyNodeOwn -> FracbyNodeOwn
      LocalMesh_CellbyNodeOwn, & ! make NodebyCellOwn_Ncpus(ip1)
      LocalMesh_NodebyNodeOwn, & ! make NodebyNodeOwn_Ncpus(ip1)
      LocalMesh_FacebyNodeOwn, & ! make FacebyNodeOwn_Ncpus(ip1), only used to make FracbyNodeOwn_Ncpus(ip1)
      LocalMesh_FracbyNodeOwn, & ! make FracbyNodeOwn_Ncpus(ip1)
      !
      LocalMesh_WellbyNodeOwn, & ! make WellInjbyNodeOwn_Ncpus(ip1)     (and WellProd)
      !
      ! CellbyFaceLocal-> CellbyFracOwn
      !                 + NodebyCellLocal -> NodebyFracOwn
      !                 + FracbyCellLocal -> FracbyFracOwn
      LocalMesh_CellbyFracOwn, & ! make CellbyFracOwn_Ncpus(ip1)
      LocalMesh_NodebyFracOwn, & ! make NodebyFracOwn_Ncpus(ip1)
      !
      LocalMesh_FracbyFracOwn, & ! make FracbyFracOwn_Ncpus(ip1)
      !
      LocalMesh_XNodeRes, & ! make XNodeRes_Ncpus(ip1)
      !
      LocalMesh_Porosite, & ! porosite
      LocalMesh_Perm, & ! Permeability
      !
      LocalMesh_NodeDatabyWellLocal, & ! make NodeDatabyWellLocal_Ncpus(ip1)
      !
      LocalMesh_FracToFace, & ! make FracToFaceLocal(ip1)
      LocalMesh_FaceToFrac

contains

   !> \brief Main soubroutine of the file
   subroutine LocalMesh_Make(ProcbyCell)

      integer, dimension(:), intent(in) :: ProcbyCell

      integer :: i
      integer, allocatable, dimension(:) :: proc_color, node_color, well_color

      !> Allocate LocalMesh entities
      ! local mesh info res(own+ghost)
      allocate (NbCellResS_Ncpus(Ncpus))
      allocate (NbFaceResS_Ncpus(Ncpus))
      allocate (NbNodeResS_Ncpus(Ncpus))
      allocate (NbFracResS_Ncpus(Ncpus))

      ! local info own
      allocate (NbCellOwnS_Ncpus(Ncpus))
      allocate (NbFaceOwnS_Ncpus(Ncpus))
      allocate (NbNodeOwnS_Ncpus(Ncpus))
      allocate (NbFracOwnS_Ncpus(Ncpus))

      ! local info well res and own
      allocate (NbWellInjResS_Ncpus(Ncpus))
      allocate (NbWellInjOwnS_Ncpus(Ncpus))
      allocate (NbWellProdResS_Ncpus(Ncpus))
      allocate (NbWellProdOwnS_Ncpus(Ncpus))

      ! local mesh and connectivities in num (global)
      allocate (FacebyCellRes_Ncpus(Ncpus))
      allocate (NodebyCellRes_Ncpus(Ncpus))
      allocate (CellbyFaceRes_Ncpus(Ncpus))
      allocate (NodebyFaceRes_Ncpus(Ncpus))
      ! local well connectivity in num (global)
      allocate (DataWellInjRes_Ncpus(NbWellInj, Ncpus))
      allocate (DataWellProdRes_Ncpus(NbWellProd, Ncpus))
      allocate (NodebyWellInjRes_Ncpus(Ncpus))
      allocate (NodebyWellProdRes_Ncpus(Ncpus))

      ! local mesh and connectivities in num (local) own+ghost
      allocate (FacebyCellLocal_Ncpus(Ncpus))
      allocate (NodebyFaceLocal_Ncpus(Ncpus))
      allocate (NodebyCellLocal_Ncpus(Ncpus))
      allocate (CellbyFaceLocal_Ncpus(Ncpus))
      allocate (FracbyCellLocal_Ncpus(Ncpus))
      ! local well connectivity in num (local)
      allocate (NodebyWellInjLocal_Ncpus(Ncpus))
      allocate (NodebyWellProdLocal_Ncpus(Ncpus))

      ! local mesh and connectivities in num (local) own
      allocate (CellbyNodeOwn_Ncpus(Ncpus))
      allocate (NodebyNodeOwn_Ncpus(Ncpus))
      allocate (FacebyNodeOwn_Ncpus(Ncpus))
      allocate (FracbyNodeOwn_Ncpus(Ncpus))

      allocate (NodebyFracOwn_Ncpus(Ncpus))
      allocate (CellbyFracOwn_Ncpus(Ncpus))
      allocate (FracbyFracOwn_Ncpus(Ncpus))

      allocate (WellInjbyNodeOwn_Ncpus(Ncpus))
      allocate (WellProdbyNodeOwn_Ncpus(Ncpus))

      allocate (XNodeRes_Ncpus(Ncpus))
      allocate (NodeFlags_Ncpus(Ncpus))
      allocate (CellFlags_Ncpus(Ncpus))
      allocate (FaceFlags_Ncpus(Ncpus))
      allocate (CellTypes_Ncpus(Ncpus))
      allocate (FaceTypes_Ncpus(Ncpus))

      allocate (NodeRocktype_Ncpus(Ncpus))
      allocate (CellRocktype_Ncpus(Ncpus))
      allocate (FracRocktype_Ncpus(Ncpus))

      ! local element by proc
      allocate (CellbyProc(Ncpus))
      allocate (FacebyProc(Ncpus))
      allocate (NodebyProc(Ncpus))
      allocate (FracbyProc(Ncpus))
      allocate (WellInjbyProc(Ncpus))
      allocate (WellProdbyProc(Ncpus))

      ! IdFace/Node
      allocate (IdFaceRes_Ncpus(Ncpus))
      allocate (IdNodeRes_Ncpus(Ncpus))
      allocate (IdCellRes_Ncpus(Ncpus))

      ! Frac to/from Face
      allocate (FracToFaceLocal_Ncpus(Ncpus))
      allocate (FaceToFracLocal_Ncpus(Ncpus))

      ! porosity
      allocate (PorositeCell_Ncpus(Ncpus))
      allocate (PorositeFrac_Ncpus(Ncpus))

      ! permeability
      allocate (PermCellLocal_Ncpus(Ncpus))
      allocate (PermFracLocal_Ncpus(Ncpus))

#ifdef _THERMIQUE_
      ! Thermal conductivity
      allocate (CondThermalCellLocal_Ncpus(Ncpus))
      allocate (CondThermalFracLocal_Ncpus(Ncpus))

      ! Thermal source
      allocate (CellThermalSource_Ncpus(Ncpus))
      allocate (FracThermalSource_Ncpus(Ncpus))
#endif

      ! Data of Node by Well
      allocate (NodeDatabyWellInjLocal_Ncpus(Ncpus))
      allocate (NodeDatabyWellProdLocal_Ncpus(Ncpus))

      ! tmp vectors used to transform from num (global) to num (local)
      allocate (localbyGlobalCell(NbCell))
      allocate (localbyGlobalFace(NbFace))
      allocate (localbyGlobalNode(NbNode))

      ! work arrays
      allocate (node_color(NbNode))
      allocate (proc_color(Ncpus))
      allocate (well_color(max(NbWellInj, NbWellProd)))

      do i = 0, Ncpus - 1

         !> Set Mesh and connectivities in num (global)
         call LocalMesh_CellbyProc(i, ProcbyCell) ! CellbyProc(ip1)
         call LocalMesh_FacebyCellRes(i)    ! FacebyCellRes_Ncpus(ip1)
         call LocalMesh_FacebyProc(i, ProcbyCell) ! FacebyProc(ip1)
         call LocalMesh_NodebyCellRes(i)    ! NodebyCellRes_Ncpus(ip1)
         call LocalMesh_NodebyProc(i, ProcbyCell) ! NodebyProc(ip1)
         call LocalMesh_NodebyFaceRes(i)    ! NodebyFaceRes_Ncpus(ip1)
         call LocalMesh_CellbyFaceRes(i)    ! CellbyFaceRes_Ncpus(ip1)

         !> Set Mesh and connectivites in num (local)

         ! Wells
         call LocalMesh_WellbyProc(i, node_color, proc_color, well_color)  ! WellInjbyProc(ip1)   (and WellProd)
         call LocalMesh_NodebyWellRes(i)  ! NodebyWellInjRes_Ncpus(ip1)   (and WellProd)

         ! vectors used to transform num (global) to num (local)
         call LocalMesh_GlobalToLocal(i)    ! localbyGlobalNode

         call LocalMesh_FacebyCellLocal(i)  ! FacebyCellLocal_Ncpus(ip1)
         call LocalMesh_NodebyCellLocal(i)  ! NodebyCellLocal_Ncpus(ip1)
         call LocalMesh_NodebyFaceLocal(i)  ! NodebyFaceLocal_Ncpus(ip1)
         call LocalMesh_CellbyFaceLocal(i)  ! CellbyNodeLocal_Ncpus(ip1)

         ! Wells
         call LocalMesh_NodebyWellLocal(i)  ! NodebyWellInjLocal_Ncpus(ip1)   (and WellProd)

         call LocalMesh_IdFaceRes(i)        ! IdFaceRes_Ncpus(ip1)
         call LocalMesh_IdNodeRes(i)        ! IdNodeRes_Ncpus(ip1)
         call LocalMesh_IdCellRes(i)        ! IdCellRes_Ncpus(ip1)

         call LocalMesh_FracbyCellLocal(i)  ! FracbyCellLocal_Ncpus(ip1)
         call LocalMesh_FracbyProc(i)       ! FracbyProc(ip1)

         call LocalMesh_CellbyNodeOwn(i)    ! NodebyCellOwn_Ncpus(ip1)
         call LocalMesh_NodebyNodeOwn(i)    ! NodebyNodeOwn_Ncpus(ip1)
         call LocalMesh_FacebyNodeOwn(i)    ! FacebyNodeOwn_Ncpus(ip1), only used to make FracbyNodeOwn_Ncpus(ip1)
         call LocalMesh_FracbyNodeOwn(i)    ! FracbyNodeOwn_Ncpus(ip1)

         ! Wells
         call LocalMesh_WellbyNodeOwn(i)    ! WellInjbyNodeOwn_Ncpus(ip1)     (and WellProd)

         call LocalMesh_CellbyFracOwn(i)    ! CellbyFracOwn_Ncpus(ip1)
         call LocalMesh_NodebyFracOwn(i)    ! NodebyFracOwn_Ncpus(ip1)
         call LocalMesh_FracbyFracOwn(i)    ! FracbyFracOwn_Ncpus(ip1)

         ! X node
         call LocalMesh_XNodeRes(i)         ! XNodeRes_Ncpus(ip1)

         ! Flags
         call LocalMesh_Flags(i)

         ! Rocktype
         call LocalMesh_Rocktype(i)

         ! porosity
         call LocalMesh_Porosite(i)         ! porosity

         ! permeability
         call LocalMesh_Perm(i)

#ifdef _THERMIQUE_
         ! permeability
         call LocalMesh_CondThermal(i)

         CALL LocalMesh_ThermalSource(i)
#endif

         ! Data Node of well
         call LocalMesh_NodeDatabyWellLocal(i)  ! NodeDatabyWellLocal_Ncpus(ip1)

         ! Face to/from Frac
         call LocalMesh_FaceToFrac(i)       ! FaceToFracLocal_Ncpus(ip1)
         call LocalMesh_FracToFace(i)       ! FracToFaceLocal_Ncpus(ip1)

         !> Free mesh and connectivities in num (global) of this proc
         call CommonType_deallocCSR(FacebyCellRes_Ncpus(i + 1))
         call CommonType_deallocCSR(NodebyCellRes_Ncpus(i + 1))
         call CommonType_deallocCSR(CellbyFaceRes_Ncpus(i + 1))
         call CommonType_deallocCSR(NodebyFaceRes_Ncpus(i + 1))
         call CommonType_deallocCSR(NodebyWellInjRes_Ncpus(i + 1))
         call CommonType_deallocCSR(NodebyWellProdRes_Ncpus(i + 1))

      end do

      !> Free mesh and connectivities in num (global)
      deallocate (FacebyCellRes_Ncpus)
      deallocate (NodebyCellRes_Ncpus)
      deallocate (CellbyFaceRes_Ncpus)
      deallocate (NodebyFaceRes_Ncpus)
      deallocate (NodebyWellInjRes_Ncpus)
      deallocate (NodebyWellProdRes_Ncpus)

      ! Free tmp vectors
      deallocate (localbyGlobalCell)
      deallocate (localbyGlobalFace)
      deallocate (localbyGlobalNode)
      deallocate (well_color)
      deallocate (proc_color)
      deallocate (node_color)

      call LocalMesh_compute_local_info(NodebyProc, NbNode, NbNodeOwnS_Ncpus, NumNodebyProc_Ncpus)
      call LocalMesh_compute_local_info(FracbyProc, NbFace, NbFracOwnS_Ncpus, NumFracbyProc_Ncpus)
      call LocalMesh_compute_local_info(WellInjbyProc, NbWellInj, NbWellInjOwnS_Ncpus, NumWellInjbyProc_Ncpus)
      call LocalMesh_compute_local_info(WellProdbyProc, NbWellProd, NbWellProdOwnS_Ncpus, NumWellProdbyProc_Ncpus)

      ! Mesh Size
      meshSizeS_xmax = Mesh_xmax
      meshSizeS_xmin = Mesh_xmin
      meshSizeS_ymax = Mesh_ymax
      meshSizeS_ymin = Mesh_ymin
      meshSizeS_zmax = Mesh_zmax
      meshSizeS_zmin = Mesh_zmin

   end subroutine LocalMesh_Make

   ! Free local elements by proc
   subroutine LocalMesh_Free

      integer :: i

      do i = 1, Ncpus
         call CommonType_deallocCSR(CellbyProc(i))
         call CommonType_deallocCSR(FacebyProc(i))
         call CommonType_deallocCSR(NodebyProc(i))
         call CommonType_deallocCSR(FracbyProc(i))
         call CommonType_deallocCSR(WellInjbyProc(i))
         call CommonType_deallocCSR(WellProdbyProc(i))
      end do

      deallocate (CellbyProc)
      deallocate (FacebyProc)
      deallocate (NodebyProc)
      deallocate (FracbyProc)
      deallocate (WellInjbyProc)
      deallocate (WellProdbyProc)
      deallocate (DataWellInjRes_Ncpus)
      deallocate (DataWellProdRes_Ncpus)

   end subroutine LocalMesh_Free

   ! Output:
   !   CellbyProc(ip)
   !   NbCellRes(ip), NbCellOwn(ip)
   ! Use:
   !   CellbyCell, ProcbyCell
   !> \brief Stores own+ghost cells of proc ip in order:
   !!   | own cells | ghost cells (own for proc i) | ghost cells (own for proc j) | ...
   subroutine LocalMesh_CellbyProc(ip, ProcbyCell)

      integer, dimension(:), intent(in) :: ProcbyCell

      ! input
      integer, intent(in) :: ip
      integer :: ip1 ! ip1 = ip + 1 ip=0,1,..., ip1 used for array

      ! tmp
      integer :: cellv, j, k, procCV
      integer :: iv
      integer :: nbProcVoisinCell

      integer, allocatable, dimension(:) :: &
         voisinTemp, colorProc, colorCell, cptCellProcVois

      ip1 = ip + 1
      allocate (colorProc(Ncpus))
      allocate (voisinTemp(Ncpus))
      allocate (colorCell(NbCell))

      nbProcVoisinCell = 0
      voisinTemp(:) = 0
      colorProc(:) = -1
      colorProc(ip + 1) = 0

      ! calcul of nbProcVoisinCell
      do k = 1, NbCell
         if (ProcbyCell(k) == ip) then
            do j = CellbyCell%Pt(k) + 1, CellbyCell%Pt(k + 1)
               procCV = ProcbyCell(CellbyCell%Num(j))
               if (colorProc(procCV + 1) == -1) then
                  nbProcVoisinCell = nbProcVoisinCell + 1
                  colorProc(procCV + 1) = 0
                  voisinTemp(nbProcVoisinCell) = procCV
               endif
            enddo
         endif
      enddo
      nbProcVoisinCell = nbProcVoisinCell + 1 ! + proc courant
      CellbyProc(ip1)%Nb = nbProcVoisinCell

      allocate (CellbyProc(ip1)%Val(CellbyProc(ip1)%Nb))

      CellbyProc(ip1)%Val(1) = ip
      CellbyProc(ip1)%Val(2:CellbyProc(ip1)%Nb) = voisinTemp(1:nbProcVoisinCell - 1)

      allocate (CellbyProc(ip1)%Pt(nbProcVoisinCell + 1))
      CellbyProc(ip1)%Pt(:) = 0
      colorCell(:) = -1

      do k = 1, NbCell
         if (ProcbyCell(k) == ip) then

            ! this cell is own for ip, this cell is stored in the first line of CellbyProc
            CellbyProc(ip1)%Pt(2:nbProcVoisinCell + 1) = CellbyProc(ip1)%Pt(2:nbProcVoisinCell + 1) + 1

            ! ghost cells
            do j = CellbyCell%Pt(k) + 1, CellbyCell%Pt(k + 1)
               cellv = CellbyCell%Num(j)
               procCV = ProcbyCell(cellv)
               if (procCV /= ip) then
                  if (colorCell(cellv) == -1) then

                     do iv = 2, nbProcVoisinCell
                        if (procCV == CellbyProc(ip1)%Val(iv)) then
                           colorCell(cellv) = 0
                           CellbyProc(ip1)%Pt(iv + 1:nbProcVoisinCell + 1) = CellbyProc(ip1)%Pt(iv + 1:nbProcVoisinCell + 1) + 1
                        endif
                     enddo

                  endif
               endif
            enddo
         endif
      enddo

      allocate (CellbyProc(ip1)%Num(CellbyProc(ip1)%Pt(nbProcVoisinCell + 1)))
      allocate (cptCellProcVois(nbProcVoisinCell))
      cptCellProcVois(:) = 0 !counter of cells per proc

      colorCell(:) = -1
      do k = 1, NbCell
         if (ProcbyCell(k) == ip) then

            ! this cell is own for ip, this cell is stored in the first line of CellbyProc
            cptCellProcVois(1) = cptCellProcVois(1) + 1
            CellbyProc(ip1)%Num(CellbyProc(ip1)%Pt(1) + cptCellProcVois(1)) = k

            ! ghost cells
            do j = CellbyCell%Pt(k) + 1, CellbyCell%Pt(k + 1)
               cellv = CellbyCell%Num(j)
               procCV = ProcbyCell(cellv)
               if (procCV /= ip) then
                  if (colorCell(cellv) == -1) then
                     do iv = 2, nbProcVoisinCell
                        if (procCV == CellbyProc(ip1)%Val(iv)) then
                           colorCell(cellv) = 0
                           cptCellProcVois(iv) = cptCellProcVois(iv) + 1
                           CellbyProc(ip1)%Num(CellbyProc(ip1)%Pt(iv) + cptCellProcVois(iv)) = cellv
                        endif
                     enddo
                  endif
               endif
            enddo
         endif
      enddo

      deallocate (cptCellProcVois)
      deallocate (voisinTemp)
      deallocate (colorProc)
      deallocate (colorCell)

      NbCellResS_Ncpus(ip1) = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1) ! NbCellResS_Ncpus
      NbCellOwnS_Ncpus(ip1) = CellbyProc(ip1)%Pt(2) - Cellbyproc(ip1)%Pt(1) ! NbCellOwnS_Ncpus

   end subroutine LocalMesh_CellbyProc

   ! Output:
   !  FacebyCellRes_Ncpus(ip)
   ! Use:
   !  CellbyProc, FacebyCell
   !> \brief Stores the global number of faces sourrounding each cell (own and ghost) of proc ip
   subroutine LocalMesh_FacebyCellRes(ip)

      integer, intent(in) :: ip
      integer :: ip1

      ! tmp
      integer :: numCellRes, Vsize
      integer :: iv, k

      ip1 = ip + 1

      FacebyCellRes_Ncpus(ip1)%Nb = NbCellResS_Ncpus(ip1)
      allocate (FacebyCellRes_Ncpus(ip1)%Pt(NbCellResS_Ncpus(ip1) + 1))

      FacebyCellRes_Ncpus(ip1)%Pt(1) = 0

      ! Counting and Filling of %Pt
      do iv = 1, CellbyProc(ip1)%Nb
         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)
            numCellRes = CellbyProc(ip1)%Num(k)

            FacebyCellRes_Ncpus(ip1)%Pt(k + 1) = FacebyCellRes_Ncpus(ip1)%Pt(k) &
                                                 + FacebyCell%Pt(numCellRes + 1) - FacebyCell%Pt(numCellRes)
         enddo
      enddo

      Vsize = FacebyCellRes_Ncpus(ip1)%Pt(NbCellResS_Ncpus(ip1) + 1)
      allocate (FacebyCellRes_Ncpus(ip1)%Num(Vsize))

      ! Filling of %Num
      do iv = 1, CellbyProc(ip1)%Nb
         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)
            numCellRes = CellbyProc(ip1)%Num(k)

            FacebyCellRes_Ncpus(ip1)%Num(FacebyCellRes_Ncpus(ip1)%Pt(k) + 1:FacebyCellRes_Ncpus(ip1)%Pt(k + 1)) &
               = FacebyCell%Num(FacebyCell%Pt(numCellRes) + 1:FacebyCell%Pt(numCellRes + 1))
         enddo
      enddo

   end subroutine LocalMesh_FacebyCellRes

   ! Output:
   !  NodebyCellRes_Ncpus(ip)
   ! Use:
   !  CellbyProc, NodebyCell
   !> \brief Stores the global number of nodes sourrounding each cell (own and ghost) of proc ip
   subroutine LocalMesh_NodebyCellRes(ip)

      integer, intent(in) :: ip
      integer :: ip1

      integer :: numCellRes, Vsize, Nb
      integer :: iv, k

      ip1 = ip + 1

      Nb = NbCellResS_Ncpus(ip1)
      NodebyCellRes_Ncpus(ip1)%Nb = Nb
      allocate (NodebyCellRes_Ncpus(ip1)%Pt(Nb + 1))

      NodebyCellRes_Ncpus(ip1)%Pt(:) = 0

      ! Counting and Filling of %Pt
      do iv = 1, CellbyProc(ip1)%Nb

         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)

            numCellRes = CellbyProc(ip1)%Num(k)
            NodebyCellRes_Ncpus(ip1)%Pt(k + 1) = NodebyCellRes_Ncpus(ip1)%Pt(k) &
                                                 + NodebyCell%Pt(numCellRes + 1) - NodebyCell%Pt(numCellRes)
         enddo

      enddo

      Vsize = NodebyCellRes_Ncpus(ip1)%Pt(NbCellResS_Ncpus(ip1) + 1)
      allocate (NodebyCellRes_Ncpus(ip1)%Num(Vsize))

      ! Filling of %Num
      do iv = 1, CellbyProc(ip1)%Nb
         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)
            numCellRes = CellbyProc(ip1)%Num(k)

            NodebyCellRes_Ncpus(ip1)%Num(NodebyCellRes_Ncpus(ip1)%Pt(k) + 1:NodebyCellRes_Ncpus(ip1)%Pt(k + 1)) &
               = NodebyCell%Num(NodebyCell%Pt(numCellRes) + 1:NodebyCell%Pt(numCellRes + 1))
         enddo
      enddo

   end subroutine LocalMesh_NodebyCellRes

   ! Output:
   !  FacebyProc(ip)
   !  NbFaceOwn(ip), NbFaceRes(ip)
   ! Use:
   !  CellbyProc(ip1),FacebyCellRes_Ncpus(ip), ProcbyCell
   !> \brief Stores own+ghost faces of proc ip in order:
   !!   | own faces | ghost faces (own for proc i) | ghost faces (own for proc j>i) | ...
   subroutine LocalMesh_FacebyProc(ip, ProcbyCell)

      integer, intent(in) :: ip
      integer, dimension(:), intent(in) :: ProcbyCell
      integer :: ip1 ! ip1 = ip + 1 ip=0,1,..., ip1 used for array

      ! tmp
      integer :: nbProcVoisinFace
      integer, allocatable, dimension(:) :: &
         colorProc, cptFaceProcVois, colorFace
      integer :: i, ind, ipf, iv, j, kf, kipf, procf, procCV, k

      ip1 = ip + 1
      allocate (colorProc(Ncpus))
      allocate (colorFace(NbFace))

      ! all faces of cells own+ghost of proc i are colored
      colorFace(:) = -1
      do iv = 1, CellbyProc(ip1)%Nb   !
         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)  ! cell (own or ghost) of proc i
            colorFace(FacebyCellRes_Ncpus(ip1)%Num(FacebyCellRes_Ncpus(ip1)%Pt(k) + 1: &
                                                   FacebyCellRes_Ncpus(ip1)%Pt(k + 1))) = ip
         enddo
      enddo

      ! Counting nbProcVoisinFace
      nbProcVoisinFace = 0
      colorProc(:) = -1
      do i = 1, NbFace
         ! concerned faces
         if (colorFace(i) == ip) then
            ipf = -1; kipf = -1
            ! all cells connected to this face
            do j = CellbyFace%Pt(i) + 1, CellbyFace%Pt(i + 1)
               kf = CellbyFace%Num(j)
               ! proc that cell is own
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo
            colorProc(procf + 1) = 1
         endif
      enddo

      do i = 1, Ncpus
         if (colorProc(i) == 1) nbProcVoisinFace = nbProcVoisinFace + 1
      enddo

      FacebyProc(ip1)%Nb = nbProcVoisinFace
      allocate (FacebyProc(ip1)%Pt(nbProcVoisinFace + 1))
      allocate (FacebyProc(ip1)%Val(nbProcVoisinFace))

      ! %Val contains the numero of proc that face is own
      nbProcVoisinFace = 1
      FacebyProc(ip1)%Val(1) = ip
      do i = 0, Ncpus - 1
         if (colorProc(i + 1) == 1 .and. ip /= i) then
            nbProcVoisinFace = nbProcVoisinFace + 1
            FacebyProc(ip1)%Val(nbProcVoisinFace) = i
         endif
      enddo

      ! Filling of %Pt
      FacebyProc(ip1)%Pt(:) = 0
      do i = 1, NbFace
         ! concerned faces
         if (colorFace(i) == ip) then
            ipf = -1; kipf = -1
            ! all cells connected to this face
            do j = CellbyFace%Pt(i) + 1, CellbyFace%Pt(i + 1)
               kf = CellbyFace%Num(j)
               ! proc that cell is own
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo
            ! looking for the proc among the neighbour proc of ip
            ind = -1
            do iv = 1, nbProcVoisinFace
               if (FacebyProc(ip1)%Val(iv) == procf) then
                  ind = 1
                  ipf = iv
               endif
            enddo
            if (ind == -1) then
               print *, ' the face has not been found among the neighbour procs of ip ', ip, kf, procf
            endif
            ! face i is own for proc ipf
            FacebyProc(ip1)%Pt(ipf + 1:nbProcVoisinFace + 1) = FacebyProc(ip1)%Pt(ipf + 1:nbProcVoisinFace + 1) + 1
         endif
      enddo

      ! Counting and filling of %Pt
      allocate (cptFaceProcVois(nbProcVoisinFace + 1))
      cptFaceProcVois(:) = 0
      allocate (FacebyProc(ip1)%Num(FacebyProc(ip1)%Pt(nbProcVoisinFace + 1)))
      do i = 1, NbFace
         if (colorFace(i) == ip) then
            ipf = -1; kipf = -1
            ! all cells connected to this face
            do j = CellbyFace%Pt(i) + 1, CellbyFace%Pt(i + 1)
               kf = CellbyFace%Num(j)
               ! proc that cell is own
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo
            ! looking for the proc among the neighbour proc of ip
            ind = -1
            do iv = 1, nbProcVoisinFace
               if (FacebyProc(ip1)%Val(iv) == procf) then
                  ind = 1
                  ipf = iv
               endif
            enddo
            if (ind == -1) then
               print *, ' the face has not been found among the neighbour procs of ip ', ip, kf, procf
            endif
            cptFaceProcVois(ipf) = cptFaceProcVois(ipf) + 1
            FacebyProc(ip1)%Num(FacebyProc(ip1)%Pt(ipf) + cptFaceProcVois(ipf)) = i
         endif
      enddo

      deallocate (cptFaceProcVois)
      deallocate (colorFace)
      deallocate (colorProc)

      NbFaceResS_Ncpus(ip1) = FacebyProc(ip1)%Pt(FacebyProc(ip1)%Nb + 1) ! NbFaceResS_Ncpus
      NbFaceOwnS_Ncpus(ip1) = FacebyProc(ip1)%Pt(2) - FacebyProc(ip1)%Pt(1) ! NbNodeResS_Ncpus

   end subroutine LocalMesh_FacebyProc

   ! Output:
   !  NodebyFaceRes_Ncpus(ip)
   ! Use:
   !  NodebyFace, FacebyProc
   subroutine LocalMesh_NodebyFaceRes(ip)

      integer, intent(in) :: ip
      integer :: ip1

      ! tmp
      integer :: iv, k, Vsize, numFaceRes

      ip1 = ip + 1

      NodebyFaceRes_Ncpus(ip1)%Nb = NbFaceResS_Ncpus(ip1)
      allocate (NodebyFaceRes_Ncpus(ip1)%Pt(NbFaceResS_Ncpus(ip1) + 1))

      NodebyFaceRes_Ncpus(ip1)%Pt(1) = 0

      ! Counting
      do iv = 1, FacebyProc(ip1)%Nb
         do k = FacebyProc(ip1)%Pt(iv) + 1, FacebyProc(ip1)%Pt(iv + 1)
            numFaceRes = FacebyProc(ip1)%Num(k)

            NodebyFaceRes_Ncpus(ip1)%Pt(k + 1) = NodebyFaceRes_Ncpus(ip1)%Pt(k) &
                                                 + NodebyFace%Pt(numFaceRes + 1) - NodebyFace%Pt(numFaceRes)
         enddo
      enddo

      ! Filling
      Vsize = NodebyFaceRes_Ncpus(ip1)%Pt(NbFaceResS_Ncpus(ip1) + 1)
      allocate (NodebyFaceRes_Ncpus(ip1)%Num(Vsize))

      do iv = 1, FacebyProc(ip1)%Nb
         do k = FacebyProc(ip1)%Pt(iv) + 1, FacebyProc(ip1)%Pt(iv + 1)
            numFaceRes = FacebyProc(ip1)%Num(k)

            NodebyFaceRes_Ncpus(ip1)%Num(NodebyFaceRes_Ncpus(ip1)%Pt(k) + 1:NodebyFaceRes_Ncpus(ip1)%Pt(k + 1)) &
               = NodebyFace%Num(NodebyFace%Pt(numFaceRes) + 1:NodebyFace%Pt(numFaceRes + 1))
         enddo
      enddo

   end subroutine LocalMesh_NodebyFaceRes

   ! Output:
   !  CellbyFaceRes_Ncpus(ip)
   ! Use:
   !  CellbyFace, FacebyProc
   subroutine LocalMesh_CellbyFaceRes(ip)

      integer, intent(in) :: ip
      integer :: ip1

      ! tmp
      integer :: iv, k, Vsize, numFaceRes

      ip1 = ip + 1

      CellbyFaceRes_Ncpus(ip1)%Nb = NbFaceResS_Ncpus(ip1)
      allocate (CellbyFaceRes_Ncpus(ip1)%Pt(NbFaceResS_Ncpus(ip1) + 1))

      CellbyFaceRes_Ncpus(ip1)%Pt(1) = 0

      ! Counting
      do iv = 1, FacebyProc(ip1)%Nb
         do k = FacebyProc(ip1)%Pt(iv) + 1, FacebyProc(ip1)%Pt(iv + 1)
            numFaceRes = FacebyProc(ip1)%Num(k)

            CellbyFaceRes_Ncpus(ip1)%Pt(k + 1) = CellbyFaceRes_Ncpus(ip1)%Pt(k) &
                                                 + CellbyFace%Pt(numFaceRes + 1) - CellbyFace%Pt(numFaceRes)
         enddo
      enddo

      ! Filling
      Vsize = CellbyFaceRes_Ncpus(ip1)%Pt(NbFaceResS_Ncpus(ip1) + 1)
      allocate (CellbyFaceRes_Ncpus(ip1)%Num(Vsize))

      do iv = 1, FacebyProc(ip1)%Nb
         do k = FacebyProc(ip1)%Pt(iv) + 1, FacebyProc(ip1)%Pt(iv + 1)
            numFaceRes = FacebyProc(ip1)%Num(k)

            CellbyFaceRes_Ncpus(ip1)%Num(CellbyFaceRes_Ncpus(ip1)%Pt(k) + 1:CellbyFaceRes_Ncpus(ip1)%Pt(k + 1)) &
               = CellbyFace%Num(CellbyFace%Pt(numFaceRes) + 1:CellbyFace%Pt(numFaceRes + 1))
         enddo
      enddo

   end subroutine LocalMesh_CellbyFaceRes

   pure subroutine LocalMesh_NodebyProc_add_wells(well_nodes, color, node_color, nb_wells, new_color)
      type(CSR), intent(in) :: well_nodes
      integer, intent(in) :: color
      integer, dimension(:), intent(in) :: node_color
      integer, intent(out) :: nb_wells
      integer, dimension(:), intent(inout) :: new_color

      integer :: wk, i, j

      nb_wells = 0
      do wk = 1, well_nodes%Nb
         ! if color is found on one node apply it to every node of the well
         do i = well_nodes%Pt(wk) + 1, well_nodes%Pt(wk + 1)
            if (node_color(well_nodes%Num(i)) == color) then
               nb_wells = nb_wells + 1
               do j = well_nodes%Pt(wk) + 1, well_nodes%Pt(wk + 1)
                  new_color(well_nodes%Num(j)) = color
               end do
               exit
            end if
         end do
      end do

   end subroutine LocalMesh_NodebyProc_add_wells

   ! Output:
   !   NodebyProc
   !   CSR with  | .... own .... | ghost from proc i | ghost from proc j (i<j) |  ... |
   ! Use:
   !   NodebyCellRes, ProcbyCell
   !> \brief Build the CSR with  | .... own .... | ghost from proc i | ghost from proc j (i<j) |  ... |
   !! The nodes Own+Ghost are the nodes of the cells Own+Ghost (contained in CellbyProc)
   !!                                       PLUS the ghost cells concerning the wells             <br>
   !! if one node of the well is contained by the proc (node own or ghost) then
   !! every node of the well must be own or ghost for this proc
   subroutine LocalMesh_NodebyProc(ip, ProcbyCell)

      integer, intent(in) :: ip
      integer, dimension(:), intent(in) :: ProcbyCell
      integer :: ip1

      ! tmp
      integer :: nbProcVoisinNode, &
                 i, ipf, kipf, ind, iv, j, kf, procCV, procf, k
      integer, allocatable, dimension(:) :: &
         colorProc, colorNode, colorNodeTmp, cptNodeProcVois

      ip1 = ip + 1
      allocate (colorProc(Ncpus))
      allocate (colorNodeTmp(NbNode))
      allocate (colorNode(NbNode))

      ! colorNodeTmp, node own+ghost of this proc (WITHOUT node ghost due to the wells)
      colorNodeTmp(:) = -1
      ! loop over cells own+ghost of this proc
      do iv = 1, CellbyProc(ip1)%Nb
         do k = CellbyProc(ip1)%Pt(iv) + 1, CellbyProc(ip1)%Pt(iv + 1)
            colorNodeTmp(NodebyCellRes_Ncpus(ip1)%Num(NodebyCellRes_Ncpus(ip1)%Pt(k) + 1: &
                                                      NodebyCellRes_Ncpus(ip1)%Pt(k + 1))) = ip
         enddo
      enddo

      ! colorNode, node own+ghost of this proc (INCLUDING node ghost due to the wells)
      colorNode(:) = colorNodeTmp(:)

      ! Wells
      ! if one node of the well is contained by the proc (node own or ghost) then
      ! every nodes of the well must be own or ghost for this proc
      call LocalMesh_NodebyProc_add_wells(NodebyWellInj, ip, colorNodeTmp, NbWellInjResS_Ncpus(ip1), colorNode)
      call LocalMesh_NodebyProc_add_wells(NodebyWellProd, ip, colorNodeTmp, NbWellProdResS_Ncpus(ip1), colorNode)

      ! Counting nbProcVoisinNode (nbr of proc from which ghost are synchronised)
      nbProcVoisinNode = 0
      ! colorProc(IdProc+1) = 1 if IdProc is neighbour of the concerned proc
      colorProc(:) = -1
      ! we include current proc because in very rare configuration it may have no own nodes
      colorProc(ip1) = 1

      do i = 1, NbNode
         if (colorNode(i) == ip) then
! the following loop is a criterium useful if node i is at the interface between 2 proc
! It distinguishes which proc will have this node as Own
! Maybe equivalent to:
!          j=CellbyNode%Pt(i+1)
!          kf = CellbyNode%Num(j)
!          ! proc containing this cell
!          procCV = ProcbyCell(kf)
!          colorProc(procCV+1) = 1

            ipf = -1; kipf = -1
            ! loop over the cells containing this node
            do j = CellbyNode%Pt(i) + 1, CellbyNode%Pt(i + 1)
               kf = CellbyNode%Num(j)

               ! proc containing this cell
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo

            colorProc(procf + 1) = 1
         endif
      enddo

      do i = 1, Ncpus
         if (colorProc(i) == 1) then
            nbProcVoisinNode = nbProcVoisinNode + 1
         end if
      enddo

      NodebyProc(ip1)%Nb = nbProcVoisinNode
      allocate (NodebyProc(ip1)%Pt(nbProcVoisinNode + 1))
      allocate (NodebyProc(ip1)%Val(nbProcVoisinNode))

      NodebyProc(ip1)%Val(1) = ip
      nbProcVoisinNode = 1
      do i = 0, Ncpus - 1
         if (colorProc(i + 1) == 1 .and. i /= ip) then
            nbProcVoisinNode = nbProcVoisinNode + 1
            NodebyProc(ip1)%Val(nbProcVoisinNode) = i
         endif
      enddo

      ! filling of Pt
      NodebyProc(ip1)%Pt(:) = 0
      do i = 1, NbNode
         ! concerned nodes
         if (colorNode(i) == ip) then
            ipf = -1; kipf = -1
            ! all cells connected to this node
            do j = CellbyNode%Pt(i) + 1, CellbyNode%Pt(i + 1)
               kf = CellbyNode%Num(j)
               ! proc that cell is own
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo
            ! looking for the proc among the neighbour proc of ip
            ind = -1
            do iv = 1, nbProcVoisinNode
               if (NodebyProc(ip1)%Val(iv) == procf) then
                  ind = 1
                  ipf = iv
               endif
            enddo
            if (ind == -1) then
               print *, ' the node has not been found among the neighbour procs of ip ', ip, kf, procf
            endif
            ! node i is own for proc ipf
            NodebyProc(ip1)%Pt(ipf + 1:nbProcVoisinNode + 1) = NodebyProc(ip1)%Pt(ipf + 1:nbProcVoisinNode + 1) + 1
         endif
      enddo

      ! Counting and Filling of %Pt
      allocate (cptNodeProcVois(nbProcVoisinNode + 1))
      cptNodeProcVois(:) = 0

      allocate (NodebyProc(ip1)%Num(NodebyProc(ip1)%Pt(nbProcVoisinNode + 1)))
      do i = 1, NbNode
         if (colorNode(i) == ip) then
            ipf = -1; kipf = -1

            ! all cells connected to this node
            do j = CellbyNode%Pt(i) + 1, CellbyNode%Pt(i + 1)
               kf = CellbyNode%Num(j)
               ! proc that cell is own
               procCV = ProcbyCell(kf)
               if ((ipf == -1) .or. (kf < kipf)) procf = procCV; kipf = kf
            enddo

            ! looking for the proc among the neighbour proc of ip
            ind = -1
            do iv = 1, nbProcVoisinNode
               if (NodebyProc(ip1)%Val(iv) == procf) then
                  ind = 1
                  ipf = iv
               endif
            enddo
            if (ind == -1) then
               print *, ' the node has not been found among the neighbour procs of ip ', ip, kf, procf
            endif
            ! node i is own for proc ipf
            cptNodeProcVois(ipf) = cptNodeProcVois(ipf) + 1
            NodebyProc(ip1)%Num(NodebyProc(ip1)%Pt(ipf) + cptNodeProcVois(ipf)) = i
         endif
      enddo

      ! write(*,*)'ip1',ip1

      ! do i=1, NodebyProc(ip1)%Nb
      !    write(*,*) getpid(), NodebyProc(ip1)%Num(NodebyProc(ip1)%Pt(i)+1:NodebyProc(ip1)%Pt(i+1))
      ! end do

      deallocate (cptNodeProcVois)
      deallocate (colorProc)
      deallocate (colorNode)
      deallocate (colorNodeTmp)

      NbNodeResS_Ncpus(ip1) = NodebyProc(ip1)%Pt(nbProcVoisinNode + 1)    ! NbNodeResS_Ncpus
      NbNodeOwnS_Ncpus(ip1) = NodebyProc(ip1)%Pt(2) - NodebyProc(ip1)%Pt(1)    ! NbNodeOwnS_Ncpus

   end subroutine LocalMesh_NodebyProc

#ifdef NDEBUG
   pure &
#endif
      subroutine LocalMesh_count_wells_in_part(well_nodes, node_proc, well_proc, proc_nb_wells)
      type(CSR), intent(in) :: well_nodes
      integer, dimension(:), intent(in) :: node_proc
      integer, dimension(:), intent(out) :: well_proc
      integer, dimension(:), intent(out) :: proc_nb_wells

      integer :: i, wk, pid, nb_wells, head_node

      nb_wells = well_nodes%Nb

#ifndef NDEBUG
      if (size(well_proc) < nb_wells) &
         call CommonMPI_abort("Inconsistent is_well_seen work array.")
      if (size(proc_nb_wells) /= Ncpus) &
         call CommonMPI_abort("Inconsistent proc_nb_wells work array.")
#endif

      well_proc = -1
      proc_nb_wells = 0
      do wk = 1, nb_wells
         head_node = well_nodes%Num(well_nodes%Pt(wk + 1))
         pid = node_proc(head_node)
         if (pid >= 0) then ! proc sees the well
#ifndef NDEBUG
            ! Check that all well nodes are in part
            do i = well_nodes%Pt(wk) + 1, well_nodes%Pt(wk + 1)
               if (node_proc(well_nodes%Num(i)) < 0) &
                  call CommonMPI_abort("All well nodes should be in part!")
            end do
#endif
            well_proc(wk) = pid
            proc_nb_wells(pid + 1) = proc_nb_wells(pid + 1) + 1
#ifndef NDEBUG
         else ! Check that all well nodes are out of part
            do i = well_nodes%Pt(wk) + 1, well_nodes%Pt(wk + 1)
               if (node_proc(well_nodes%Num(i)) >= 0) &
                  call CommonMPI_abort("All well nodes should be out of part!")
            end do
#endif
         end if
      end do

   end subroutine LocalMesh_count_wells_in_part

#ifdef NDEBUG
   pure &
#endif
      subroutine LocalMesh_WellbyProc_part(well_nodes, procs_in_part, node_proc, well_proc, proc_nb_wells, WellbyProc)
      type(CSR), intent(in) :: well_nodes ! all wells
      integer, dimension(:), intent(in) :: procs_in_part
      integer, dimension(:), intent(in) :: node_proc ! -1 if node is not in part, else id of proc that owns the node
      integer, dimension(:), intent(out) :: well_proc
      integer, dimension(:), intent(out) :: proc_nb_wells ! work array
      type(CSR), intent(out) :: WellbyProc

      integer :: i, wk, pid, nb_wells, nb_procs_in_part

      nb_procs_in_part = size(procs_in_part)

      call LocalMesh_count_wells_in_part(well_nodes, node_proc, well_proc, proc_nb_wells)

#ifndef NDEBUG
      if (nb_procs_in_part < count(proc_nb_wells > 0)) &
         call CommonMPI_abort("Not enough proc!")
      if (sum(proc_nb_wells) /= count(well_proc >= 0)) &
         call CommonMPI_abort("LocalMesh_count_wells_in_part returned inconsistent result.")
#endif

      WellbyProc%Nb = nb_procs_in_part

      ! copy the pid of part procs
      ! FIXME: this data should be shared and NOT copied
      allocate (WellbyProc%Val(nb_procs_in_part))
      WellbyProc%Val = procs_in_part

      allocate (WellbyProc%Pt(nb_procs_in_part + 1))
      WellbyProc%Pt(:) = 0
      do i = 1, nb_procs_in_part
         pid = procs_in_part(i)
         WellbyProc%Pt(i + 1) = WellbyProc%Pt(i) + proc_nb_wells(pid + 1)
      enddo

      allocate (WellbyProc%Num(sum(proc_nb_wells)))
      proc_nb_wells = 0
      nb_wells = well_nodes%Nb
      do wk = 1, nb_wells
         pid = well_proc(wk)
         if (pid < 0) cycle
         ! find local id of proc
         do i = 1, nb_procs_in_part
            if (procs_in_part(i) == pid) then
               proc_nb_wells(pid + 1) = proc_nb_wells(pid + 1) + 1
#ifndef NDEBUG
               if (WellbyProc%Pt(i) + proc_nb_wells(pid + 1) > WellbyProc%Pt(i + 1)) &
                  call CommonMPI_abort("Well proc inconsistency.")
#endif
               ! FIXME: we could use the well id here
               !        which could be a way to have all wells stored in the same structure
               WellbyProc%Num(WellbyProc%Pt(i) + proc_nb_wells(pid + 1)) = wk
               exit
            end if
         end do
      end do

   end subroutine LocalMesh_WellbyProc_part

   pure subroutine LocalMesh_WellbyProc_copy_welldata(part_wells, all_data, part_data)
      integer, dimension(:), intent(in) :: part_wells
      type(WellData_type), dimension(:), intent(in) :: all_data
      type(WellData_type), dimension(:), intent(out) :: part_data

      integer :: i

      do i = 1, size(part_wells)
         part_data(i) = all_data(part_wells(i))
      enddo

   end subroutine LocalMesh_WellbyProc_copy_welldata

#ifdef NDEBUG
   pure &
#endif
      subroutine LocalMesh_color_global_nodes(selection, color)
      type(CSR), intent(in) :: selection
      integer, dimension(:), intent(out) :: color

      integer :: i, j

#ifndef NDEBUG
      if (size(color) /= NbNode) &
         call CommonMPI_abort("Inconsistent color work array.")
#endif

      color(:) = -1
      do i = 1, selection%Nb
         do j = selection%Pt(i) + 1, selection%Pt(i + 1)
            color(selection%Num(j)) = selection%Val(i)
         enddo
      enddo

   end subroutine LocalMesh_color_global_nodes

   ! Output:
   !  NbWellOwnS_Ncpus
   !  WellbyProc
   !  DataWellRes_Ncpus
   ! Use:
   !  NodebyProc, CellbyProc, NodebyWell, DataWell
   !> \brief Stores the global number of the wells (own or ghost) of proc ip in a CSR:                   <br>
   !!   The number of rows is the same as that of NodebyProc.                  <br>
   !!   If there is no well in proc i, %Pt(i+1)=%Pt(i).                  <br>
   !!  | wells own | wells ghost (which are own for proc i) | wells ghost (which are own for proc j>i) | ...                  <br>
   !! A well is own if the head Node of the well is own                  <br>
   !! Carreful if (at least) two wells share a node
   subroutine LocalMesh_WellbyProc(ip, node_proc, proc_nb_wells, well_proc)
      integer, intent(in) :: ip
      integer, dimension(:), intent(out) :: node_proc, proc_nb_wells, well_proc ! work arrays

      integer :: ip1

      ip1 = ip + 1

      call LocalMesh_color_global_nodes(NodebyProc(ip1), node_proc)

      call LocalMesh_WellbyProc_part(NodebyWellInj, NodebyProc(ip1)%Val, node_proc, well_proc, proc_nb_wells, WellInjbyProc(ip1))
      call LocalMesh_WellbyProc_copy_welldata(WellInjbyProc(ip1)%Num, DataWellInj, DataWellInjRes_Ncpus(:, ip1))
      NbWellInjOwnS_Ncpus(ip1) = WellInjbyProc(ip1)%Pt(2)

      call LocalMesh_WellbyProc_part(NodebyWellProd, NodebyProc(ip1)%Val, node_proc, well_proc, proc_nb_wells, WellProdbyProc(ip1))
      call LocalMesh_WellbyProc_copy_welldata(WellProdbyProc(ip1)%Num, DataWellProd, DataWellProdRes_Ncpus(:, ip1))
      NbWellProdOwnS_Ncpus(ip1) = WellProdbyProc(ip1)%Pt(2)

   end subroutine LocalMesh_WellbyProc

   subroutine LocalMesh_Flags(ip)
      integer, intent(in) :: ip
      integer :: k, n, ip1

      ip1 = ip + 1

      n = size(NodebyProc(ip1)%Num)
      allocate (NodeFlags_Ncpus(ip1)%Val(n))
      do k = 1, n
         NodeFlags_Ncpus(ip1)%Val(k) = NodeFlags(NodebyProc(ip1)%Num(k))
      end do

      n = size(CellbyProc(ip1)%Num)
      allocate (CellFlags_Ncpus(ip1)%Val(n))
      do k = 1, n
         CellFlags_Ncpus(ip1)%Val(k) = CellFlags(CellbyProc(ip1)%Num(k))
      end do

      n = size(FacebyProc(ip1)%Num)
      allocate (FaceFlags_Ncpus(ip1)%Val(n))
      do k = 1, n
         FaceFlags_Ncpus(ip1)%Val(k) = FaceFlags(FacebyProc(ip1)%Num(k))
      end do

      n = size(CellbyProc(ip1)%Num)
      allocate (CellTypes_Ncpus(ip1)%Val(n))
      do k = 1, n
         CellTypes_Ncpus(ip1)%Val(k) = CellTypes(CellbyProc(ip1)%Num(k))
      end do

      n = size(FacebyProc(ip1)%Num)
      allocate (FaceTypes_Ncpus(ip1)%Val(n))
      do k = 1, n
         FaceTypes_Ncpus(ip1)%Val(k) = FaceTypes(FacebyProc(ip1)%Num(k))
      end do

   end subroutine LocalMesh_Flags

   !> \brief Initialize Rocktypes of the local from the info of the global mesh
   subroutine LocalMesh_Rocktype(ip)
      integer, intent(in) :: ip
      integer :: k, n, ip1

      ip1 = ip + 1

      n = size(NodebyProc(ip1)%Num)
      allocate (NodeRocktype_Ncpus(ip1)%Array2d(IndThermique + 1, n))
      do k = 1, n
         NodeRocktype_Ncpus(ip1)%Array2d(:, k) = NodeRocktype(:, NodebyProc(ip1)%Num(k))
      end do

      n = size(CellbyProc(ip1)%Num)
      allocate (CellRocktype_Ncpus(ip1)%Array2d(IndThermique + 1, n))
      do k = 1, n
         CellRocktype_Ncpus(ip1)%Array2d(:, k) = CellRocktype(:, CellbyProc(ip1)%Num(k))
      end do

      n = size(FracbyProc(ip1)%Num)
      allocate (FracRocktype_Ncpus(ip1)%Array2d(IndThermique + 1, n))
      do k = 1, n
         FracRocktype_Ncpus(ip1)%Array2d(:, k) = FracRocktype(:, FracbyProc(ip1)%Num(k))
      end do

   end subroutine LocalMesh_Rocktype

   ! Output:
   !  XNodeRes
   ! Use:
   !  XNode, NodebyProc
   subroutine LocalMesh_XNodeRes(ip)

      integer, intent(in) :: ip
      integer :: ip1

      ! tmp
      integer :: cpt, iv, k, numNodeRes

      ip1 = ip + 1
      allocate (XNodeRes_Ncpus(ip1)%Array2d(3, NbNodeResS_Ncpus(ip1)))

      cpt = 0
      do iv = 1, NodebyProc(ip1)%Nb
         do k = NodebyProc(ip1)%Pt(iv) + 1, NodebyProc(ip1)%Pt(iv + 1)
            numNodeRes = NodebyProc(ip1)%Num(k)
            cpt = cpt + 1
            XNodeRes_Ncpus(ip1)%Array2d(:, cpt) = XNode(:, numNodeRes)
         enddo
      enddo

      if (cpt /= NbNodeResS_Ncpus(ip1)) then
         print *, 'pb : cpt = ', cpt, ' vs NbNodeResS_Ncpus = ', NbNodeResS_Ncpus(ip1)
      endif

   end subroutine LocalMesh_XNodeRes

   ! Output:
   !  NodebyWellRes_Ncpus(ip1)
   ! Use:
   !  WellbyProc, NodebyWell
   !> \brief Nodes of local wells (global number)                  <br>
   !! 1 local well = 1 line in NodebyWellRes_Ncpus (CSR) in order:                  <br>
   !!   | nodes of well own 1 | nodes of well own 2 | ... | nodes of well ghost 1 | nodes of well ghost 2 | ...
   subroutine LocalMesh_NodebyWellRes(ip)

      integer, intent(in) :: ip
      integer :: ip1, k, numWellRes, Vsize, cpt

      ip1 = ip + 1

      !! INJ WELL
      NodebyWellInjRes_Ncpus(ip1)%Nb = NbWellInjResS_Ncpus(ip1)
      allocate (NodebyWellInjRes_Ncpus(ip1)%Pt(NodebyWellInjRes_Ncpus(ip1)%Nb + 1))

      NodebyWellInjRes_Ncpus(ip1)%Pt(:) = 0

      ! filling of Pt
      do k = 1, WellInjbyProc(ip1)%Pt(WellInjbyProc(ip1)%Nb + 1)
         numWellRes = WellInjbyProc(ip1)%Num(k)
         ! number of nodes in well numWellRes
         NodebyWellInjRes_Ncpus(ip1)%Pt(k + 1) = NodebyWellInjRes_Ncpus(ip1)%Pt(k) &
                                                 + NodebyWellInj%Pt(numWellRes + 1) - NodebyWellInj%Pt(numWellRes)
      enddo

      Vsize = NodebyWellInjRes_Ncpus(ip1)%Pt(NodebyWellInjRes_Ncpus(ip1)%Nb + 1)
      allocate (NodebyWellInjRes_Ncpus(ip1)%Num(Vsize))

      ! filling of Num
      do k = 1, WellInjbyProc(ip1)%Pt(WellInjbyProc(ip1)%Nb + 1)

         numWellRes = WellInjbyProc(ip1)%Num(k)
         do cpt = 1, NodebyWellInj%Pt(numWellRes + 1) - NodebyWellInj%Pt(numWellRes)
            NodebyWellInjRes_Ncpus(ip1)%Num(NodebyWellInjRes_Ncpus(ip1)%Pt(k) + cpt) = &
               NodebyWellInj%Num(NodebyWellInj%Pt(numWellRes) + cpt)
         enddo
      enddo

      !! PROD WELL
      NodebyWellProdRes_Ncpus(ip1)%Nb = NbWellProdResS_Ncpus(ip1)
      allocate (NodebyWellProdRes_Ncpus(ip1)%Pt(NodebyWellProdRes_Ncpus(ip1)%Nb + 1))

      NodebyWellProdRes_Ncpus(ip1)%Pt(:) = 0

      ! filling of Pt
      do k = 1, WellProdbyProc(ip1)%Pt(WellProdbyProc(ip1)%Nb + 1)

         numWellRes = WellProdbyProc(ip1)%Num(k)
         ! number of nodes in well numWellRes
         NodebyWellProdRes_Ncpus(ip1)%Pt(k + 1) = NodebyWellProdRes_Ncpus(ip1)%Pt(k) &
                                                  + NodebyWellProd%Pt(numWellRes + 1) - NodebyWellProd%Pt(numWellRes)
      enddo

      Vsize = NodebyWellProdRes_Ncpus(ip1)%Pt(NodebyWellProdRes_Ncpus(ip1)%Nb + 1)
      allocate (NodebyWellProdRes_Ncpus(ip1)%Num(Vsize))

      ! filling of Num
      do k = 1, WellProdbyProc(ip1)%Pt(WellProdbyProc(ip1)%Nb + 1)

         numWellRes = WellProdbyProc(ip1)%Num(k)
         do cpt = 1, NodebyWellProd%Pt(numWellRes + 1) - NodebyWellProd%Pt(numWellRes)
            NodebyWellProdRes_Ncpus(ip1)%Num(NodebyWellProdRes_Ncpus(ip1)%Pt(k) + cpt) = &
               NodebyWellProd%Num(NodebyWellProd%Pt(numWellRes) + cpt)
         enddo
      enddo

   end subroutine LocalMesh_NodebyWellRes

   ! Output:
   !  localbyGlobalCell/Face/Node
   ! Use:
   !  CellbyProc, FacebyProc, NodebyProc
   ! Size NbCell/NbFace/NbNode
   !> \brief   Contains the local number of Cell/Face/Node given the global number for proc ip
   subroutine LocalMesh_GlobalToLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, i, j

      ip1 = ip + 1

      ! initialization to zero
      localbyGlobalCell(:) = 0
      localbyGlobalFace(:) = 0
      localbyGlobalNode(:) = 0

      ! global to local
      do i = 1, CellbyProc(ip1)%Nb
         do j = CellbyProc(ip1)%Pt(i) + 1, CellbyProc(ip1)%Pt(i + 1)
            localbyGlobalCell(CellbyProc(ip1)%Num(j)) = j
         enddo
      enddo

      do i = 1, FacebyProc(ip1)%Nb
         do j = FacebyProc(ip1)%Pt(i) + 1, FacebyProc(ip1)%Pt(i + 1)
            localbyGlobalFace(FacebyProc(ip1)%Num(j)) = j
         enddo
      enddo

      do i = 1, NodebyProc(ip1)%Nb
         do j = NodebyProc(ip1)%Pt(i) + 1, NodebyProc(ip1)%Pt(i + 1)
            localbyGlobalNode(NodebyProc(ip1)%Num(j)) = j
         enddo
      enddo

   end subroutine LocalMesh_GlobalToLocal

   ! Output:
   !  NodebyWellLocal_Ncpus(ip))
   ! Use:
   !  NodebyWellRes_Ncpus(ip), localbyGlobalNode
   !> \brief Nodes (local number) of wells (own+ghost) of proc ip
   subroutine LocalMesh_NodebyWellLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, Nnnz, i

      ip1 = ip + 1

      ! INJ WELLS
      Nb = NodebyWellInjRes_Ncpus(ip1)%Nb
      NodebyWellInjLocal_Ncpus(ip1)%Nb = Nb
      Nnnz = NodebyWellInjRes_Ncpus(ip1)%Pt(Nb + 1)
      allocate (NodebyWellInjLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodebyWellInjLocal_Ncpus(ip1)%Num(Nnnz))

      do i = 1, Nb + 1
         NodebyWellInjLocal_Ncpus(ip1)%Pt(i) = NodebyWellInjRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         NodebyWellInjLocal_Ncpus(ip1)%Num(i) = localbyGlobalNode(NodebyWellInjRes_Ncpus(ip1)%Num(i))
      enddo

      ! PROD WELLS
      Nb = NodebyWellProdRes_Ncpus(ip1)%Nb
      NodebyWellProdLocal_Ncpus(ip1)%Nb = Nb
      Nnnz = NodebyWellProdRes_Ncpus(ip1)%Pt(Nb + 1)
      allocate (NodebyWellProdLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodebyWellProdLocal_Ncpus(ip1)%Num(Nnnz))

      do i = 1, Nb + 1
         NodebyWellProdLocal_Ncpus(ip1)%Pt(i) = NodebyWellProdRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         NodebyWellProdLocal_Ncpus(ip1)%Num(i) = localbyGlobalNode(NodebyWellProdRes_Ncpus(ip1)%Num(i))
      enddo

   end subroutine LocalMesh_NodebyWellLocal

   ! Output:
   !  FacebyCellLocal_Ncpus(ip)
   ! Use:
   !  FacebyCellRes_Ncpus(ip), localbyGlobalCell
   !> \brief Faces (local number) of cells (own+ghost) of proc ip
   subroutine LocalMesh_FacebyCellLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, Nnnz, i

      ip1 = ip + 1

      ! FacebyCellLocal_Ncpus
      Nb = FacebyCellRes_Ncpus(ip1)%Nb
      Nnnz = FacebyCellRes_Ncpus(ip1)%Pt(Nb + 1)
      allocate (FacebyCellLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (FacebyCellLocal_Ncpus(ip1)%Num(Nnnz))

      FacebyCellLocal_Ncpus(ip1)%Nb = Nb
      do i = 1, Nb + 1
         FacebyCellLocal_Ncpus(ip1)%Pt(i) = FacebyCellRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         FacebyCellLocal_Ncpus(ip1)%Num(i) = localbyGlobalFace(FacebyCellRes_Ncpus(ip1)%Num(i))
      enddo

   end subroutine LocalMesh_FacebyCellLocal

   ! Output:
   !  NodebyCellLocal_Ncpus(ip)
   ! Use:
   !  NodebyCellRes_Ncpus(ip), localbyGlobalNode
   !> \brief Nodes (local number) of cells (own+ghost) of proc ip
   subroutine LocalMesh_NodebyCellLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, Nnnz, i

      ip1 = ip + 1

      Nb = NodebyCellRes_Ncpus(ip1)%Nb
      Nnnz = NodebyCellRes_Ncpus(ip1)%Pt(Nb + 1)
      allocate (NodebyCellLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodebyCellLocal_Ncpus(ip1)%Num(Nnnz))

      NodebyCellLocal_Ncpus(ip1)%Nb = Nb
      do i = 1, Nb + 1
         NodebyCellLocal_Ncpus(ip1)%Pt(i) = NodebyCellRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         NodebyCellLocal_Ncpus(ip1)%Num(i) = localbyGlobalNode(NodebyCellRes_Ncpus(ip1)%Num(i))
      enddo

   end subroutine LocalMesh_NodebyCellLocal

   ! Output:
   !  NodebyFaceLocal_Ncpus(ip)
   ! Use:
   !  NodebyFaceRes_Ncpus(ip), localbyGlobalNode
   !> \brief Nodes (local number) of faces (own+ghost) of proc ip
   subroutine LocalMesh_NodebyFaceLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, Nnnz, i

      ip1 = ip + 1

      Nb = NodebyFaceRes_Ncpus(ip1)%Nb
      Nnnz = NodebyFaceRes_Ncpus(ip1)%Pt(Nb + 1)
      allocate (NodebyFaceLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodebyFaceLocal_Ncpus(ip1)%Num(Nnnz))

      NodebyFaceLocal_Ncpus(ip1)%Nb = Nb
      do i = 1, Nb + 1
         NodebyFaceLocal_Ncpus(ip1)%Pt(i) = NodebyFaceRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         NodebyFaceLocal_Ncpus(ip1)%Num(i) = localbyGlobalNode(NodebyFaceRes_Ncpus(ip1)%Num(i))
      enddo

   end subroutine LocalMesh_NodebyFaceLocal

   ! Output:
   !  CellbyFaceLocal_Ncpus(ip)
   ! Use:
   !  CellbyFaceRes_Ncpus(ip)
   subroutine LocalMesh_CellbyFaceLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, Nnnz, i

      ip1 = ip + 1

      Nb = CellbyFaceRes_Ncpus(ip1)%Nb
      Nnnz = CellbyFaceRes_Ncpus(ip1)%Pt(Nb + 1)

      allocate (CellbyFaceLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (CellbyFaceLocal_Ncpus(ip1)%Num(Nnnz))

      CellbyFaceLocal_Ncpus(ip1)%Nb = Nb
      do i = 1, Nb + 1
         CellbyFaceLocal_Ncpus(ip1)%Pt(i) = CellbyFaceRes_Ncpus(ip1)%Pt(i)
      end do
      do i = 1, Nnnz
         CellbyFaceLocal_Ncpus(ip1)%Num(i) = localbyGlobalCell(CellbyFaceRes_Ncpus(ip1)%Num(i))
      enddo

   end subroutine LocalMesh_CellbyFaceLocal

   ! Output:
   !  IdCellRes_Ncpus(ip)
   ! Use:
   !  IdCell, CellbyProc
   !> \brief Id Cell of cells own+ghost of proc ip stored in %Val
   subroutine LocalMesh_IdCellRes(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1)
      allocate (IdCellRes_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         IdCellRes_Ncpus(ip1)%Val(i) = IdCell(CellbyProc(ip1)%Num(i))
      end do

   end subroutine LocalMesh_IdCellRes

   ! Output:
   !  IdFaceRes_Ncpus(ip)
   ! Use:
   !  IdFace, FacebyProc
   !> \brief Id Face of faces own+ghost of proc ip stored in %Val
   subroutine LocalMesh_IdFaceRes(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = FacebyProc(ip1)%Pt(FacebyProc(ip1)%Nb + 1)
      allocate (IdFaceRes_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         IdFaceRes_Ncpus(ip1)%Val(i) = IdFace(FacebyProc(ip1)%Num(i))
      end do

   end subroutine LocalMesh_IdFaceRes

   ! Output:
   !  IdNodeRes_Ncpus(ip)
   ! Use:
   !  IdNode, NodebyProc
   !> \brief Id Node of nodes own+ghost of proc ip stored in %Val
   subroutine LocalMesh_IdNodeRes(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = NbNodeResS_Ncpus(ip1)
      allocate (IdNodeRes_Ncpus(ip1)%Val(Nb))

      do i = 1, NbNodeOwnS_Ncpus(ip1)
         IdNodeRes_Ncpus(ip1)%Val(i)%proc = "o" ! own node
         IdNodeRes_Ncpus(ip1)%Val(i)%Frac = IdNode(NodebyProc(ip1)%Num(i))%Frac
         IdNodeRes_Ncpus(ip1)%Val(i)%P = IdNode(NodebyProc(ip1)%Num(i))%P
         IdNodeRes_Ncpus(ip1)%Val(i)%T = IdNode(NodebyProc(ip1)%Num(i))%T
      end do

      do i = NbNodeOwnS_Ncpus(ip1) + 1, NbNodeResS_Ncpus(ip1)
         IdNodeRes_Ncpus(ip1)%Val(i)%proc = "g" ! ghost node
         IdNodeRes_Ncpus(ip1)%Val(i)%Frac = IdNode(NodebyProc(ip1)%Num(i))%Frac
         IdNodeRes_Ncpus(ip1)%Val(i)%P = IdNode(NodebyProc(ip1)%Num(i))%P
         IdNodeRes_Ncpus(ip1)%Val(i)%T = IdNode(NodebyProc(ip1)%Num(i))%T
      end do

   end subroutine LocalMesh_IdNodeRes

   ! Output:
   !  PorositeCell_Ncpus(ip), PorositeFrac_Ncpus(ip)
   ! Use:
   !  CellbyProc(ip), PorositeCell,
   !  FracbyProc(ip), PorositeFrac
   subroutine LocalMesh_Porosite(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1)
      allocate (PorositeCell_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         PorositeCell_Ncpus(ip1)%Val(i) = PorositeCell(CellbyProc(ip1)%Num(i))
      end do

      Nb = FracbyProc(ip1)%Pt(FracbyProc(ip1)%Nb + 1)
      !print *, "DEBUG - Copying", Nb, "fracture/face global porosities from:", PorositeFrac
      allocate (PorositeFrac_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         !print *, "DEBUG - Allocating fracture", i, &
         !         "which is face", FracbyProc(ip1)%Num(i), &
         !         !"with porosity value", PorositeFrac( FracbyProc(ip1)%Num(i))
         !         "with porosity value", PorositeFrac(i)
         !PorositeFrac_Ncpus(ip1)%Val(i) = PorositeFrac(FracbyProc(ip1)%Num(i))
         PorositeFrac_Ncpus(ip1)%Val(i) = PorositeFrac(FracbyProc(ip1)%Num(i))
      end do

   end subroutine LocalMesh_Porosite

   ! Output:
   !  PermCellLocal_Ncpus(ip), PermFracLocal_Ncpus(ip)
   ! Use:
   !  CellbyProc(ip), PermCell,
   !  FracbyProc(ip), PermFace
   subroutine LocalMesh_Perm(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1)
      allocate (PermCellLocal_Ncpus(ip1)%Array3d(3, 3, Nb))
      do i = 1, Nb
         PermCellLocal_Ncpus(ip1)%Array3d(:, :, i) = PermCell(:, :, CellbyProc(ip1)%Num(i))
      end do

      Nb = FracbyProc(ip1)%Pt(FracbyProc(ip1)%Nb + 1)
      allocate (PermFracLocal_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         PermFracLocal_Ncpus(ip1)%Val(i) = PermFrac(FracbyProc(ip1)%Num(i))
      end do

   end subroutine LocalMesh_Perm

#ifdef _THERMIQUE_

   ! Output:
   !  CondThermalCellLocal_Ncpus(ip), CondThermalFracLocal_Ncpus(ip)
   ! Use:
   !  CellbyProc(ip), CondThermalCell,
   !  FracbyProc(ip), CondThermalFace
   subroutine LocalMesh_CondThermal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1)
      allocate (CondThermalCellLocal_Ncpus(ip1)%Array3d(3, 3, Nb))
      do i = 1, Nb
         CondThermalCellLocal_Ncpus(ip1)%Array3d(:, :, i) = CondThermalCell(:, :, CellbyProc(ip1)%Num(i))
      end do

      Nb = FracbyProc(ip1)%Pt(FracbyProc(ip1)%Nb + 1)
      allocate (CondThermalFracLocal_Ncpus(ip1)%Val(Nb))
      do i = 1, Nb
         CondThermalFracLocal_Ncpus(ip1)%Val(i) = CondThermalFrac(FracbyProc(ip1)%Num(i))
      end do

   end subroutine LocalMesh_CondThermal

   ! Output:
   !  CellThermalSource_Ncpus(ip), FracThermalSource_Ncpus(ip)
   ! Use:
   !  CellbyProc(ip), CellThermalSource,
   !  FracbyProc(ip), FracThermalSource,
   SUBROUTINE LocalMesh_ThermalSource(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, i

      ip1 = ip + 1

      Nb = CellbyProc(ip1)%Pt(CellbyProc(ip1)%Nb + 1)
      ALLOCATE (CellThermalSource_Ncpus(ip1)%Val(Nb))
      DO i = 1, Nb
         CellThermalSource_Ncpus(ip1)%Val(i) = CellThermalSource(CellbyProc(ip1)%Num(i))
      END DO

      Nb = FracbyProc(ip1)%Pt(FracbyProc(ip1)%Nb + 1)
      ALLOCATE (FracThermalSource_Ncpus(ip1)%Val(Nb))
      DO i = 1, Nb
         FracThermalSource_Ncpus(ip1)%Val(i) = FracThermalSource(FracbyProc(ip1)%Num(i))
      END DO
   END SUBROUTINE LocalMesh_ThermalSource

#endif

   ! Output:
   !  NodeDatabyWellLocal_Ncpus(ip)
   ! Use:
   !  WellbyProc, NodebyWellLocal_Ncpus, NodebyWell, NodeDatabyWell
   !> \brief Contains the local data for the wells (own+ghost) of proc ip
   !! in particular %Parent contains the local number of the node Parent                   <br>
   !!               %Parent remains -1 if head node of the well
! FIXME: Refactor with the same function to treat successively injectors and producers
   subroutine LocalMesh_NodeDatabyWellLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, Nb, nnz, k, cpt, Data_cpt, numWellRes, num_parent

      ip1 = ip + 1

      !! INJ WELL
      Nb = NodebyWellInjLocal_Ncpus(ip1)%Nb
      NodeDatabyWellInjLocal_Ncpus(ip1)%Nb = Nb
      nnz = NodebyWellInjLocal_Ncpus(ip1)%Pt(Nb + 1)

      allocate (NodeDatabyWellInjLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodeDatabyWellInjLocal_Ncpus(ip1)%Num(nnz))
      allocate (NodeDatabyWellInjLocal_Ncpus(ip1)%Val(nnz))

      ! copy CSR offsets
      ! FIXME: would be shared using (shared... ?) pointers...
      NodeDatabyWellInjLocal_Ncpus(ip1)%Pt(:) = NodebyWellInjLocal_Ncpus(ip1)%Pt(:)
      NodeDatabyWellInjLocal_Ncpus(ip1)%Num(:) = NodebyWellInjLocal_Ncpus(ip1)%Num(:)

      ! NodeDatabyWell contains the data for every wells, whereas
      ! NodeDatabyWellLocal_Ncpus contains the data for the wells (own+ghost) of proc ip
      do k = 1, WellInjbyProc(ip1)%Pt(WellInjbyProc(ip1)%Nb + 1)
         numWellRes = WellInjbyProc(ip1)%Num(k)

         do cpt = 1, NodebyWellInj%Pt(numWellRes + 1) - NodebyWellInj%Pt(numWellRes)
            Data_cpt = NodeDatabyWellInjLocal_Ncpus(ip1)%Pt(k) + cpt
            NodeDatabyWellInjLocal_Ncpus(ip1)%Val(Data_cpt) = &
               NodeDatabyWellInj%Val(NodebyWellInj%Pt(numWellRes) + cpt)
            ! Change %Parent to have the local number of the node
            ! if Parent = -1 (head of the well), does not change it
            if (NodeDatabyWellInjLocal_Ncpus(ip1)%Val(Data_cpt)%Parent > -1) then
               NodeDatabyWellInjLocal_Ncpus(ip1)%Val(Data_cpt)%Parent = &
                  localbyGlobalNode(NodeDatabyWellInjLocal_Ncpus(ip1)%Val(Data_cpt)%Parent)
            endif
         enddo

#ifndef NDEBUG
         if (NodeDatabyWellInjLocal_Ncpus(ip1)%Val(NodebyWellInjLocal_Ncpus(ip1)%Pt(k + 1))%Parent /= -1) &
            call CommonMPI_abort("Inconsistent injector head")
#endif
         ! Find pointer of parent and fill ...%Val%PtParent (except for head node)
         do cpt = NodebyWellInjLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellInjLocal_Ncpus(ip1)%Pt(k + 1) - 1
            ! reservoir node index
            num_parent = NodeDatabyWellInjLocal_Ncpus(ip1)%Val(cpt)%Parent
#ifndef NDEBUG
            if (num_parent == -1) &
               call CommonMPI_abort("Inconsistent parent vertex in injector")
#endif
            ! look for parent (stored after son)
            do Data_cpt = cpt + 1, NodebyWellInjLocal_Ncpus(ip1)%Pt(k + 1)
               if (NodeDatabyWellInjLocal_Ncpus(ip1)%Num(Data_cpt) == num_parent) then
                  NodeDatabyWellInjLocal_Ncpus(ip1)%Val(cpt)%PtParent = Data_cpt
                  NodeDatabyWellInjLocal_Ncpus(ip1)%Val(cpt)%RelParent = &
                     Data_cpt - NodebyWellInjLocal_Ncpus(ip1)%Pt(k)
                  exit
               endif
            enddo
         enddo

      enddo

      !! PROD WELL
      Nb = NodebyWellProdLocal_Ncpus(ip1)%Nb
      NodeDatabyWellProdLocal_Ncpus(ip1)%Nb = Nb
      nnz = NodebyWellProdLocal_Ncpus(ip1)%Pt(Nb + 1)

      allocate (NodeDatabyWellProdLocal_Ncpus(ip1)%Pt(Nb + 1))
      allocate (NodeDatabyWellProdLocal_Ncpus(ip1)%Num(nnz))
      allocate (NodeDatabyWellProdLocal_Ncpus(ip1)%Val(nnz))

      ! copy CSR offsets
      ! FIXME: would be shared using (shared... ?) pointers...
      NodeDatabyWellProdLocal_Ncpus(ip1)%Pt(:) = NodebyWellProdLocal_Ncpus(ip1)%Pt(:)
      NodeDatabyWellProdLocal_Ncpus(ip1)%Num(:) = NodebyWellProdLocal_Ncpus(ip1)%Num(:)

      ! NodeDatabyWell contains the data for every wells, whereas
      ! NodeDatabyWellLocal_Ncpus contains the data for the wells (own+ghost) of proc ip
      do k = 1, WellProdbyProc(ip1)%Pt(WellProdbyProc(ip1)%Nb + 1)
         numWellRes = WellProdbyProc(ip1)%Num(k)

         do cpt = 1, NodebyWellProd%Pt(numWellRes + 1) - NodebyWellProd%Pt(numWellRes)
            Data_cpt = NodeDatabyWellProdLocal_Ncpus(ip1)%Pt(k) + cpt
            NodeDatabyWellProdLocal_Ncpus(ip1)%Val(Data_cpt) = &
               NodeDatabyWellProd%Val(NodebyWellProd%Pt(numWellRes) + cpt)
            ! Change %Parent to have the local number of the node
            ! if Parent = -1 (head of the well), does not change it
            if (NodeDatabyWellProdLocal_Ncpus(ip1)%Val(Data_cpt)%Parent > -1) then
               NodeDatabyWellProdLocal_Ncpus(ip1)%Val(Data_cpt)%Parent = &
                  localbyGlobalNode(NodeDatabyWellProdLocal_Ncpus(ip1)%Val(Data_cpt)%Parent)
            endif
         enddo

#ifndef NDEBUG
         if (NodeDatabyWellProdLocal_Ncpus(ip1)%Val(NodebyWellProdLocal_Ncpus(ip1)%Pt(k + 1))%Parent /= -1) &
            call CommonMPI_abort("Inconsistent producer head")
#endif

         ! Find pointer of parent and fill ...%Val%PtParent (except for head node)
         do cpt = NodebyWellProdLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellProdLocal_Ncpus(ip1)%Pt(k + 1) - 1
            ! reservoir node index
            num_parent = NodeDatabyWellProdLocal_Ncpus(ip1)%Val(cpt)%Parent
#ifndef NDEBUG
            if (num_parent == -1) &
               call CommonMPI_abort("Inconsistent parent vertex in producer")
#endif
            ! look for parent (stored after son)
            do Data_cpt = cpt + 1, NodebyWellProdLocal_Ncpus(ip1)%Pt(k + 1)
               if (NodeDatabyWellProdLocal_Ncpus(ip1)%Num(Data_cpt) == num_parent) then
                  NodeDatabyWellProdLocal_Ncpus(ip1)%Val(cpt)%PtParent = Data_cpt
                  NodeDatabyWellProdLocal_Ncpus(ip1)%Val(cpt)%RelParent = &
                     Data_cpt - NodebyWellProdLocal_Ncpus(ip1)%Pt(k)
                  exit
               endif
            enddo
         enddo

      enddo

   end subroutine LocalMesh_NodeDatabyWellLocal

   ! Output:
   !   WellbyNodeOwn_Ncpus(ip)
   ! Use:
   !   NodebyWellLocal_Ncpus(ip),
   !   NbNodeOwnS_Ncpus(ip), NbWellLocal_Ncpus(ip)
   !> \brief Fill WellbyNodeOwn_Ncpus(ip)
   !!     with the local number of well for each own node (local number) of proc ip.
   !!     The number of rows is the number own node of proc i.
   !!     If there is no well in own node i, juste take %Pt(i+1)=%Pt(i).
   subroutine LocalMesh_WellbyNodeOwn(ip)

      integer, intent(in) :: ip
      integer :: ip1, k, n, npt
      integer, dimension(:), allocatable :: tabNbWellbyNode

      ip1 = ip + 1

      ! INJ WELL

      ! %Nb
      WellInjbyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
      allocate (WellInjbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1))
      WellInjbyNodeOwn_Ncpus(ip1)%Pt(:) = 0

      allocate (tabNbWellbyNode(NbNodeOwnS_Ncpus(ip1)))

      ! Counting
      tabNbWellbyNode(:) = 0
      do k = 1, NbWellInjResS_Ncpus(ip1)
         do npt = NodebyWellInjLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellInjLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyWellInjLocal_Ncpus(ip1)%Num(npt)   ! local num of node
            if (n <= NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
               tabNbWellbyNode(n) = tabNbWellbyNode(n) + 1
            endif
         enddo
      enddo

      ! Filling of %Pt
      do n = 1, NbNodeOwnS_Ncpus(ip1)
         WellInjbyNodeOwn_Ncpus(ip1)%Pt(n + 1) = WellInjbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbWellbyNode(n)
      enddo

      ! Filling of %Num
      allocate (WellInjbyNodeOwn_Ncpus(ip1)%Num(WellInjbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1)))
      tabNbWellbyNode(:) = 0
      do k = 1, NbWellInjResS_Ncpus(ip1)
         do npt = NodebyWellInjLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellInjLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyWellInjLocal_Ncpus(ip1)%Num(npt)   ! local num of node
            if (n <= NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
               tabNbWellbyNode(n) = tabNbWellbyNode(n) + 1
               WellInjbyNodeOwn_Ncpus(ip1)%Num(WellInjbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbWellbyNode(n)) = k
            endif
         enddo
      enddo

      ! PROD WELL

      ! %Nb
      WellProdbyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
      allocate (WellProdbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1))
      WellProdbyNodeOwn_Ncpus(ip1)%Pt(:) = 0

      ! Counting
      tabNbWellbyNode(:) = 0
      do k = 1, NbWellProdResS_Ncpus(ip1)
         do npt = NodebyWellProdLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellProdLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyWellProdLocal_Ncpus(ip1)%Num(npt)   ! local num of node
            if (n <= NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
               tabNbWellbyNode(n) = tabNbWellbyNode(n) + 1
            endif
         enddo
      enddo

      ! Filling of %Pt
      do n = 1, NbNodeOwnS_Ncpus(ip1)
         WellProdbyNodeOwn_Ncpus(ip1)%Pt(n + 1) = WellProdbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbWellbyNode(n)
      enddo

      ! Filling of %Num
      allocate (WellProdbyNodeOwn_Ncpus(ip1)%Num(WellProdbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1)))
      tabNbWellbyNode(:) = 0
      do k = 1, NbWellProdResS_Ncpus(ip1)
         do npt = NodebyWellProdLocal_Ncpus(ip1)%Pt(k) + 1, NodebyWellProdLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyWellProdLocal_Ncpus(ip1)%Num(npt)   ! local num of node
            if (n <= NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
               tabNbWellbyNode(n) = tabNbWellbyNode(n) + 1
               WellProdbyNodeOwn_Ncpus(ip1)%Num(WellProdbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbWellbyNode(n)) = k
            endif
         enddo
      enddo

      deallocate (tabNbWellbyNode)

   end subroutine LocalMesh_WellbyNodeOwn

   ! Output:
   !   CellbyNodeOwn_Ncpus
   ! Use:
   !   NodebyCellLocal_Ncpus(ip)
   !   NbNodeOwnS_Ncpus(ip), NbCellLocal_Ncpus(ip)
   !> \brief Store the local number of cells surrounding own nodes of proc ip
   subroutine LocalMesh_CellbyNodeOwn(ip)

      integer, intent(in) :: ip
      integer :: ip1, k, npt, n

      integer, dimension(:), allocatable :: tabNbCellbyNode

      ip1 = ip + 1

      allocate (tabNbCellbyNode(NbNodeOwnS_Ncpus(ip1)))

      ! %Nb
      CellbyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
      allocate (CellbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1))
      CellbyNodeOwn_Ncpus(ip1)%Pt(:) = 0

      ! Counting
      tabNbCellbyNode(:) = 0
      do k = 1, NbCellResS_Ncpus(ip1)
         do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
            if (n <= NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
               tabNbCellbyNode(n) = tabNbCellbyNode(n) + 1
            endif
         enddo
      enddo

      ! Filling of %Num
      do n = 1, NbNodeOwnS_Ncpus(ip1)
         CellbyNodeOwn_Ncpus(ip1)%Pt(n + 1) = CellbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbCellbyNode(n)
      enddo

      allocate (CellbyNodeOwn_Ncpus(ip1)%Num(CellbyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1)))

      ! Filling
      tabNbCellbyNode(:) = 0
      do k = 1, NbCellResS_Ncpus(ip1)
         do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
            n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
            if (n <= NbNodeOwnS_Ncpus(ip1)) then
               tabNbCellbyNode(n) = tabNbCellbyNode(n) + 1
               CellbyNodeOwn_Ncpus(ip1)%Num(CellbyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbCellbyNode(n)) = k
            endif
         enddo
      enddo

      deallocate (tabNbCellbyNode)

   end subroutine LocalMesh_CellbyNodeOwn

   ! Output:
   !  NodebyNodeOwn(ip)
   ! Use:
   !  NodebyCellLocal_Ncpus(ip), CellbyNodeOwn_Ncpus(ip)
   !> \brief Store the neighbour nodes (local) of own nodes
   !!    (we know neighbour nodes are all own or ghost for this proc)
   !!    CAREFUL: the nodes have no particulary order (own nodes are not at the beginning !)
   subroutine LocalMesh_NodebyNodeOwn(ip)

      integer, intent(in) :: ip

      integer :: ip1

      integer m, k, kpt, n, npt
      integer, dimension(:), allocatable :: &
         colorNodeLocal, & ! used to mark if the node is considered or not
         tabNbNodebyNode

      ip1 = ip + 1

      allocate (colorNodeLocal(NbNodeResS_Ncpus(ip1)))
      allocate (tabNbNodebyNode(NbNodeOwnS_Ncpus(ip1)))

      colorNodeLocal(:) = 0
      tabNbNodebyNode(:) = 0

      ! own nodes of proc ip1
      do m = 1, NbNodeOwnS_Ncpus(ip1)

         ! loop over cells surrounding cell m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            ! loop over nodes of cell k
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               colorNodeLocal(n) = m
            enddo
         enddo

         ! number of neighbour nodes of the node own m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorNodeLocal(n) == m) then
                  tabNbNodebyNode(m) = tabNbNodebyNode(m) + 1
                  colorNodeLocal(n) = 0
               endif
            enddo
         enddo

      enddo

      NodebyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
      allocate (NodebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1))

      ! Filling of %Pt
      NodebyNodeOwn_Ncpus(ip1)%Pt(1) = 0
      do m = 1, NbNodeOwnS_Ncpus(ip1)
         NodebyNodeOwn_Ncpus(ip1)%Pt(m + 1) = NodebyNodeOwn_Ncpus(ip1)%Pt(m) + tabNbNodebyNode(m)
      enddo

      ! Filling of %Num
      allocate (NodebyNodeOwn_Ncpus(ip1)%Num(NodebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1)))

      colorNodeLocal(:) = 0
      tabNbNodebyNode(:) = 0

      do m = 1, NbNodeOwnS_Ncpus(ip1)

         ! loop over cells surrounding cell m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            ! loop over nodes of cell k
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               colorNodeLocal(n) = m
            enddo
         enddo
         ! filling of number of neighbour nodes of node own m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorNodeLocal(n) == m) then
                  tabNbNodebyNode(m) = tabNbNodebyNode(m) + 1
                  colorNodeLocal(n) = 0
                  NodebyNodeOwn_Ncpus(ip1)%Num(NodebyNodeOwn_Ncpus(ip1)%Pt(m) + tabNbNodebyNode(m)) = n
               endif
            enddo
         enddo
      enddo

      deallocate (colorNodeLocal)
      deallocate (tabNbNodebyNode)

   end subroutine LocalMesh_NodebyNodeOwn

   ! Output:
   !  FracbyProc_Ncpus(ip)
   ! Use:
   !  FacebyProc_Ncpus(ip),IdFace
   !> \brief Compress FacebyProc_Ncpus(ip) using IdFace
   !! FracbyProc is considered as a subset of FacebyProc, csr.
   !!   The number of rows are same as that of FacebyProc.
   !!   If there is no frac in row i, juste take %Pt(i+1)=%Pt(i).
   !!   It's not same with *byProc, no empty row in, *=cell/face/node
   subroutine LocalMesh_FracbyProc(ip)

      integer, intent(in) :: ip

      integer :: ip1, nf, j, i, Nb
      ip1 = ip + 1

      ! FracbyProc%Nb
      Nb = FacebyProc(ip1)%Nb
      FracbyProc(ip1)%Nb = Nb

      ! allocate
      allocate (FracbyProc(ip1)%Pt(Nb + 1)) ! allo %Pt of proc1
      FracbyProc(ip1)%Pt(1) = 0

      ! loop csr FracbyProc%Pt
      ! nf: nb of frac, used to make FracbyProc(ip)%Num
      nf = 0
      do i = 1, FacebyProc(ip1)%Nb
         do j = FacebyProc(ip1)%Pt(i) + 1, FacebyProc(ip1)%Pt(i + 1)
            ! if FacebyProc(ip1)%Num(j) is frac
            if (IdFace(FacebyProc(ip1)%Num(j)) == -2) then
               nf = nf + 1
            end if
         end do
         FracbyProc(ip1)%Pt(i + 1) = nf
      end do

      !print *, "DEBUG - Registered", nf, "fractures"
      ! rq: if no frac, nf=0, %Pt=0

      ! loop csr for FracbyProc(ip1)%Num
      allocate (FracbyProc(ip1)%Num(nf)) ! here nf is equal to the size of %Num

      nf = 0
      do i = 1, FacebyProc(ip1)%Nb
         do j = FacebyProc(ip1)%Pt(i) + 1, FacebyProc(ip1)%Pt(i + 1)
            ! if FacebyProc(ip1)%Num(j) is frac
            if (IdFace(FacebyProc(ip1)%Num(j)) == -2) then
               nf = nf + 1
               FracbyProc(ip1)%Num(nf) = FacebyProc(ip1)%Num(j)
            end if
         end do
      end do

      ! Not all elements of %val(:) are useful. Since in some rows, no frac
      allocate (FracbyProc(ip1)%Val(Nb))
      FracbyProc(ip1)%Val = FacebyProc(ip1)%Val

      ! NbFracResS_Ncpus(ip1), NbFracOwnS_Ncpus(ip1), rq: =0 if no frac
      NbFracResS_Ncpus(ip1) = FracbyProc(ip1)%Pt(FracbyProc(ip1)%Nb + 1) - FracbyProc(ip1)%Pt(1)
      NbFracOwnS_Ncpus(ip1) = FracbyProc(ip1)%Pt(2) - FracbyProc(ip1)%Pt(1)

   end subroutine LocalMesh_FracbyProc

   ! Output:
   !  NodebyFracOwn_Npus(ip)
   ! Use:
   !  CellbyFracOwn_Ncpus(ip), NodebyCellLocal
   !> \brief Store the neighbour cells (local) of own frac
   subroutine LocalMesh_NodebyFracOwn(ip)

      ! Calcule les noeuds voisin de chaque frac own (on sait qu'ils sont tous own ou ghost au proc)
      ! voisins ordonnes dans un ordre quelconque (pas les own en premier!)

      integer, intent(in) :: ip

      integer :: ip1

      integer m, k, kpt, n, npt
      integer, dimension(:), allocatable :: &
         colorNodeLocal, & ! used to mark if the node is considered or not
         tabNbNodebyFrac

      ip1 = ip + 1

      allocate (colorNodeLocal(NbNodeResS_Ncpus(ip1)))
      allocate (tabNbNodebyFrac(NbFracOwnS_Ncpus(ip1)))

      colorNodeLocal(:) = 0
      tabNbNodebyFrac(:) = 0

      ! boucle sur les noeuds own
      do m = 1, NbFracOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les noeuds n de k
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               colorNodeLocal(n) = m ! this can be set more than once
            enddo
         enddo

         ! number of neighbour nodes of the node own m
         ! we have this loop because colorNodeLocal can be set more than once
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorNodeLocal(n) == m) then
                  tabNbNodebyFrac(m) = tabNbNodebyFrac(m) + 1
                  colorNodeLocal(n) = 0
               endif
            enddo
         enddo

      enddo

      NodebyFracOwn_Ncpus(ip1)%Nb = NbFracOwnS_Ncpus(ip1)
      allocate (NodebyFracOwn_Ncpus(ip1)%Pt(NbFracOwnS_Ncpus(ip1) + 1))

      ! Filling of %Pt
      NodebyFracOwn_Ncpus(ip1)%Pt(1) = 0
      do m = 1, NbFracOwnS_Ncpus(ip1)
         NodebyFracOwn_Ncpus(ip1)%Pt(m + 1) = NodebyFracOwn_Ncpus(ip1)%Pt(m) + tabNbNodebyFrac(m)
      enddo

      ! Filling of %Num
      allocate (NodebyFracOwn_Ncpus(ip1)%Num(NodebyFracOwn_Ncpus(ip1)%Pt(NbFracOwnS_Ncpus(ip1) + 1)))

      colorNodeLocal(:) = 0
      tabNbNodebyFrac(:) = 0

      do m = 1, NbFracOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les noeuds n de k
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               colorNodeLocal(n) = m
            enddo
         enddo
         ! remplissage du nombre de noeuds voisins pour le noeud own m
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = NodebyCellLocal_Ncpus(ip1)%Pt(k) + 1, NodebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = NodebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorNodeLocal(n) == m) then
                  tabNbNodebyFrac(m) = tabNbNodebyFrac(m) + 1
                  colorNodeLocal(n) = 0
                  NodebyFracOwn_Ncpus(ip1)%Num(NodebyFracOwn_Ncpus(ip1)%Pt(m) + tabNbNodebyFrac(m)) = n
               endif
            enddo
         enddo
      enddo

      deallocate (colorNodeLocal)
      deallocate (tabNbNodebyFrac)

   end subroutine LocalMesh_NodebyFracOwn

   ! Output:
   !   CellbyFracOwn_Npus(ip)
   ! Use:
   !   CellbyFaceLocal_Ncpus(ip), IdFaceRes_Ncpus(ip)
   subroutine LocalMesh_CellbyFracOwn(ip)

      integer, intent(in) :: ip
      integer :: ip1

      integer :: Nb, nf, m, Nnz, &
                 startfrac, endfrac, startcell, endcell

      ip1 = ip + 1

      ! %Nb
      Nb = NbFracOwnS_Ncpus(ip1)
      CellbyFracOwn_Ncpus(ip1)%Nb = Nb

      ! %Pt
      allocate (CellbyFracOwn_Ncpus(ip1)%Pt(Nb + 1))

      nf = 1
      CellbyFracOwn_Ncpus(ip1)%Pt(1) = 0

      do m = 1, NbFaceOwnS_Ncpus(ip1) ! face own

         if (IdFaceRes_Ncpus(ip1)%Val(m) == -2) then ! frac,

            CellbyFracOwn_Ncpus(ip1)%Pt(nf + 1) = CellbyFracOwn_Ncpus(ip1)%Pt(nf) &
                                                  + CellbyFaceLocal_Ncpus(ip1)%Pt(m + 1) &
                                                  - CellbyFaceLocal_Ncpus(ip1)%Pt(m)
            nf = nf + 1
         end if
      end do

      ! %Num
      Nnz = CellbyFracOwn_Ncpus(ip1)%Pt(Nb + 1)
      allocate (CellbyFracOwn_Ncpus(ip1)%Num(Nnz))

      nf = 1
      do m = 1, NbFaceOwnS_Ncpus(ip1) ! face own

         if (IdFaceRes_Ncpus(ip1)%Val(m) == -2) then ! frac, copy %Num

            ! number of non zero this line
            startfrac = CellbyFracOwn_Ncpus(ip1)%Pt(nf) + 1
            endfrac = CellbyFracOwn_Ncpus(ip1)%Pt(nf + 1)

            startcell = CellbyFaceLocal_Ncpus(ip1)%Pt(m) + 1
            endcell = CellbyFaceLocal_Ncpus(ip1)%Pt(m + 1)

            CellbyFracOwn_Ncpus(ip1)%Num(startfrac:endfrac) = &
               CellbyFaceLocal_Ncpus(ip1)%Num(startcell:endcell)

            nf = nf + 1
         end if
      end do

   end subroutine LocalMesh_CellbyFracOwn

   ! Output:
   !   FracbyCellLocal_Ncpu(ip)
   ! Use:
   !   FacebyCellLocal_Ncpus(ip), IdFaceRes_Ncpus(ip)
   subroutine LocalMesh_FracbyCellLocal(ip)

      integer, intent(in) :: ip
      integer :: ip1, nf, j, i, Nb

      ip1 = ip + 1

      ! %Nb
      Nb = FacebyCellLocal_Ncpus(ip1)%Nb
      FracbyCellLocal_Ncpus(ip1)%Nb = Nb

      allocate (FracbyCellLocal_Ncpus(ip1)%Pt(Nb + 1)) ! allo %Pt of proc1
      FracbyCellLocal_Ncpus(ip1)%Pt(1) = 0

      ! loop csr FracbyCellLocal(ip1)
      ! nf: nb of frac, used to make FracbyProc(ip)%Num
      nf = 0
      do i = 1, FacebyCellLocal_Ncpus(ip1)%Nb

         do j = FacebyCellLocal_Ncpus(ip1)%Pt(i) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(i + 1)

            ! if FacebyCellLocal_Ncpus(ip1)%Num(j) is frac
            if (IdFaceRes_Ncpus(ip1)%Val(FacebyCellLocal_Ncpus(ip1)%Num(j)) == -2) then
               nf = nf + 1
            end if
         end do

         FracbyCellLocal_Ncpus(ip1)%Pt(i + 1) = nf
      end do

      ! loop csr for FracbyProc(ip1)%Num
      allocate (FracbyCellLocal_Ncpus(ip1)%Num(nf)) ! here nf is equal to the size of %Num

      nf = 0
      do i = 1, FacebyCellLocal_Ncpus(ip1)%Nb

         do j = FacebyCellLocal_Ncpus(ip1)%Pt(i) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(i + 1)

            ! if FacebyProc(ip1)%Num(j) is frac
            if (IdFaceRes_Ncpus(ip1)%Val(FacebyCellLocal_Ncpus(ip1)%Num(j)) == -2) then
               nf = nf + 1
               FracbyCellLocal_Ncpus(ip1)%Num(nf) = FacebyCellLocal_Ncpus(ip1)%Num(j)
            end if
         end do
      end do

   end subroutine LocalMesh_FracbyCellLocal

   ! Output:
   !  FacebyNodeOwn(ip)
   ! Use:
   !  FacebyCellLocal_Ncpus(ip), CellbyNodeOwn_Ncpus(ip)
   subroutine LocalMesh_FacebyNodeOwn(ip)

      ! Calcule les faces voisin de chaque noeud own (on sait qu'ils sont tous own ou ghost au proc)
      ! voisins ordonnes dans un ordre quelconque (pas les own en premier!)

      integer, intent(in) :: ip

      integer :: ip1

      integer m, k, kpt, n, npt
      integer, dimension(:), allocatable :: &
         colorFaceLocal, & ! used to mark if the node is considered or not
         tabNbFacebyNode

      ip1 = ip + 1

      allocate (colorFaceLocal(NbFaceResS_Ncpus(ip1)))
      allocate (tabNbFacebyNode(NbNodeOwnS_Ncpus(ip1)))

      colorFaceLocal(:) = 0
      tabNbFacebyNode(:) = 0

      ! boucle sur les noeuds own
      do m = 1, NbNodeOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les noeuds n de k
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = FacebyCellLocal_Ncpus(ip1)%Pt(k) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FacebyCellLocal_Ncpus(ip1)%Num(npt)
               colorFaceLocal(n) = m
            enddo
         enddo

         ! number of neighbour faces of the node own m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = FacebyCellLocal_Ncpus(ip1)%Pt(k) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FacebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorFaceLocal(n) == m) then
                  tabNbFacebyNode(m) = tabNbFacebyNode(m) + 1
                  colorFaceLocal(n) = 0
               endif
            enddo
         enddo

      enddo

      FacebyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
      allocate (FacebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1))

      ! Filling of %Pt
      FacebyNodeOwn_Ncpus(ip1)%Pt(1) = 0
      do m = 1, NbNodeOwnS_Ncpus(ip1)
         FacebyNodeOwn_Ncpus(ip1)%Pt(m + 1) = FacebyNodeOwn_Ncpus(ip1)%Pt(m) + tabNbFacebyNode(m)
      enddo

      ! Filling of %Num
      allocate (FacebyNodeOwn_Ncpus(ip1)%Num(FacebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1) + 1)))

      colorFaceLocal(:) = 0
      tabNbFacebyNode(:) = 0

      do m = 1, NbNodeOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les noeuds n de k
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = FacebyCellLocal_Ncpus(ip1)%Pt(k) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FacebyCellLocal_Ncpus(ip1)%Num(npt)
               colorFaceLocal(n) = m
            enddo
         enddo
         ! remplissage du nombre de noeuds voisins pour le noeud own m
         do kpt = CellbyNodeOwn_Ncpus(ip1)%Pt(m) + 1, CellbyNodeOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyNodeOwn_Ncpus(ip1)%Num(kpt)
            do npt = FacebyCellLocal_Ncpus(ip1)%Pt(k) + 1, FacebyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FacebyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorFaceLocal(n) == m) then
                  tabNbFacebyNode(m) = tabNbFacebyNode(m) + 1
                  colorFaceLocal(n) = 0
                  FacebyNodeOwn_Ncpus(ip1)%Num(FacebyNodeOwn_Ncpus(ip1)%Pt(m) + tabNbFacebyNode(m)) = n
               endif
            enddo
         enddo
      enddo

      deallocate (colorFaceLocal)
      deallocate (tabNbFacebyNode)

   end subroutine LocalMesh_FacebyNodeOwn

   ! ! Rq important:
   ! !   The folllowing way to make FacebyNodeOwn is wrong
   ! !   ex: node 1 is in cell 2, the following way gives the faces which contains the node 1. Not all the faces in cell 2 could be contained.

   ! ! Output:
   ! !   FacebyNodeOwn_Ncpus(ip)
   ! ! Use:
   ! !   NodebyFaceLocal_Ncpus(ip), NbNodeOwnS_Ncpus(ip)
   ! subroutine LocalMesh_FacebyNodeOwn(ip)

   !   integer, intent(in) :: ip

   !   integer :: &
   !        ip1, & ! = ip+1
   !        k, npt, n

   !   integer, dimension(:), allocatable :: &
   !        tabNbFacebyNode

   !   ip1 = ip + 1

   !   allocate(tabNbFacebyNode(NbNodeOwnS_Ncpus(ip1)))

   !   ! %Nb
   !   FacebyNodeOwn_Ncpus(ip1)%Nb = NbNodeOwnS_Ncpus(ip1)
   !   allocate(FacebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1)+1))

   !   ! Counting
   !   tabNbFacebyNode(:) = 0
   !   do k = 1,NbFaceResS_Ncpus(ip1)
   !      do npt = NodebyFaceLocal_Ncpus(ip1)%Pt(k)+1,NodebyFaceLocal_Ncpus(ip1)%Pt(k+1)
   !         n = NodebyFaceLocal_Ncpus(ip1)%Num(npt)
   !         if (n<=NbNodeOwnS_Ncpus(ip1)) then ! node = (node own, node ghost)
   !            tabNbFacebyNode(n) = tabNbFacebyNode(n) + 1
   !         endif
   !      enddo
   !   enddo

   !   FacebyNodeOwn_Ncpus(ip1)%Pt(1) = 0
   !   do n = 1,NbNodeOwnS_Ncpus(ip1)
   !      FacebyNodeOwn_Ncpus(ip1)%Pt(n+1) = FacebyNodeOwn_Ncpus(ip1)%Pt(n) + tabNbFacebyNode(n)
   !   enddo

   !   allocate(FacebyNodeOwn_Ncpus(ip1)%Num(FacebyNodeOwn_Ncpus(ip1)%Pt(NbNodeOwnS_Ncpus(ip1)+1)))

   !   ! Filling
   !   tabNbFacebyNode(:) = 0
   !   do k = 1,NbFaceResS_Ncpus(ip1)
   !      do npt = NodebyFaceLocal_Ncpus(ip1)%Pt(k)+1,NodebyFaceLocal_Ncpus(ip1)%Pt(k+1)
   !         n = NodebyFaceLocal_Ncpus(ip1)%Num(npt)
   !         if (n<=NbNodeOwnS_Ncpus(ip1)) then
   !            tabNbFacebyNode(n) = tabNbFacebyNode(n) + 1
   !            FacebyNodeOwn_Ncpus(ip1)%Num(FacebyNodeOwn_Ncpus(ip1)%Pt(n)+tabNbFacebyNode(n)) = k
   !         endif
   !      enddo
   !   enddo

   !   deallocate(tabNbFacebyNode )

   ! end subroutine LocalMesh_FacebyNodeOwn

   ! Output:
   !   FracbyNodeOwn_Ncpus
   ! Use:
   !   FacebyNodeOwn_Ncpus, IdFaceRes_Ncpus(ip)
   !> \brief Contains the local number of frac for each own node (local number) of proc ip.
   !!     The number of rows is the number of own nodes of proc i.
   !!     If there is no frac in own node i, juste take %Pt(i+1)=%Pt(i).
   subroutine LocalMesh_FracbyNodeOwn(ip)

      integer, intent(in) :: ip
      integer :: ip1, nf, j, i, Nb

      ip1 = ip + 1

      ! %Nb
      Nb = FacebyNodeOwn_Ncpus(ip1)%Nb
      FracbyNodeOwn_Ncpus(ip1)%Nb = Nb

      allocate (FracbyNodeOwn_Ncpus(ip1)%Pt(Nb + 1)) ! allo %Pt of proc1
      FracbyNodeOwn_Ncpus(ip1)%Pt(1) = 0

      ! loop csr FacebyNodeOwn(ip1)
      ! nf: nb of frac, used to make FracbyProc(ip)%Num and %Pt
      nf = 0
      do i = 1, FacebyNodeOwn_Ncpus(ip1)%Nb

         do j = FacebyNodeOwn_Ncpus(ip1)%Pt(i) + 1, FacebyNodeOwn_Ncpus(ip1)%Pt(i + 1)

            ! if FracbyNodeOwn_Ncpus(ip1)%Num(j) is frac
            if (IdFaceRes_Ncpus(ip1)%Val(FacebyNodeOwn_Ncpus(ip1)%Num(j)) == -2) then
               nf = nf + 1
            end if
         end do

         FracbyNodeOwn_Ncpus(ip1)%Pt(i + 1) = nf
      end do

      ! loop csr for FracbyProc(ip1)%Num
      allocate (FracbyNodeOwn_Ncpus(ip1)%Num(nf)) ! here nf is equal to the size of %Num

      nf = 0
      do i = 1, FacebyNodeOwn_Ncpus(ip1)%Nb

         do j = FacebyNodeOwn_Ncpus(ip1)%Pt(i) + 1, FacebyNodeOwn_Ncpus(ip1)%Pt(i + 1)

            ! if frac
            if (IdFaceRes_Ncpus(ip1)%Val(FacebyNodeOwn_Ncpus(ip1)%Num(j)) == -2) then
               nf = nf + 1
               FracbyNodeOwn_Ncpus(ip1)%Num(nf) = FacebyNodeOwn_Ncpus(ip1)%Num(j)
            end if
         end do
      end do

   end subroutine LocalMesh_FracbyNodeOwn

   ! Output:
   !  FracbyFracOwn_Ncpus(ip1)
   ! Use:
   !  CellbyFracOwn_Ncpus(ip1), FracbyCellLocal_Ncpus(ip1)
   !> \brief Store the neighbour frac (local number) of own frac.
   !!     The number of rows is the number of own frac of proc i.
   subroutine LocalMesh_FracbyFracOwn(ip)

      integer, intent(in) :: ip

      integer :: ip1

      integer m, k, kpt, n, npt
      integer, dimension(:), allocatable :: &
         colorFracLocal, & ! used to mark if the frac is considered or not
         tabNbFracbyFrac

      ip1 = ip + 1

      allocate (colorFracLocal(NbFaceResS_Ncpus(ip1)))
      allocate (tabNbFracbyFrac(NbFracOwnS_Ncpus(ip1)))

      colorFracLocal(:) = 0
      tabNbFracbyFrac(:) = 0

      ! boucle sur les frac own
      do m = 1, NbFracOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les frac n de k
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = FracbyCellLocal_Ncpus(ip1)%Pt(k) + 1, FracbyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FracbyCellLocal_Ncpus(ip1)%Num(npt)
               colorFracLocal(n) = m
            enddo
         enddo

         ! number of neighbour frac of the frac own m
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = FracbyCellLocal_Ncpus(ip1)%Pt(k) + 1, FracbyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FracbyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorFracLocal(n) == m) then
                  tabNbFracbyFrac(m) = tabNbFracbyFrac(m) + 1
                  colorFracLocal(n) = 0
               endif
            enddo
         enddo

      enddo

      FracbyFracOwn_Ncpus(ip1)%Nb = NbFracOwnS_Ncpus(ip1)
      allocate (FracbyFracOwn_Ncpus(ip1)%Pt(NbFracOwnS_Ncpus(ip1) + 1))

      FracbyFracOwn_Ncpus(ip1)%Pt(1) = 0
      do m = 1, NbFracOwnS_Ncpus(ip1)
         FracbyFracOwn_Ncpus(ip1)%Pt(m + 1) = FracbyFracOwn_Ncpus(ip1)%Pt(m) + tabNbFracbyFrac(m)
      enddo

      ! Filling of %Num
      allocate (FracbyFracOwn_Ncpus(ip1)%Num(FracbyFracOwn_Ncpus(ip1)%Pt(NbFracOwnS_Ncpus(ip1) + 1)))

      colorFracLocal(:) = 0
      tabNbFracbyFrac(:) = 0

      do m = 1, NbFracOwnS_Ncpus(ip1)

         ! boucle sur les mailles k voisines de m puis sur les frac n de k
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = FracbyCellLocal_Ncpus(ip1)%Pt(k) + 1, FracbyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FracbyCellLocal_Ncpus(ip1)%Num(npt)
               colorFracLocal(n) = m
            enddo
         enddo
         ! remplissage du nombre de noeuds voisins pour le frac own m
         do kpt = CellbyFracOwn_Ncpus(ip1)%Pt(m) + 1, CellbyFracOwn_Ncpus(ip1)%Pt(m + 1)
            k = CellbyFracOwn_Ncpus(ip1)%Num(kpt)
            do npt = FracbyCellLocal_Ncpus(ip1)%Pt(k) + 1, FracbyCellLocal_Ncpus(ip1)%Pt(k + 1)
               n = FracbyCellLocal_Ncpus(ip1)%Num(npt)
               if (colorFracLocal(n) == m) then
                  tabNbFracbyFrac(m) = tabNbFracbyFrac(m) + 1
                  colorFracLocal(n) = 0
                  FracbyFracOwn_Ncpus(ip1)%Num(FracbyFracOwn_Ncpus(ip1)%Pt(m) + tabNbFracbyFrac(m)) = n
               endif
            enddo
         enddo
      enddo

      deallocate (colorFracLocal)
      deallocate (tabNbFracbyFrac)

   end subroutine LocalMesh_FracbyFracOwn

   ! FracToFaceLocal_Ncpus(ip1)
   ! Use IdFaceRes_Ncpus(ip1)
   subroutine LocalMesh_FracToFace(ip)

      integer, intent(in) :: ip
      integer :: ip1, j, i, Nb

      ip1 = ip + 1

      Nb = NbFracResS_Ncpus(ip1)
      allocate (FracToFaceLocal_Ncpus(ip1)%Val(Nb))
      FracToFaceLocal_Ncpus(ip1)%Val(:) = 0

      j = 1
      do i = 1, NbFaceResS_Ncpus(ip1) ! loop of face
         if (IdFaceRes_Ncpus(ip1)%Val(i) == -2) then ! if frac
            FracToFaceLocal_Ncpus(ip1)%Val(j) = i
            j = j + 1
         end if
      end do

   end subroutine LocalMesh_FracToFace

   ! FaceToFracLocal_Ncpus(ip1)
   ! Use IdFaceRes_Ncpus(ip1)
   subroutine LocalMesh_FaceToFrac(ip)

      integer, intent(in) :: ip
      integer :: ip1, j, i, Nb

      ip1 = ip + 1

      Nb = NbFaceResS_Ncpus(ip1)
      allocate (FaceToFracLocal_Ncpus(ip1)%Val(Nb))
      FaceToFracLocal_Ncpus(ip1)%Val(:) = 0

      j = 1
      do i = 1, Nb ! loop of face
         if (IdFaceRes_Ncpus(ip1)%Val(i) == -2) then ! if frac
            FaceToFracLocal_Ncpus(ip1)%Val(i) = j
            j = j + 1
         end if
      end do

   end subroutine LocalMesh_FaceToFrac

   ! ************************************************************** !

#ifdef NDEBUG
   pure &
#endif
      subroutine LocalMesh_compute_local_info(part_global_info, total_nb_elements, part_nb_owns, part_local_info)
      type(CSR), dimension(:), target, intent(in) :: part_global_info
      integer, intent(in) :: total_nb_elements
      integer, intent(in) :: part_nb_owns(:)
      type(FamilyDOFIdCOC), allocatable, dimension(:), target, intent(out) :: part_local_info
      integer :: nb_parts, nb_neighbors, nb_nodes, i, k, proc, proc_k
      type(FamilyDOFId), allocatable, dimension(:) :: GtoL_map
      integer(c_int), pointer, dimension(:) :: offsets

      allocate (GtoL_map(total_nb_elements))
      call LocalMesh_make_GtoL_map(part_global_info, part_nb_owns, GtoL_map)

      nb_parts = size(part_global_info)
      allocate (part_local_info(nb_parts))

      do proc = 1, nb_parts

         call copy_sparsity_pattern(part_global_info(proc), part_local_info(proc))

         nb_nodes = size(part_local_info(proc)%ids)
         do i = 1, nb_nodes
            part_local_info(proc)%ids(i)%local_id = GtoL_map(part_global_info(proc)%Num(i))%local_id
         end do

         ! FIXME: part_global_info(proc)%Val has not the same size as NodebuProc%Num what is misleading!
         !        part_global_info(proc)%Val actually stores the real proc id in a small vector (cf. code below)
         offsets => part_local_info(proc)%offsets
         nb_neighbors = part_global_info(proc)%Nb ! including proc
         do k = 1, nb_neighbors
            proc_k = part_global_info(proc)%Val(k)
            do i = offsets(k) + 1, offsets(k + 1)
               part_local_info(proc)%ids(i)%proc = proc_k
#ifndef NDEBUG
               if (GtoL_map(part_global_info(proc)%Num(i))%proc /= proc_k) &
                  call CommonMPI_abort("Inconsistent proc id")
#endif
            end do
         end do

      end do

      deallocate (GtoL_map)

   end subroutine LocalMesh_compute_local_info

#ifdef NDEBUG
   pure &
#endif
      subroutine LocalMesh_make_GtoL_map(partition, nb_owns, GtoL_map)
      type(CSR), dimension(:), intent(in) :: partition
      integer, intent(in) :: nb_owns(:)
      type(FamilyDOFId), allocatable, dimension(:), intent(inout) :: GtoL_map

      integer :: i, proc, global_id
      type(FamilyDOFId) :: default_dof_id

#ifndef NDEBUG
      if (.not. allocated(GtoL_map)) &
         call CommonMPI_abort("GtoL_map should already be a valid pointer")
#endif

      default_dof_id%proc = -1
      default_dof_id%local_id = 0
      do i = 1, size(GtoL_map)
         GtoL_map(i) = default_dof_id
      end do

      do proc = 1, size(partition)
#ifndef NDEBUG
         if (nb_owns(proc) > partition(proc)%Pt(partition(proc)%Nb + 1)) then
            write (*, *) "Weird partition: proc", proc, "nb owns", nb_owns(proc), "csr Nb", partition(proc)%Nb
            call CommonMPI_abort("Inconsistent partition")
         end if
#endif
         do i = 1, nb_owns(proc)
            global_id = partition(proc)%Num(i)
            GtoL_map(global_id)%proc = proc - 1 ! C indexing, master proc rank is 0
            ! the ith DOF is owned by proc (i <= nb_owns(proc))
            GtoL_map(global_id)%local_id = i
         end do
      end do

   end subroutine LocalMesh_make_GtoL_map

end module LocalMesh
