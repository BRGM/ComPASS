!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module MeshSchema
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding, only: c_int, c_double, c_size_t, c_null_ptr, c_loc, c_ptr, c_int8_t, c_bool
   use InteroperabilityStructures, only: cpp_array_wrapper, retrieve_csr_matrix, csr_matrix_wrapper
   use mpi
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort

   use CommonType
   use DefModel
   use DefWell
   use DefMSWell
   use MeshInfo
   use LocalMesh
#ifdef _WITH_FREEFLOW_STRUCTURES_
   use FreeFlowTypes
#endif

#else
   use iso_c_binding, only: c_int, c_double, c_size_t, c_null_ptr, c_loc, c_ptr, c_int8_t, c_bool
   use InteroperabilityStructures, only: cpp_array_wrapper, retrieve_csr_matrix, csr_matrix_wrapper
   use mpi
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort

   use CommonType, only: &
      Type_IdNode, CSR, &
      VALSIZE_NB, VALSIZE_NNZ, VALSIZE_ZERO, &
      FamilyDOFId, FamilyDOFIdCOC, check_FamilyDOFIdCOC, create_FamilyDOFId_MPI_struct, &
      CommonType_deallocCSR, CommonType_csrcopy, free_ComPASS_struct
   use DefModel, only: IndThermique, NbPhase, NbComp

   use DefWell, only: &
      TYPE_CSRDataNodeWell, PerforationDataCSR_wrapper, WellData_type, &
      DefWell_mpi_register_well_data_description, DefWell_csrdatawellcopy, &
      DefWell_deallocCSRDataWell

   use DefMSWell, only: &
      MSWellData_type, DefMSWell_mpi_register_mswell_data_description
   use MeshInfo, only: &
      PartInfo

   use LocalMesh, only: &
      XNodeRes_Ncpus, &
      meshSizeS_xmin, meshSizeS_xmax, meshSizeS_ymin, meshSizeS_ymax, meshSizeS_zmin, meshSizeS_zmax, &
      NodeFlags_Ncpus, CellFlags_Ncpus, FaceFlags_Ncpus, &
      CellTypes_Ncpus, FaceTypes_Ncpus, &
      NodeRocktype_Ncpus, CellRocktype_Ncpus, FracRocktype_Ncpus, &
      NbCellResS_Ncpus, NbFaceResS_Ncpus, NbNodeResS_Ncpus, NbFracResS_Ncpus, &
      NbCellOwnS_Ncpus, NbFaceOwnS_Ncpus, NbNodeOwnS_Ncpus, NbFracOwnS_Ncpus, &
      NbWellInjResS_Ncpus, NbWellProdResS_Ncpus, &
      NbWellInjOwnS_Ncpus, NbWellProdOwnS_Ncpus, &
      NbMSWellResS_Ncpus, NbMSWellOwnS_Ncpus, &
      DataWellInjRes_Ncpus, DataWellProdRes_Ncpus, &
      DataMSWellRes_Ncpus, &
      IdCellRes_Ncpus, IdFaceRes_Ncpus, IdNodeRes_Ncpus, &
      FracToFaceLocal_Ncpus, FaceToFracLocal_Ncpus, &
      PorositeCell_Ncpus, PorositeFrac_Ncpus, &
      PermCellLocal_Ncpus, PermFracLocal_Ncpus, &
      CondThermalCellLocal_Ncpus, CondThermalFracLocal_Ncpus, &
      CellThermalSource_Ncpus, FracThermalSource_Ncpus, &
      NodebyNodeOwn_Ncpus, FracbyNodeOwn_Ncpus, CellbyNodeOwn_Ncpus, &
      NodebyFracOwn_Ncpus, CellbyFracOwn_Ncpus, FracbyFracOwn_Ncpus, &
      FacebyCellLocal_Ncpus, FracbyCellLocal_Ncpus, NodebyCellLocal_Ncpus, NodebyFaceLocal_Ncpus, &
      WellInjbyNodeOwn_Ncpus, WellProdbyNodeOwn_Ncpus, MSWellbyNodeOwn_Ncpus, &
      NodebyWellInjLocal_Ncpus, NodebyWellProdLocal_Ncpus, NodebyMSWellLocal_Ncpus, &
      NumNodebyProc_Ncpus, NumFracbyProc_Ncpus, &
      NumWellInjbyProc_Ncpus, NumWellProdbyProc_Ncpus, &
      NodeDatabyWellInjLocal_Ncpus, NodeDatabyWellProdLocal_Ncpus, NodeDatabyMSWellLocal_Ncpus

#ifdef _WITH_FREEFLOW_STRUCTURES_
   use FreeFlowTypes, only: TYPE_FFfarfield
#endif

#endif
   implicit none

   ! 1. Mesh Info
   integer(c_int), allocatable, dimension(:), target, public :: &
      NbCellLocal_Ncpus, NbCellOwn_Ncpus, &
      NbFaceLocal_Ncpus, NbFaceOwn_Ncpus, &
      NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbFracLocal_Ncpus, NbFracOwn_Ncpus

   double precision, protected :: &
      meshSize_xmax, & !< Global size of mesh xmax, known by all CPUs
      meshSize_xmin, & !< Global size of mesh xmin, known by all CPUs
      meshSize_ymax, & !< Global size of mesh ymax, known by all CPUs
      meshSize_ymin, & !< Global size of mesh ymin, known by all CPUs
      meshSize_zmax, & !< Global size of mesh zmax, known by all CPUs
      meshSize_zmin    !< Global size of mesh zmin, known by all CPUs

   ! 2. Mesh and connectivity
   !    Used to (1) calcul Transmissivities: TkLocal and TkFracLocal
   !            (2) Assembly
   type(CSR), public :: &
      NodebyNodeOwn, &   ! (1)
      FracbyNodeOwn, &   ! (1)
      CellbyNodeOwn, &   ! (1)
      !
      NodebyFracOwn, &   ! (1)
      ! WARNING these are not the fracture nodes but the set of nodes
      !         of the two cells on each side of the fracture
      CellbyFracOwn, &   ! (1)
      FracbyFracOwn, &   ! (1)
      !
      FacebyCellLocal, & ! (1,2)
      FracbyCellLocal, & ! (1,2)
      NodebyCellLocal, & ! (1,2)
      !
      NodebyFaceLocal, & ! (1,2)
      NodebyFractureLocal

   ! Number of Edges by Well
   integer, allocatable, dimension(:), protected :: &
      NbEdgebyWellInjLocal, &
      NbEdgebyWellProdLocal, &
      NbEdgebyMSWellLocal

   ! Num (local) of Nodes by Edge by Well local (own+ghost)
   integer, allocatable, dimension(:, :, :), protected :: &
      NumNodebyEdgebyWellInjLocal, &
      NumNodebyEdgebyWellProdLocal, &
      NumNodebyEdgebyMSWellLocal

   ! 3. X node
   double precision, allocatable, dimension(:, :), target :: &
      XNodeLocal

   integer(c_int), allocatable, dimension(:), target :: &
      NodeFlagsLocal, &
      CellFlagsLocal, &
      FaceFlagsLocal

   integer(c_int8_t), allocatable, dimension(:), target :: &
      CellTypesLocal, &
      FaceTypesLocal

   integer(c_int), allocatable, dimension(:), target, public :: &
      AllDarcyRocktypesLocal
   integer(c_int), dimension(:), pointer, public :: &
      NodeDarcyRocktypesLocal, &
      CellDarcyRocktypesLocal, &
      FracDarcyRocktypesLocal

#ifdef _THERMIQUE_
   integer(c_int), allocatable, dimension(:), target, public :: &
      AllFourierRocktypesLocal
   integer(c_int), dimension(:), pointer, public :: &
      NodeFourierRocktypesLocal, &
      CellFourierRocktypesLocal, &
      FracFourierRocktypesLocal
#endif

   ! 4. IdCell/IdFace/IdNode
   integer, allocatable, dimension(:), protected :: &
      IdCellLocal, &
      IdFaceLocal

   type(Type_IdNode), allocatable, dimension(:), target :: &
      IdNodeLocal

#ifdef _WITH_FREEFLOW_STRUCTURES_
   logical(c_bool), allocatable, dimension(:), target :: &
      IsFreeflowNode    !< Boolean to identify the nodes Freeflow Boundary Condition (atmospheric BC)
   type(TYPE_FFfarfield), allocatable, dimension(:), target :: &
      AtmState    !< contains the atm values at each FF node
#endif

   ! Well and MSWells
   integer, dimension(:), allocatable, target, protected :: &
      NbWellInjLocal_Ncpus, NbWellInjOwn_Ncpus, &
      NbWellProdLocal_Ncpus, NbWellProdOwn_Ncpus, &
      NbMSWellLocal_Ncpus, NbMSWellOwn_Ncpus

   ! Well connectivity in local
   type(CSR), protected :: &
      NodebyWellInjLocal, &
      NodebyWellProdLocal, &
      NodebyMSWellLocal

   type(WellData_type), allocatable, dimension(:), target, public :: &
      DataWellInjLocal, & !< Data of injection well (Radius,...) for local well
      DataWellProdLocal !< Data of production well (Radius,...) for local well

   type(MSWellData_type), allocatable, dimension(:), target, public :: &
      DataMSWellLocal   !< Data of mswell (Radius,...) for local mswell

   ! Data in nodes of wells and mswells
   type(TYPE_CSRDataNodeWell), protected :: &
      NodeDatabyWellInjLocal, &
      NodeDatabyWellProdLocal, &
      NodeDatabyMSWellLocal

   !! The follwoing vectors are used for the strucutre of Jacobian
   type(CSR), protected :: &
      WellInjbyNodeOwn, &     ! numero (local) of well inj connected to this node own
      WellProdbyNodeOwn, &    ! numero (local) of well prod connected to this node own
      MSWellbyNodeOwn         ! numero (local) of mswell  connected to this node own

   ! 5. Frac to Face, Face to Frac
   integer(c_int), allocatable, dimension(:), target :: &
      FracToFaceLocal, &
      FaceToFracLocal

   ! 6. NumNodebyProc, NumFracbyProc
   type(FamilyDOFIdCOC), protected:: &
      NumNodebyProc, &
      NumFracbyProc, &
      NumWellInjbyProc, &
      NumWellProdbyProc

   ! 7. XCellLocal, XFaceLocal
   real(c_double), dimension(:, :), allocatable, public, target :: &
      XCellLocal, &     ! center of cell
      XFaceLocal       ! center of frac

   ! 8. VolCellLocal, SurfFreeFlowLocal, SurfFracLocal
   double precision, dimension(:), allocatable, protected :: &
      VolCellLocal, &      ! vol of cell
      SurfFracLocal        ! surf of frac face
   real(c_double), dimension(:), allocatable, target :: &
      SurfFreeFlowLocal ! area of faces allocated to each freeflow node (size is nb_nodes)

   ! 9. max number of nodes/frac in a cell
   integer, protected :: &
      NbNodeCellMax, &
      NbFracCellMax, &
      NbNodeFaceMax

   ! 10. Porosity
   real(c_double), dimension(:), allocatable, target :: &
      PorositeCellLocal, &
      PorositeFracLocal
   ! permeability
   real(c_double), dimension(:, :, :), allocatable, target :: &
      PermCellLocal
   real(c_double), dimension(:), allocatable, target :: &
      PermFracLocal

#ifdef _THERMIQUE_
   ! Thermal conductivity
   real(c_double), dimension(:, :, :), allocatable, target :: &
      CondThermalCellLocal
   real(c_double), dimension(:), allocatable, target :: &
      CondThermalFracLocal
#endif

#ifdef _THERMIQUE_
   ! Thermal source
   real(c_double), DIMENSION(:), ALLOCATABLE, PUBLIC, target :: CellThermalSourceLocal
   real(c_double), DIMENSION(:), ALLOCATABLE, PUBLIC, target :: FracThermalSourceLocal
#endif

   ! MPI TYPE for DataNodewell: MPI_DATANODEWELL
   integer, private :: MPI_DATANODEWELL

   type SubArraySizes
      integer(c_size_t) :: nodes, fractures, cells
   end type SubArraySizes

   type SubArrayOffsets
      integer(c_size_t) :: nodes, fractures, cells
   end type SubArrayOffsets

   type SubArrayInfo
      type(SubArraySizes) :: nb
      type(SubArrayOffsets) :: offset
   end type SubArrayInfo

   type SubArrayView
      real(c_double), pointer, dimension(:) :: nodes, fractures, cells
   end type SubArrayView

   type DOFFamilyArray
      real(c_double), allocatable, dimension(:) :: values
      real(c_double), pointer, dimension(:) :: nodes, fractures, cells
   end type DOFFamilyArray

   ! FIXME: use parametrized types
   type PhaseDOFFamilyArray
      real(c_double), allocatable, dimension(:, :) :: values
      real(c_double), pointer, dimension(:, :) :: nodes, fractures, cells
   end type PhaseDOFFamilyArray

   ! FIXME: use parametrized types
   type CompPhaseDOFFamilyArray
      real(c_double), allocatable, dimension(:, :, :) :: values
      real(c_double), pointer, dimension(:, :, :) :: nodes, fractures, cells
   end type CompPhaseDOFFamilyArray

   private :: &
      MeshSchema_csrsend, &   ! send csr
      MeshSchema_csrrecv, &   ! recv csr
      MeshSchema_csrdatawellsend, & ! send csrdatawell
      MeshSchema_csrdatawellrecv, & ! recv csrdatawell
      !
      MeshSchema_sendrecv, &
      MeshSchema_NumNodebyEdgebyWellLocal, &
      MeshSchema_XCellLocal, &
      MeshSchema_VolCellLocal, &
      MeshSchema_XFaceLocal, &
      MeshSchema_SurfFracLocal, &
      MeshSchema_NbNodeCellMax, &
      MeshSchema_NbNodeFaceMax, &
      MeshSchema_local_face_surface_from_nodes

   public :: &
      MeshSchema_make, &
      MeshSchema_Free, &
      MeshSchema_triangle_area, &
      MeshSchema_local_face_surface, &
      get_injectors_data, nb_injectors, &
      get_producers_data, nb_producers, &
      MeshSchema_subarrays_sizes, &
      MeshSchema_subarrays_offsets, &
      MeshSchema_subarrays_info, &
      MeshSchema_subarrays_views, &
      MeshSchema_allocate_DOFFamilyArray, MeshSchema_free_DOFFamilyArray, &
      MeshSchema_allocate_PhaseDOFFamilyArray, MeshSchema_free_PhaseDOFFamilyArray, &
      MeshSchema_allocate_CompPhaseDOFFamilyArray, MeshSchema_free_CompPhaseDOFFamilyArray, &
      MeshSchema_part_info

contains

   subroutine retrieve_CSRDataNodeWell(data, wrapper)
      type(TYPE_CSRDataNodeWell), intent(in) :: data
      type(PerforationDataCSR_wrapper), intent(inout) :: wrapper

      if (.not. (associated(data%Pt) .and. associated(data%Num) .and. associated(data%Val))) then
         write (*, *) ""
         write (*, *) " WARNING - Some well pointers are not associated!"
         if (.not. (associated(data%Pt))) &
            write (*, *) "    Pt is not associated!"
         if (.not. (associated(data%Num))) &
            write (*, *) "    Num is not associated!"
         if (.not. (associated(data%Val))) &
            write (*, *) "    Val is not associated!"
         write (*, *) ""
         wrapper%nb_wells = 0
         wrapper%well_offset = c_null_ptr
         wrapper%node_vertex = c_null_ptr
         wrapper%data = c_null_ptr
         return
      endif

      if (data%Nb + 1 /= size(data%Pt)) then
         write (*, *) " WARNING - Inconsistent Pt in CSRDataNodeWell!"
         write (*, *) " Nb =", data%Nb, " Pt size =", size(data%Pt)
         call CommonMPI_abort("Inconsistent Pt in CSRDataNodeWell")
      end if
      if (data%Pt(data%Nb + 1) /= size(data%Num)) &
         call CommonMPI_abort("Inconsistent Num in CSRDataNodeWell")
      if (data%Pt(data%Nb + 1) /= size(data%Val)) &
         call CommonMPI_abort("Inconsistent Val in CSRDataNodeWell")

      wrapper%nb_wells = data%Nb
      if (data%Nb == 0) then
         wrapper%well_offset = c_null_ptr
         wrapper%node_vertex = c_null_ptr
         wrapper%data = c_null_ptr
      else
         wrapper%well_offset = c_loc(data%Pt(1))
         wrapper%node_vertex = c_loc(data%Num(1))
         wrapper%data = c_loc(data%Val(1))
      end if

   end subroutine retrieve_CSRDataNodeWell

   subroutine retrieve_CSRData_producer(wrapper) &
      bind(C, name="retrieve_CSRData_producer")
      type(PerforationDataCSR_wrapper), intent(inout) :: wrapper

      call retrieve_CSRDataNodeWell(NodeDatabyWellProdLocal, wrapper)

   end subroutine retrieve_CSRData_producer

   subroutine retrieve_CSRData_injector(wrapper) &
      bind(C, name="retrieve_CSRData_injector")
      type(PerforationDataCSR_wrapper), intent(inout) :: wrapper

      call retrieve_CSRDataNodeWell(NodeDatabyWellInjLocal, wrapper)

   end subroutine retrieve_CSRData_injector

   subroutine MeshSchema_part_info_by_rank(info, rank) &
      bind(C, name="MeshSchema_part_info_by_rank")
      type(PartInfo), intent(inout) :: info
      integer(c_int), intent(in) :: rank

      if (rank >= Ncpus) then
         call CommonMPI_abort("Index Error : Retrieving PartInfo for a"// &
                              "proc rank greater than number of procs")
      end if

      info%ncpus = Ncpus
      info%rank = rank
      info%nodes%nb_owns = NbNodeOwn_Ncpus(rank + 1)
      info%nodes%nb = NbNodeLocal_Ncpus(rank + 1)
      info%fractures%nb_owns = NbFracOwn_Ncpus(rank + 1)
      info%fractures%nb = NbFracLocal_Ncpus(rank + 1)
      info%injectors%nb_owns = NbWellInjOwn_Ncpus(rank + 1)
      info%injectors%nb = NbWellInjLocal_Ncpus(rank + 1)
      info%producers%nb_owns = NbWellProdOwn_Ncpus(rank + 1)
      info%producers%nb = NbWellProdLocal_Ncpus(rank + 1)

   end subroutine MeshSchema_part_info_by_rank

   subroutine MeshSchema_part_info(info) &
      bind(C, name="MeshSchema_part_info")
      type(PartInfo), intent(out) :: info

      info%ncpus = Ncpus
      info%rank = commRank
      info%nodes%nb_owns = NbNodeOwn_Ncpus(commRank + 1)
      info%nodes%nb = NbNodeLocal_Ncpus(commRank + 1)
      info%fractures%nb_owns = NbFracOwn_Ncpus(commRank + 1)
      info%fractures%nb = NbFracLocal_Ncpus(commRank + 1)
      info%injectors%nb_owns = NbWellInjOwn_Ncpus(commRank + 1)
      info%injectors%nb = NbWellInjLocal_Ncpus(commRank + 1)
      info%producers%nb_owns = NbWellProdOwn_Ncpus(commRank + 1)
      info%producers%nb = NbWellProdLocal_Ncpus(commRank + 1)

   end subroutine MeshSchema_part_info

   subroutine MeshSchema_subarrays_sizes(sizes)
      type(SubArraySizes), intent(out) :: sizes

      sizes%nodes = NbNodeLocal_Ncpus(commRank + 1)
      sizes%fractures = NbFracLocal_Ncpus(commRank + 1)
      sizes%cells = NbCellLocal_Ncpus(commRank + 1)

   end subroutine MeshSchema_subarrays_sizes

   subroutine MeshSchema_subarrays_offsets(offsets)
      type(SubArrayOffsets), intent(out) :: offsets
      type(SubArraySizes) :: nb

      call MeshSchema_subarrays_compute_info(nb, offsets)

   end subroutine MeshSchema_subarrays_offsets

   subroutine MeshSchema_subarrays_compute_info(nb, offsets)
      type(SubArraySizes), intent(out) :: nb
      type(SubArrayOffsets), intent(out) :: offsets

      call MeshSchema_subarrays_sizes(nb)
      offsets%nodes = 1
      offsets%fractures = offsets%nodes + nb%nodes
      offsets%cells = offsets%fractures + nb%fractures

   end subroutine MeshSchema_subarrays_compute_info

   subroutine MeshSchema_subarrays_info(info)
      type(SubArrayInfo), intent(out) :: info

      call MeshSchema_subarrays_compute_info(info%nb, info%offset)

   end subroutine MeshSchema_subarrays_info

   subroutine MeshSchema_subarrays_views(a, views)
      real(c_double), allocatable, dimension(:), target, intent(in) :: a
      type(SubArrayView), intent(out) :: views

      type(SubArraySizes) :: nb
      type(SubArrayOffsets) :: offset

      call MeshSchema_subarrays_compute_info(nb, offset)

      if (.not. allocated(a)) &
         call CommonMPI_abort('MeshSchema_subarrays_views: unallocated array')
      if (size(a) /= nb%nodes + nb%fractures + nb%cells) &
         call CommonMPI_abort('MeshSchema_subarrays_views: unconsistent sizes')

      views%nodes => a(offset%nodes:offset%nodes - 1 + nb%nodes)
      views%fractures => a(offset%fractures:offset%fractures - 1 + nb%fractures)
      views%cells => a(offset%cells:offset%cells - 1 + nb%cells)

   end subroutine MeshSchema_subarrays_views

   subroutine MeshSchema_retrieve_local_cell_permeability(buffer) &
      bind(C, name="retrieve_cell_permeability")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(PermCellLocal))) &
         call CommonMPI_abort("PermCellLocal is not allocated")
      buffer%p = c_loc(PermCellLocal(1, 1, 1))
      buffer%n = size(PermCellLocal, 3, c_size_t)

   end subroutine MeshSchema_retrieve_local_cell_permeability

   subroutine MeshSchema_retrieve_local_fracture_permeability(buffer) &
      bind(C, name="retrieve_fracture_permeability")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(PermFracLocal))) &
         call CommonMPI_abort("PermFracLocal is not allocated")
      buffer%p = c_loc(PermFracLocal(1))
      buffer%n = size(PermFracLocal, 1, c_size_t)

   end subroutine MeshSchema_retrieve_local_fracture_permeability

#ifdef _THERMIQUE_

   subroutine MeshSchema_retrieve_local_cell_thermal_conductivity(buffer) &
      bind(C, name="retrieve_cell_thermal_conductivity")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(CondThermalCellLocal))) &
         call CommonMPI_abort("CondThermalCellLocal is not allocated")
      buffer%p = c_loc(CondThermalCellLocal(1, 1, 1))
      buffer%n = size(CondThermalCellLocal, 3, c_size_t)

   end subroutine MeshSchema_retrieve_local_cell_thermal_conductivity

   subroutine MeshSchema_retrieve_local_fracture_thermal_conductivity(buffer) &
      bind(C, name="retrieve_fracture_thermal_conductivity")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(CondThermalFracLocal))) &
         call CommonMPI_abort("CondThermalFracLocal is not allocated")
      buffer%p = c_loc(CondThermalFracLocal(1))
      buffer%n = size(CondThermalFracLocal, 1, c_size_t)

   end subroutine MeshSchema_retrieve_local_fracture_thermal_conductivity

#endif

   subroutine MeshSchema_retrieve_local_cell_porosity(buffer) &
      bind(C, name="retrieve_cell_porosity")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(PorositeCellLocal))) &
         call CommonMPI_abort("PorositeCellLocal is not allocated")
      buffer%p = c_loc(PorositeCellLocal(1))
      buffer%n = size(PorositeCellLocal, 1, c_size_t)

   end subroutine MeshSchema_retrieve_local_cell_porosity

   subroutine MeshSchema_retrieve_local_fracture_porosity(buffer) &
      bind(C, name="retrieve_fracture_porosity")
      type(cpp_array_wrapper), intent(inout) :: buffer

      if (.not. (allocated(PorositeFracLocal))) &
         call CommonMPI_abort("PorositeFracLocal is not allocated")
      buffer%p = c_loc(PorositeFracLocal(1))
      buffer%n = size(PorositeFracLocal, 1, c_size_t)

   end subroutine MeshSchema_retrieve_local_fracture_porosity

   subroutine MeshSchema_allocate_DOFFamilyArray(a)
      type(DOFFamilyArray), target, intent(inout) :: a

      type(SubArraySizes) :: nb
      type(SubArrayOffsets) :: offset

      call MeshSchema_subarrays_compute_info(nb, offset)

      if (allocated(a%values)) then
         if (size(a%values) /= (nb%nodes + nb%fractures + nb%cells)) &
            call CommonMPI_abort('MeshSchema_allocate_DOFFamilyArray: already allocated with incompatible size.')
      else
         allocate (a%values(nb%nodes + nb%fractures + nb%cells))
      endif

      nullify (a%nodes, a%fractures, a%cells)
      a%nodes => a%values(offset%nodes:offset%nodes - 1 + nb%nodes)
      a%fractures => a%values(offset%fractures:offset%fractures - 1 + nb%fractures)
      a%cells => a%values(offset%cells:offset%cells - 1 + nb%cells)

   end subroutine MeshSchema_allocate_DOFFamilyArray

   subroutine MeshSchema_free_DOFFamilyArray(a)
      type(DOFFamilyArray), intent(inout) :: a

      if (allocated(a%values)) deallocate (a%values)
      nullify (a%nodes, a%fractures, a%cells)

   end subroutine MeshSchema_free_DOFFamilyArray

   subroutine MeshSchema_allocate_PhaseDOFFamilyArray(a)
      type(PhaseDOFFamilyArray), target, intent(inout) :: a

      type(SubArraySizes) :: nb
      type(SubArrayOffsets) :: offset

      if (allocated(a%values)) &
         call CommonMPI_abort('MeshSchema_allocate_DOFFamilyArray: already allocated')

      call MeshSchema_subarrays_compute_info(nb, offset)

      allocate (a%values(NbPhase, nb%nodes + nb%fractures + nb%cells))
      nullify (a%nodes, a%fractures, a%cells)
      a%nodes => a%values(:, offset%nodes:offset%nodes - 1 + nb%nodes)
      a%fractures => a%values(:, offset%fractures:offset%fractures - 1 + nb%fractures)
      a%cells => a%values(:, offset%cells:offset%cells - 1 + nb%cells)

   end subroutine MeshSchema_allocate_PhaseDOFFamilyArray

   subroutine MeshSchema_free_PhaseDOFFamilyArray(a)
      type(PhaseDOFFamilyArray), intent(inout) :: a

      if (allocated(a%values)) deallocate (a%values)
      nullify (a%nodes, a%fractures, a%cells)

   end subroutine MeshSchema_free_PhaseDOFFamilyArray

   subroutine MeshSchema_allocate_CompPhaseDOFFamilyArray(a)
      type(CompPhaseDOFFamilyArray), target, intent(inout) :: a

      type(SubArraySizes) :: nb
      type(SubArrayOffsets) :: offset

      if (allocated(a%values)) &
         call CommonMPI_abort('MeshSchema_allocate_DOFFamilyArray: already allocated')

      call MeshSchema_subarrays_compute_info(nb, offset)

      allocate (a%values(NbComp, NbPhase, nb%nodes + nb%fractures + nb%cells))
      nullify (a%nodes, a%fractures, a%cells)
      a%nodes => a%values(:, :, offset%nodes:offset%nodes - 1 + nb%nodes)
      a%fractures => a%values(:, :, offset%fractures:offset%fractures - 1 + nb%fractures)
      a%cells => a%values(:, :, offset%cells:offset%cells - 1 + nb%cells)

   end subroutine MeshSchema_allocate_CompPhaseDOFFamilyArray

   subroutine MeshSchema_free_CompPhaseDOFFamilyArray(a)
      type(CompPhaseDOFFamilyArray), intent(inout) :: a

      if (allocated(a%values)) deallocate (a%values)
      nullify (a%nodes, a%fractures, a%cells)

   end subroutine MeshSchema_free_CompPhaseDOFFamilyArray

   function get_injectors_data() result(p) &
      bind(C, name="get_injectors_data")

      type(c_ptr) :: p

      if (allocated(DataWellInjLocal) .and. size(DataWellInjLocal, 1) > 0) then
         p = c_loc(DataWellInjLocal(1))
      else
         p = c_null_ptr
      end if

   end function get_injectors_data

   function nb_injectors() result(n) &
      bind(C, name="nb_injectors")

      integer(c_size_t) :: n

      if (allocated(DataWellInjLocal)) then
         if (.not. allocated(NbWellInjLocal_Ncpus)) &
            call CommonMPI_abort("NbWellInjLocal_Ncpus not allocated.")
         if (NbWellInjLocal_Ncpus(commRank + 1) /= size(DataWellInjLocal, 1)) &
            call CommonMPI_abort("NbWellInjLocal_Ncpus inconsistency")
         n = size(DataWellInjLocal, 1)
      else
         n = 0
      end if

   end function nb_injectors

   function get_producers_data() result(p) &
      bind(C, name="get_producers_data")

      type(c_ptr) :: p

      if (allocated(DataWellProdLocal) .and. size(DataWellProdLocal, 1) > 0) then
         p = c_loc(DataWellProdLocal(1))
      else
         p = c_null_ptr
      end if

   end function get_producers_data

   function nb_producers() result(n) &
      bind(C, name="nb_producers")

      integer(c_size_t) :: n

      if (allocated(DataWellProdLocal)) then
         if (.not. allocated(NbWellProdLocal_Ncpus)) &
            call CommonMPI_abort("NbWellProdLocal_Ncpus not allocated.")
         if (NbWellProdLocal_Ncpus(commRank + 1) /= size(DataWellProdLocal, 1)) &
            call CommonMPI_abort("NbWellProdLocal_Ncpus inconsistency")
         n = size(DataWellProdLocal, 1)
      else
         n = 0
      end if

   end function nb_producers

   pure function number_of_nodes() result(n) &
      bind(C, name="number_of_nodes")
      integer(c_size_t) :: n
      n = NbNodeLocal_Ncpus(commRank + 1)
   end function number_of_nodes

   pure function number_of_own_nodes() result(n) &
      bind(C, name="number_of_own_nodes")
      integer(c_size_t) :: n
      n = NbNodeOwn_Ncpus(commRank + 1)
   end function number_of_own_nodes

   pure function number_of_cells() result(n) &
      bind(C, name="number_of_cells")
      integer(c_size_t) :: n
      n = NbCellLocal_Ncpus(commRank + 1)
   end function number_of_cells

   pure function number_of_own_cells() result(n) &
      bind(C, name="number_of_own_cells")
      integer(c_size_t) :: n
      n = NbCellOwn_Ncpus(commRank + 1)
   end function number_of_own_cells

   pure function number_of_faces() result(n) &
      bind(C, name="number_of_faces")
      integer(c_size_t) :: n
      n = NbFaceLocal_Ncpus(commRank + 1)
   end function number_of_faces

   pure function number_of_own_faces() result(n) &
      bind(C, name="number_of_own_faces")
      integer(c_size_t) :: n
      n = NbFaceOwn_Ncpus(commRank + 1)
   end function number_of_own_faces

   pure function number_of_fractures() result(n) &
      bind(C, name="number_of_fractures")
      integer(c_size_t) :: n
      n = NbFracLocal_Ncpus(commRank + 1)
   end function number_of_fractures

   pure function number_of_own_fractures() result(n) &
      bind(C, name="number_of_own_fractures")
      integer(c_size_t) :: n
      n = NbFracOwn_Ncpus(commRank + 1)
   end function number_of_own_fractures

   pure function number_of_own_injectors() result(n) &
      bind(C, name="number_of_own_injectors")
      integer(c_size_t) :: n
      n = NbWellInjOwn_Ncpus(commRank + 1)
   end function number_of_own_injectors

   pure function number_of_own_producers() result(n) &
      bind(C, name="number_of_own_producers")
      integer(c_size_t) :: n
      n = NbWellProdOwn_Ncpus(commRank + 1)
   end function number_of_own_producers

   subroutine MeshSchema_make

      call MeshSchema_sendrecv
      call MeshSchema_collect_fracture_nodes

      ! List of Edge (local number of nodes) by local Well (own+ghost)
      call MeshSchema_NumNodebyEdgebyWellLocal(NodeDatabyWellInjLocal, &
                                               NbEdgebyWellInjLocal, &
                                               NumNodebyEdgebyWellInjLocal)

      call MeshSchema_NumNodebyEdgebyWellLocal(NodeDatabyWellProdLocal, &
                                               NbEdgebyWellProdLocal, &
                                               NumNodebyEdgebyWellProdLocal)

      ! List of Edge (local number of nodes) by local MSWell (own+ghost)
      call MeshSchema_NumNodebyEdgebyWellLocal(NodeDatabyMSWellLocal, &
                                               NbEdgebyMSWellLocal, &
                                               NumNodebyEdgebyMSWellLocal)

      ! XCellLocal and XFaceLocal
      call MeshSchema_XCellLocal
      call MeshSchema_XFaceLocal

      ! VolCellLocal and SurfFraclocal
      call MeshSchema_VolCellLocal
      call MeshSchema_SurfFracLocal

      call MeshSchema_NbNodeCellMax ! max nb of nodes in a cell
      call MeshSchema_NbFracCellMax ! max nb of frac in a cell
      call MeshSchema_NbNodeFaceMax ! max nb of nodes in a frac

   end subroutine MeshSchema_make

   ! send/receive (commRank>=1)
   ! or copy (commRank=0)
   subroutine MeshSchema_sendrecv

      integer :: dest, Ierr, i, j, Nb
      integer stat(MPI_STATUS_SIZE)

      integer :: blen(1), offsets(1), oldtypes(1), MPI_IDNODE
      integer :: blocklen(5), arraytype(5)
      integer(kind=MPI_ADDRESS_KIND) ::disp(5)

      integer :: MPI_WELLDATA_ID, MPI_MSWELLDATA_ID, nb_nodes, nb_fracs, nb_cells

      ! ************************************* !

      ! Send Nb*

      if (commRank == 0) then

         do dest = 1, Ncpus - 1

            ! Nb*ResS_Ncpus
            call MPI_Send(NbCellResS_Ncpus, Ncpus, MPI_INTEGER, dest, 11, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbFaceResS_Ncpus, Ncpus, MPI_INTEGER, dest, 12, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbNodeResS_Ncpus, Ncpus, MPI_INTEGER, dest, 13, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbFracResS_Ncpus, Ncpus, MPI_INTEGER, dest, 14, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbWellInjResS_Ncpus, Ncpus, MPI_INTEGER, dest, 141, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbWellProdResS_Ncpus, Ncpus, MPI_INTEGER, dest, 142, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbMSWellResS_Ncpus, Ncpus, MPI_INTEGER, dest, 143, ComPASS_COMM_WORLD, Ierr)

            ! Nb*OwnS_Ncpus
            call MPI_Send(NbCellOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 15, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbFaceOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 16, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbNodeOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 17, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbFracOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 18, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbWellInjOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 181, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbWellProdOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 182, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(NbMSWellOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 183, ComPASS_COMM_WORLD, Ierr)

         end do

         allocate (NbCellLocal_Ncpus(Ncpus))
         allocate (NbFaceLocal_Ncpus(Ncpus))
         allocate (NbNodeLocal_Ncpus(Ncpus))
         allocate (NbFracLocal_Ncpus(Ncpus))
         allocate (NbWellInjLocal_Ncpus(Ncpus))
         allocate (NbWellProdLocal_Ncpus(Ncpus))
         allocate (NbMSWellLocal_Ncpus(Ncpus))

         allocate (NbCellOwn_Ncpus(Ncpus))
         allocate (NbFaceOwn_Ncpus(Ncpus))
         allocate (NbNodeOwn_Ncpus(Ncpus))
         allocate (NbFracOwn_Ncpus(Ncpus))
         allocate (NbWellInjOwn_Ncpus(Ncpus))
         allocate (NbWellProdOwn_Ncpus(Ncpus))
         allocate (NbMSWellOwn_Ncpus(Ncpus))

         NbCellLocal_Ncpus(:) = NbCellResS_Ncpus(:)
         NbFaceLocal_Ncpus(:) = NbFaceResS_Ncpus(:)
         NbNodeLocal_Ncpus(:) = NbNodeResS_Ncpus(:)
         NbFracLocal_Ncpus(:) = NbFracResS_Ncpus(:)
         NbWellInjLocal_Ncpus(:) = NbWellInjResS_Ncpus(:)
         NbWellProdLocal_Ncpus(:) = NbWellProdResS_Ncpus(:)
         NbMSWellLocal_Ncpus(:) = NbMSWellResS_Ncpus(:)

         NbCellOwn_Ncpus(:) = NbCellOwnS_Ncpus(:)
         NbFaceOwn_Ncpus(:) = NbFaceOwnS_Ncpus(:)
         NbNodeOwn_Ncpus(:) = NbNodeOwnS_Ncpus(:)
         NbFracOwn_Ncpus(:) = NbFracOwnS_Ncpus(:)
         NbWellInjOwn_Ncpus(:) = NbWellInjOwnS_Ncpus(:)
         NbWellProdOwn_Ncpus(:) = NbWellProdOwnS_Ncpus(:)
         NbMSWellOwn_Ncpus(:) = NbMSWellOwnS_Ncpus(:)

      else
         allocate (NbCellLocal_Ncpus(Ncpus))
         allocate (NbFaceLocal_Ncpus(Ncpus))
         allocate (NbNodeLocal_Ncpus(Ncpus))
         allocate (NbFracLocal_Ncpus(Ncpus))
         allocate (NbWellInjLocal_Ncpus(Ncpus))
         allocate (NbWellProdLocal_Ncpus(Ncpus))
         allocate (NbMSWellLocal_Ncpus(Ncpus))

         allocate (NbCellOwn_Ncpus(Ncpus))
         allocate (NbFaceOwn_Ncpus(Ncpus))
         allocate (NbNodeOwn_Ncpus(Ncpus))
         allocate (NbFracOwn_Ncpus(Ncpus))
         allocate (NbWellInjOwn_Ncpus(Ncpus))
         allocate (NbWellProdOwn_Ncpus(Ncpus))
         allocate (NbMSWellOwn_Ncpus(Ncpus))

         ! Nb*ResS_Ncpus
         call MPI_recv(NbCellLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbFaceLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbNodeLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbFracLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbWellInjLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 141, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbWellProdLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 142, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbMSWellLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 143, ComPASS_COMM_WORLD, stat, Ierr)

         ! Nb*OwnS_Ncpus
         call MPI_Recv(NbCellOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 15, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbFaceOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 16, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbNodeOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 17, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbFracOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 18, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbWellInjOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 181, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbWellProdOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 182, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_Recv(NbMSWellOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 183, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      !write(*,*) 'on proc', commRank, ':', NbWellInjOwn_Ncpus, 'injectors'
      !write(*,*) 'on proc', commRank, ':', NbWellProdOwn_Ncpus, 'producers'

      if (commRank == 0) then
         deallocate (NbCellResS_Ncpus)
         deallocate (NbFaceResS_Ncpus)
         deallocate (NbNodeResS_Ncpus)
         deallocate (NbFracResS_Ncpus)
         deallocate (NbWellInjResS_Ncpus)
         deallocate (NbWellProdResS_Ncpus)
         deallocate (NbMSWellResS_Ncpus)
         deallocate (NbCellOwnS_Ncpus)
         deallocate (NbFaceOwnS_Ncpus)
         deallocate (NbNodeOwnS_Ncpus)
         deallocate (NbFracOwnS_Ncpus)
         deallocate (NbWellInjOwnS_Ncpus)
         deallocate (NbWellProdOwnS_Ncpus)
         deallocate (NbMSWellOwnS_Ncpus)
      end if

      ! Send mesh size
      if (commRank == 0) then

         do dest = 1, Ncpus - 1

            ! meshSize_*max, meshSize_*min
            call MPI_Send(meshSizeS_xmax, 1, MPI_DOUBLE, dest, 11, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(meshSizeS_xmin, 1, MPI_DOUBLE, dest, 12, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(meshSizeS_ymax, 1, MPI_DOUBLE, dest, 13, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(meshSizeS_ymin, 1, MPI_DOUBLE, dest, 14, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(meshSizeS_zmax, 1, MPI_DOUBLE, dest, 15, ComPASS_COMM_WORLD, Ierr)
            call MPI_Send(meshSizeS_zmin, 1, MPI_DOUBLE, dest, 16, ComPASS_COMM_WORLD, Ierr)
         end do

         meshSize_xmax = meshSizeS_xmax
         meshSize_xmin = meshSizeS_xmin
         meshSize_ymax = meshSizeS_ymax
         meshSize_ymin = meshSizeS_ymin
         meshSize_zmax = meshSizeS_zmax
         meshSize_zmin = meshSizeS_zmin
      else

         call MPI_recv(meshSize_xmax, 1, MPI_DOUBLE, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_recv(meshSize_xmin, 1, MPI_DOUBLE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_recv(meshSize_ymax, 1, MPI_DOUBLE, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_recv(meshSize_ymin, 1, MPI_DOUBLE, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_recv(meshSize_zmax, 1, MPI_DOUBLE, 0, 15, ComPASS_COMM_WORLD, stat, Ierr)
         call MPI_recv(meshSize_zmin, 1, MPI_DOUBLE, 0, 16, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! **************************************** !

      ! Send XNodeRes_Ncpus to XNodeLocal
      if (commRank == 0) then

         do i = 1, Ncpus - 1
            call MPI_Send(XNodeRes_Ncpus(i + 1)%Array2d, NbNodeLocal_Ncpus(i + 1)*3, MPI_DOUBLE, i, 21, ComPASS_COMM_WORLD, Ierr)
         end do

         allocate (XNodeLocal(3, NbNodeLocal_Ncpus(1)))
         do j = 1, NbNodeLocal_Ncpus(1)
            XNodeLocal(1, j) = XNodeRes_Ncpus(1)%Array2d(1, j)
            XNodeLocal(2, j) = XNodeRes_Ncpus(1)%Array2d(2, j)
            XNodeLocal(3, j) = XNodeRes_Ncpus(1)%Array2d(3, j)
         end do

      else
         allocate (XNodeLocal(3, NbNodeLocal_Ncpus(commRank + 1)))
         call MPI_Recv(XNodeLocal, NbNodeLocal_Ncpus(commRank + 1)*3, MPI_DOUBLE, 0, 21, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (XNodeRes_Ncpus(i)%Array2d)
         end do
         deallocate (XNodeRes_Ncpus)
      end if

      ! Send flags

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send(NodeFlags_Ncpus(i + 1)%Val, NbNodeLocal_Ncpus(i + 1), MPI_INTEGER, i, 22, ComPASS_COMM_WORLD, Ierr)
         end do
         allocate (NodeFlagsLocal(NbNodeLocal_Ncpus(1)))
         NodeFlagsLocal = NodeFlags_Ncpus(1)%Val
      else
         allocate (NodeFlagsLocal(NbNodeLocal_Ncpus(commRank + 1)))
         call MPI_Recv(NodeFlagsLocal, NbNodeLocal_Ncpus(commRank + 1), MPI_INTEGER, 0, 22, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send(CellFlags_Ncpus(i + 1)%Val, NbCellLocal_Ncpus(i + 1), MPI_INTEGER, i, 23, ComPASS_COMM_WORLD, Ierr)
         end do
         allocate (CellFlagsLocal(NbCellLocal_Ncpus(1)))
         CellFlagsLocal = CellFlags_Ncpus(1)%Val
      else
         allocate (CellFlagsLocal(NbCellLocal_Ncpus(commRank + 1)))
         call MPI_Recv(CellFlagsLocal, NbCellLocal_Ncpus(commRank + 1), MPI_INTEGER, 0, 23, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send(FaceFlags_Ncpus(i + 1)%Val, NbFaceLocal_Ncpus(i + 1), MPI_INTEGER, i, 24, ComPASS_COMM_WORLD, Ierr)
         end do
         allocate (FaceFlagsLocal(NbFaceLocal_Ncpus(1)))
         FaceFlagsLocal = FaceFlags_Ncpus(1)%Val
      else
         allocate (FaceFlagsLocal(NbFaceLocal_Ncpus(commRank + 1)))
         call MPI_Recv(FaceFlagsLocal, NbFaceLocal_Ncpus(commRank + 1), MPI_INTEGER, 0, 24, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send(CellTypes_Ncpus(i + 1)%Val, NbCellLocal_Ncpus(i + 1), MPI_INTEGER1, i, 25, ComPASS_COMM_WORLD, Ierr)
         end do
         allocate (CellTypesLocal(NbCellLocal_Ncpus(1)))
         CellTypesLocal = CellTypes_Ncpus(1)%Val
      else
         allocate (CellTypesLocal(NbCellLocal_Ncpus(commRank + 1)))
         call MPI_Recv(CellTypesLocal, NbCellLocal_Ncpus(commRank + 1), MPI_INTEGER1, 0, 25, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send(FaceTypes_Ncpus(i + 1)%Val, NbFaceLocal_Ncpus(i + 1), MPI_INTEGER1, i, 26, ComPASS_COMM_WORLD, Ierr)
         end do
         allocate (FaceTypesLocal(NbFaceLocal_Ncpus(1)))
         FaceTypesLocal = FaceTypes_Ncpus(1)%Val
      else
         allocate (FaceTypesLocal(NbFaceLocal_Ncpus(commRank + 1)))
         call MPI_Recv(FaceTypesLocal, NbFaceLocal_Ncpus(commRank + 1), MPI_INTEGER1, 0, 26, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (NodeFlags_Ncpus(i)%Val)
            deallocate (CellFlags_Ncpus(i)%Val)
            deallocate (FaceFlags_Ncpus(i)%Val)
            deallocate (CellTypes_Ncpus(i)%Val)
            deallocate (FaceTypes_Ncpus(i)%Val)
         end do
         deallocate (NodeFlags_Ncpus)
         deallocate (CellFlags_Ncpus)
         deallocate (FaceFlags_Ncpus)
         deallocate (FaceTypes_Ncpus)
      end if

      ! Send Rocktype
      nb_nodes = NbNodeLocal_Ncpus(commRank + 1)
      nb_cells = NbCellLocal_Ncpus(commRank + 1)
      nb_fracs = NbFracLocal_Ncpus(commRank + 1)
      allocate (AllDarcyRocktypesLocal(nb_nodes + nb_fracs + nb_cells))
      NodeDarcyRocktypesLocal => AllDarcyRocktypesLocal(1:nb_nodes)
      FracDarcyRocktypesLocal => AllDarcyRocktypesLocal(nb_nodes + 1:nb_nodes + nb_fracs)
      CellDarcyRocktypesLocal => AllDarcyRocktypesLocal(nb_nodes + nb_fracs + 1:nb_nodes + nb_fracs + nb_cells)

#ifdef _THERMIQUE_
      allocate (AllFourierRocktypesLocal(nb_nodes + nb_fracs + nb_cells))
      NodeFourierRocktypesLocal => AllFourierRocktypesLocal(1:nb_nodes)
      FracFourierRocktypesLocal => AllFourierRocktypesLocal(nb_nodes + 1:nb_nodes + nb_fracs)
      CellFourierRocktypesLocal => AllFourierRocktypesLocal(nb_nodes + nb_fracs + 1:nb_nodes + nb_fracs + nb_cells)
#endif

      ! Nodes
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               NodeRocktype_Ncpus(i + 1)%Array2d(1, :), NbNodeLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 25, ComPASS_COMM_WORLD, Ierr)
         end do
         NodeDarcyRocktypesLocal = NodeRocktype_Ncpus(1)%Array2d(1, :)
      else
         call MPI_Recv( &
            NodeDarcyRocktypesLocal, nb_nodes, &
            MPI_INTEGER, 0, 25, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! Fractures
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               FracRocktype_Ncpus(i + 1)%Array2d(1, :), NbFracLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 25, ComPASS_COMM_WORLD, Ierr)
         end do
         FracDarcyRocktypesLocal = FracRocktype_Ncpus(1)%Array2d(1, :)
      else
         call MPI_Recv( &
            FracDarcyRocktypesLocal, nb_fracs, &
            MPI_INTEGER, 0, 25, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! Cells
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               CellRocktype_Ncpus(i + 1)%Array2d(1, :), NbCellLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 25, ComPASS_COMM_WORLD, Ierr)
         end do
         CellDarcyRocktypesLocal = CellRocktype_Ncpus(1)%Array2d(1, :)
      else
         call MPI_Recv( &
            CellDarcyRocktypesLocal, nb_cells, &
            MPI_INTEGER, 0, 25, ComPASS_COMM_WORLD, stat, Ierr)
      end if

#ifdef _THERMIQUE_

      ! Nodes
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               NodeRocktype_Ncpus(i + 1)%Array2d(2, :), NbNodeLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 28, ComPASS_COMM_WORLD, Ierr)
         end do
         NodeFourierRocktypesLocal = NodeRocktype_Ncpus(1)%Array2d(2, :)
      else
         call MPI_Recv( &
            NodeFourierRocktypesLocal, nb_nodes, &
            MPI_INTEGER, 0, 28, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! Fractures
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               FracRocktype_Ncpus(i + 1)%Array2d(2, :), NbFracLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 29, ComPASS_COMM_WORLD, Ierr)
         end do
         FracFourierRocktypesLocal = FracRocktype_Ncpus(1)%Array2d(2, :)
      else
         call MPI_Recv( &
            FracFourierRocktypesLocal, nb_fracs, &
            MPI_INTEGER, 0, 29, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! Cells
      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MPI_Send( &
               CellRocktype_Ncpus(i + 1)%Array2d(2, :), NbCellLocal_Ncpus(i + 1), &
               MPI_INTEGER, i, 30, ComPASS_COMM_WORLD, Ierr)
         end do
         CellFourierRocktypesLocal = CellRocktype_Ncpus(1)%Array2d(2, :)
      else
         call MPI_Recv( &
            CellFourierRocktypesLocal, nb_cells, &
            MPI_INTEGER, 0, 30, ComPASS_COMM_WORLD, stat, Ierr)
      end if

#endif

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (NodeRocktype_Ncpus(i)%Array2d)
            deallocate (CellRocktype_Ncpus(i)%Array2d)
            deallocate (FracRocktype_Ncpus(i)%Array2d)
         end do
         deallocate (NodeRocktype_Ncpus)
         deallocate (CellRocktype_Ncpus)
         deallocate (FracRocktype_Ncpus)
      end if

      ! ************************************** !

      ! Send mesh and connectivities
      if (commRank == 0) then

         ! for proc >=1, send
         do i = 1, Ncpus - 1

            ! NodebyNodeOwn
            call MeshSchema_csrsend(NodebyNodeOwn_Ncpus(i + 1), i, 100, VALSIZE_ZERO)
            ! FracbyNodeOwn
            call MeshSchema_csrsend(FracbyNodeOwn_Ncpus(i + 1), i, 200, VALSIZE_ZERO)
            ! CellbyNodeOwn
            call MeshSchema_csrsend(CellbyNodeOwn_Ncpus(i + 1), i, 300, VALSIZE_ZERO)

            ! NodebyFracOwn
            call MeshSchema_csrsend(NodebyFracOwn_Ncpus(i + 1), i, 400, VALSIZE_ZERO)
            ! CellbyFracOwn
            call MeshSchema_csrsend(CellbyFracOwn_Ncpus(i + 1), i, 500, VALSIZE_ZERO)
            ! FracbyFracOwn
            call MeshSchema_csrsend(FracbyFracOwn_Ncpus(i + 1), i, 600, VALSIZE_ZERO)
            ! WellInjbyNodeOwn
            call MeshSchema_csrsend(WellInjbyNodeOwn_Ncpus(i + 1), i, 610, VALSIZE_ZERO)
            ! WellProdbyNodeOwn
            call MeshSchema_csrsend(WellProdbyNodeOwn_Ncpus(i + 1), i, 620, VALSIZE_ZERO)
            ! MSWells
            call MeshSchema_csrsend(MSWellbyNodeOwn_Ncpus(i + 1), i, 630, VALSIZE_ZERO)

            ! FacebyCellLocal
            call MeshSchema_csrsend(FacebyCellLocal_Ncpus(i + 1), i, 700, VALSIZE_ZERO)
            ! FracbyCellLocal
            call MeshSchema_csrsend(FracbyCellLocal_Ncpus(i + 1), i, 800, VALSIZE_ZERO)
            ! NodebyCellLocal
            call MeshSchema_csrsend(NodebyCellLocal_Ncpus(i + 1), i, 900, VALSIZE_ZERO)

            ! NodebyFaceLocal
            call MeshSchema_csrsend(NodebyFaceLocal_Ncpus(i + 1), i, 1000, VALSIZE_ZERO)

            ! NodebyWellInjLocal
            call MeshSchema_csrsend(NodebyWellInjLocal_Ncpus(i + 1), i, 1100, VALSIZE_ZERO)
            ! NodebyWellProdLocal
            call MeshSchema_csrsend(NodebyWellProdLocal_Ncpus(i + 1), i, 1200, VALSIZE_ZERO)
            ! NodebyMSWellProdLocal
            call MeshSchema_csrsend(NodebyMSWellLocal_Ncpus(i + 1), i, 1300, VALSIZE_ZERO)

         end do

         ! for proc 0, copy
         call CommonType_csrcopy(NodebyNodeOwn_Ncpus(1), NodebyNodeOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(FracbyNodeOwn_Ncpus(1), FracbyNodeOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(CellbyNodeOwn_Ncpus(1), CellbyNodeOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(WellInjbyNodeOwn_Ncpus(1), WellInjbyNodeOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(WellProdbyNodeOwn_Ncpus(1), WellProdbyNodeOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(MSWellbyNodeOwn_Ncpus(1), MSWellbyNodeOwn, VALSIZE_ZERO)

         call CommonType_csrcopy(NodebyFracOwn_Ncpus(1), NodebyFracOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(CellbyFracOwn_Ncpus(1), CellbyFracOwn, VALSIZE_ZERO)
         call CommonType_csrcopy(FracbyFracOwn_Ncpus(1), FracbyFracOwn, VALSIZE_ZERO)

         call CommonType_csrcopy(FacebyCellLocal_Ncpus(1), FacebyCellLocal, VALSIZE_ZERO)
         call CommonType_csrcopy(FracbyCellLocal_Ncpus(1), FracbyCellLocal, VALSIZE_ZERO)
         call CommonType_csrcopy(NodebyCellLocal_Ncpus(1), NodebyCellLocal, VALSIZE_ZERO)
         call CommonType_csrcopy(NodebyFaceLocal_Ncpus(1), NodebyFaceLocal, VALSIZE_ZERO)

         call CommonType_csrcopy(NodebyWellInjLocal_Ncpus(1), NodebyWellInjLocal, VALSIZE_ZERO)
         call CommonType_csrcopy(NodebyWellProdLocal_Ncpus(1), NodebyWellProdLocal, VALSIZE_ZERO)
         call CommonType_csrcopy(NodebyMSWellLocal_Ncpus(1), NodebyMSWellLocal, VALSIZE_ZERO)

      else

         ! NodebyNodeOwn
         call MeshSchema_csrrecv(NodebyNodeOwn, 0, 100, VALSIZE_ZERO)
         ! FracbyNodeOwn
         call MeshSchema_csrrecv(FracbyNodeOwn, 0, 200, VALSIZE_ZERO)
         ! CellbyNodeOwn
         call MeshSchema_csrrecv(CellbyNodeOwn, 0, 300, VALSIZE_ZERO)

         ! NodebyFracOwn
         call MeshSchema_csrrecv(NodebyFracOwn, 0, 400, VALSIZE_ZERO)
         ! CellbyFracOwn
         call MeshSchema_csrrecv(CellbyFracOwn, 0, 500, VALSIZE_ZERO)
         ! FracbyFracOwn
         call MeshSchema_csrrecv(FracbyFracOwn, 0, 600, VALSIZE_ZERO)

         ! WellInjbyNodeOwn
         call MeshSchema_csrrecv(WellInjbyNodeOwn, 0, 610, VALSIZE_ZERO)
         ! WellProdbyNodeOwn
         call MeshSchema_csrrecv(WellProdbyNodeOwn, 0, 620, VALSIZE_ZERO)
         ! MSWellbyNodeOwn
         call MeshSchema_csrrecv(MSWellbyNodeOwn, 0, 630, VALSIZE_ZERO)

         ! FacebyCellLocal
         call MeshSchema_csrrecv(FacebyCellLocal, 0, 700, VALSIZE_ZERO)
         ! FracbyCellLocal
         call MeshSchema_csrrecv(FracbyCellLocal, 0, 800, VALSIZE_ZERO)
         ! NodebyCellLocal
         call MeshSchema_csrrecv(NodebyCellLocal, 0, 900, VALSIZE_ZERO)
         ! NodebyFaceLocal
         call MeshSchema_csrrecv(NodebyFaceLocal, 0, 1000, VALSIZE_ZERO)

         ! NodebyWellInjLocal
         call MeshSchema_csrrecv(NodebyWellInjLocal, 0, 1100, VALSIZE_ZERO)
         ! NodebyWellProdLocal
         call MeshSchema_csrrecv(NodebyWellProdLocal, 0, 1200, VALSIZE_ZERO)
         ! NodebyMSWellLocal
         call MeshSchema_csrrecv(NodebyMSWellLocal, 0, 1300, VALSIZE_ZERO)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            call CommonType_deallocCSR(NodebyNodeOwn_Ncpus(i))
            call CommonType_deallocCSR(FracbyNodeOwn_Ncpus(i))
            call CommonType_deallocCSR(CellbyNodeOwn_Ncpus(i))

            call CommonType_deallocCSR(NodebyFracOwn_Ncpus(i))
            call CommonType_deallocCSR(CellbyFracOwn_Ncpus(i))
            call CommonType_deallocCSR(FracbyFracOwn_Ncpus(i))

            call CommonType_deallocCSR(WellInjbyNodeOwn_Ncpus(i))
            call CommonType_deallocCSR(WellProdbyNodeOwn_Ncpus(i))
            call CommonType_deallocCSR(MSWellbyNodeOwn_Ncpus(i))

            call CommonType_deallocCSR(FacebyCellLocal_Ncpus(i))
            call CommonType_deallocCSR(FracbyCellLocal_Ncpus(i))
            call CommonType_deallocCSR(NodebyCellLocal_Ncpus(i))
            call CommonType_deallocCSR(NodebyFaceLocal_Ncpus(i))

            call CommonType_deallocCSR(NodebyWellInjLocal_Ncpus(i))
            call CommonType_deallocCSR(NodebyWellProdLocal_Ncpus(i))
            call CommonType_deallocCSR(NodebyMSWellLocal_Ncpus(i))
         end do

         deallocate (NodebyNodeOwn_Ncpus)
         deallocate (FracbyNodeOwn_Ncpus)
         deallocate (CellbyNodeOwn_Ncpus)

         deallocate (NodebyFracOwn_Ncpus)
         deallocate (CellbyFracOwn_Ncpus)
         deallocate (FracbyFracOwn_Ncpus)

         deallocate (WellInjbyNodeOwn_Ncpus)
         deallocate (WellProdbyNodeOwn_Ncpus)
         deallocate (MSWellbyNodeOwn_Ncpus)

         deallocate (FacebyCellLocal_Ncpus)
         deallocate (FracbyCellLocal_Ncpus)
         deallocate (NodebyCellLocal_Ncpus)
         deallocate (NodebyFaceLocal_Ncpus)

         deallocate (NodebyWellInjLocal_Ncpus)
         deallocate (NodebyWellProdLocal_Ncpus)
         deallocate (NodebyMSWellLocal_Ncpus)
      end if

      ! Well data type

      call DefWell_mpi_register_well_data_description(MPI_WELLDATA_ID)

      ! Injectors
      if (commRank == 0) then
         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbWellInjLocal_Ncpus(i + 1)
            call MPI_Send(DataWellInjRes_Ncpus(1:Nb, i + 1), Nb, MPI_WELLDATA_ID, i, 210, ComPASS_COMM_WORLD, Ierr)
         end do
      end if
      Nb = NbWellInjLocal_Ncpus(commRank + 1)
      allocate (DataWellInjLocal(Nb))
      if (commRank == 0) then
         DataWellInjLocal(:) = DataWellInjRes_Ncpus(1:Nb, 1) ! proc=0, copy
      else
         call MPI_Recv(DataWellInjLocal, Nb, MPI_WELLDATA_ID, 0, 210, ComPASS_COMM_WORLD, stat, Ierr) ! proc >=1, receive
      end if

      ! Producers
      if (commRank == 0) then
         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbWellProdLocal_Ncpus(i + 1)
            call MPI_Send(DataWellProdRes_Ncpus(1:Nb, i + 1), Nb, MPI_WELLDATA_ID, i, 211, ComPASS_COMM_WORLD, Ierr)
         end do
      end if
      Nb = NbWellProdLocal_Ncpus(commRank + 1)
      allocate (DataWellProdLocal(Nb))
      if (commRank == 0) then
         DataWellProdLocal(:) = DataWellProdRes_Ncpus(1:Nb, 1) ! proc 0, copy
      else
         call MPI_Recv(DataWellProdLocal, Nb, MPI_WELLDATA_ID, 0, 211, ComPASS_COMM_WORLD, stat, Ierr) ! proc >=1, receive
      end if

      call MPI_Type_free(MPI_WELLDATA_ID, Ierr)

      ! MSWell data type

      call DefMSWell_mpi_register_mswell_data_description(MPI_MSWELLDATA_ID)

      ! MSWells
      if (commRank == 0) then
         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbMSWellLocal_Ncpus(i + 1)
            call MPI_Send(DataMSWellRes_Ncpus(1:Nb, i + 1), Nb, MPI_MSWELLDATA_ID, i, 212, ComPASS_COMM_WORLD, Ierr)
         end do
      end if
      Nb = NbMSWellLocal_Ncpus(commRank + 1)
      allocate (DataMSWellLocal(Nb))
      if (commRank == 0) then
         DataMSWellLocal(:) = DataMSWellRes_Ncpus(1:Nb, 1) ! proc 0, copy
      else
         call MPI_Recv(DataMSWellLocal, Nb, MPI_MSWELLDATA_ID, 0, 212, ComPASS_COMM_WORLD, stat, Ierr) ! proc >=1, receive
      end if

      call MPI_Type_free(MPI_MSWELLDATA_ID, Ierr)

      ! ************************************* !

      ! Send IdCellLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbCellLocal_Ncpus(i + 1)
            call MPI_Send(IdCellRes_Ncpus(i + 1)%Val, Nb, MPI_INTEGER, i, 11, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbCellLocal_Ncpus(1)
         allocate (IdCellLocal(Nb))
         IdCellLocal(:) = IdCellRes_Ncpus(1)%Val(:)

      else
         Nb = NbCellLocal_Ncpus(commRank + 1)
         allocate (IdCellLocal(Nb))
         call MPI_Recv(IdCellLocal, Nb, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (IdCellRes_Ncpus(i)%Val)
         end do
         deallocate (IdCellRes_Ncpus)
      end if

      ! Send IdFaceLocal
      if (commRank == 0) then

         ! send to proc >=1
         do i = 1, Ncpus - 1
            Nb = NbFaceLocal_Ncpus(i + 1)
            call MPI_Send(IdFaceRes_Ncpus(i + 1)%Val, Nb, MPI_INTEGER, i, 11, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbFaceLocal_Ncpus(1)
         allocate (IdFaceLocal(Nb))
         IdFaceLocal(:) = IdFaceRes_Ncpus(1)%Val(:)

      else
         Nb = NbFaceLocal_Ncpus(commRank + 1)
         allocate (IdFaceLocal(Nb))
         call MPI_Recv(IdFaceLocal, Nb, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (IdFaceRes_Ncpus(i)%Val)
         end do
         deallocate (IdFaceRes_Ncpus)
      end if

      ! ************************************* !

      ! new MPI type: MPI_IDNODE
      blen(1) = 4
      offsets(1) = 0
      oldtypes(1) = MPI_CHARACTER

      call MPI_Type_struct(1, blen, offsets, oldtypes, MPI_IDNODE, Ierr)
      call MPI_Type_commit(MPI_IDNODE, Ierr)

      ! Send IdNodeLocal
      if (commRank == 0) then

         ! send to proc >=1
         do i = 1, Ncpus - 1
            Nb = NbNodeLocal_Ncpus(i + 1)
            call MPI_Send(IdNodeRes_Ncpus(i + 1)%Val, Nb, MPI_IDNODE, i, 12, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbNodeLocal_Ncpus(1)
         allocate (IdNodeLocal(Nb))
         IdNodeLocal(:) = IdNodeRes_Ncpus(1)%Val(:)

      else   ! not master proc
         Nb = NbNodeLocal_Ncpus(commRank + 1)
         allocate (IdNodeLocal(Nb))
         call MPI_Recv(IdNodeLocal, Nb, MPI_IDNODE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      ! free MPI_IDNODE
      call MPI_Type_free(MPI_IDNODE, Ierr)

      ! free IdNode
      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (IdNodeRes_Ncpus(i)%Val)
         end do
         deallocate (IdNodeRes_Ncpus)
      end if

      ! ************************************ !

      ! Send FracToFace and FaceToFrac
      if (commRank == 0) then

         ! proc>=1, send
         do i = 1, Ncpus - 1
            Nb = NbFracLocal_Ncpus(i + 1)
            call MPI_Send(FracToFaceLocal_Ncpus(i + 1)%Val, Nb, MPI_INTEGER, i, 21, ComPASS_COMM_WORLD, Ierr) ! FracToFace

            Nb = NbFaceLocal_Ncpus(i + 1)
            call MPI_Send(FaceToFracLocal_Ncpus(i + 1)%Val, Nb, MPI_INTEGER, i, 22, ComPASS_COMM_WORLD, Ierr) ! FracToFace
         end do

         ! proc=0, copy
         Nb = NbFracLocal_Ncpus(1)
         allocate (FracToFaceLocal(Nb))
         FracToFaceLocal(:) = FracToFaceLocal_Ncpus(1)%Val(:) ! FracToFace

         Nb = NbFaceLocal_Ncpus(1)
         allocate (FaceToFracLocal(Nb))
         FaceToFracLocal(:) = FaceToFracLocal_Ncpus(1)%Val(:) ! FaceToFrac

      else
         Nb = NbFracLocal_Ncpus(commRank + 1)
         allocate (FracToFaceLocal(Nb))
         call MPI_Recv(FracToFaceLocal, Nb, MPI_INTEGER, 0, 21, ComPASS_COMM_WORLD, stat, Ierr) ! FracToFace

         Nb = NbFaceLocal_Ncpus(commRank + 1)
         allocate (FaceToFracLocal(Nb))
         call MPI_Recv(FaceToFracLocal, Nb, MPI_INTEGER, 0, 22, ComPASS_COMM_WORLD, stat, Ierr) ! FaceToFrac
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (FracToFaceLocal_Ncpus(i)%Val)
            deallocate (FaceToFracLocal_Ncpus(i)%Val)
         end do
         deallocate (FracToFaceLocal_Ncpus)
         deallocate (FaceToFracLocal_Ncpus)
      end if

      call MeshSchema_send_recv_FamilyDOFIdCOC

      ! send PorositeCell
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbCellLocal_Ncpus(i + 1)
            call MPI_Send(PorositeCell_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 11, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbCellLocal_Ncpus(1)
         allocate (PorositeCellLocal(Nb))
         PorositeCellLocal(:) = PorositeCell_Ncpus(1)%Val(:)

      else
         Nb = NbCellLocal_Ncpus(commRank + 1)
         allocate (PorositeCellLocal(Nb))
         call MPI_Recv(PorositeCellLocal, Nb, MPI_DOUBLE, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (PorositeCell_Ncpus(i)%Val)
         end do
         deallocate (PorositeCell_Ncpus)
      end if

      ! send PorositeFrac
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbFracLocal_Ncpus(i + 1)
            call MPI_Send(PorositeFrac_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 12, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         !print *, "DEBUG - Copying", Nb, "fracture porosity values:", PorositeFrac_Ncpus(1)%Val(:)
         Nb = NbFracLocal_Ncpus(1)
         allocate (PorositeFracLocal(Nb))
         PorositeFracLocal(:) = PorositeFrac_Ncpus(1)%Val(:)

      else
         Nb = NbFracLocal_Ncpus(commRank + 1)
         allocate (PorositeFracLocal(Nb))
         call MPI_Recv(PorositeFracLocal, Nb, MPI_DOUBLE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (PorositeFrac_Ncpus(i)%Val)
         end do
         deallocate (PorositeFrac_Ncpus)
      end if

      ! send PermCellLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbCellLocal_Ncpus(i + 1)
            call MPI_Send(PermCellLocal_Ncpus(i + 1)%Array3d, Nb*9, MPI_DOUBLE, i, 13, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbCellLocal_Ncpus(1)
         allocate (PermCellLocal(3, 3, Nb))
         PermCellLocal(:, :, :) = PermCellLocal_Ncpus(1)%Array3d(:, :, :)

      else
         Nb = NbCellLocal_Ncpus(commRank + 1)
         allocate (PermCellLocal(3, 3, Nb))
         call MPI_Recv(PermCellLocal, Nb*9, MPI_DOUBLE, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (PermCellLocal_Ncpus(i)%Array3d)
         end do
         deallocate (PermCellLocal_Ncpus)
      end if

      ! send PermFracLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbFracLocal_Ncpus(i + 1)
            call MPI_Send(PermFracLocal_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 14, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbFracLocal_Ncpus(1)
         allocate (PermFracLocal(Nb))
         PermFracLocal(:) = PermFracLocal_Ncpus(1)%Val(:)

      else
         Nb = NbFracLocal_Ncpus(commRank + 1)
         allocate (PermFracLocal(Nb))
         call MPI_Recv(PermFracLocal, Nb, MPI_DOUBLE, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (PermFracLocal_Ncpus(i)%Val)
         end do
         deallocate (PermFracLocal_Ncpus)
      end if

#ifdef _THERMIQUE_
      ! send CondThermalCellLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbCellLocal_Ncpus(i + 1)
            call MPI_Send(CondThermalCellLocal_Ncpus(i + 1)%Array3d, Nb*9, MPI_DOUBLE, i, 13, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbCellLocal_Ncpus(1)
         allocate (CondThermalCellLocal(3, 3, Nb))
         CondThermalCellLocal(:, :, :) = CondThermalCellLocal_Ncpus(1)%Array3d(:, :, :)

      else
         Nb = NbCellLocal_Ncpus(commRank + 1)
         allocate (CondThermalCellLocal(3, 3, Nb))
         call MPI_Recv(CondThermalCellLocal, Nb*9, MPI_DOUBLE, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (CondThermalCellLocal_Ncpus(i)%Array3d)
         end do
         deallocate (CondThermalCellLocal_Ncpus)
      end if

      ! send CondThermalFracLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbFracLocal_Ncpus(i + 1)
            call MPI_Send(CondThermalFracLocal_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 14, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbFracLocal_Ncpus(1)
         allocate (CondThermalFracLocal(Nb))
         CondThermalFracLocal(:) = CondThermalFracLocal_Ncpus(1)%Val(:)

      else
         Nb = NbFracLocal_Ncpus(commRank + 1)
         allocate (CondThermalFracLocal(Nb))
         call MPI_Recv(CondThermalFracLocal, Nb, MPI_DOUBLE, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (CondThermalFracLocal_Ncpus(i)%Val)
         end do
         deallocate (CondThermalFracLocal_Ncpus)
      end if

      ! send CellThermalSourceLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbCellLocal_Ncpus(i + 1)
            call MPI_Send(CellThermalSource_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 15, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbCellLocal_Ncpus(1)
         allocate (CellThermalSourceLocal(Nb))
         CellThermalSourceLocal = CellThermalSource_Ncpus(1)%Val

      else
         Nb = NbCellLocal_Ncpus(commRank + 1)
         allocate (CellThermalSourceLocal(Nb))
         call MPI_Recv(CellThermalSourceLocal, Nb, MPI_DOUBLE, 0, 15, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (CellThermalSource_Ncpus(i)%Val)
         end do
         deallocate (CellThermalSource_Ncpus)
      end if

      ! send FracThermalSourceLocal
      if (commRank == 0) then

         ! proc >=1, send
         do i = 1, Ncpus - 1
            Nb = NbFracLocal_Ncpus(i + 1)
            call MPI_Send(FracThermalSource_Ncpus(i + 1)%Val, Nb, MPI_DOUBLE, i, 16, ComPASS_COMM_WORLD, Ierr)
         end do

         ! proc=0, copy
         Nb = NbFracLocal_Ncpus(1)
         allocate (FracThermalSourceLocal(Nb))
         FracThermalSourceLocal = FracThermalSource_Ncpus(1)%Val

      else
         Nb = NbFracLocal_Ncpus(commRank + 1)
         allocate (FracThermalSourceLocal(Nb))
         call MPI_Recv(FracThermalSourceLocal, Nb, MPI_DOUBLE, 0, 16, ComPASS_COMM_WORLD, stat, Ierr)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            deallocate (FracThermalSource_Ncpus(i)%Val)
         end do
         deallocate (FracThermalSource_Ncpus)
      end if
#endif

      ! MPI TYPE for DataNodewell: MPI_DATANODEWELL (also used for MSWELLS)
      blocklen = 1
      arraytype(1:3) = MPI_INTEGER
      arraytype(4:5) = MPI_DOUBLE
      disp(1) = 0 !  = 0
      disp(2) = 4 !  + integer
      disp(3) = 8 !  + integer
      disp(4) = 12 !  + integer
      disp(5) = 20 ! + double

      ! Create and commit
      call MPI_Type_Create_Struct(5, blocklen, disp, arraytype, MPI_DATANODEWELL, Ierr)
      call MPI_Type_commit(MPI_DATANODEWELL, Ierr)

      ! send NodeDatabyWell
      if (commRank == 0) then

         do i = 1, Ncpus - 1
            call MeshSchema_csrdatawellsend(NodeDatabyWellInjLocal_Ncpus(i + 1), i, 100)
            call MeshSchema_csrdatawellsend(NodeDatabyWellProdLocal_Ncpus(i + 1), i, 200)
            call MeshSchema_csrdatawellsend(NodeDatabyMSWellLocal_Ncpus(i + 1), i, 300)
         end do

         call DefWell_csrdatawellcopy(NodeDatabyWellInjLocal_Ncpus(1), NodeDatabyWellInjLocal)
         call DefWell_csrdatawellcopy(NodeDatabyWellProdLocal_Ncpus(1), NodeDatabyWellProdLocal)
         call DefWell_csrdatawellcopy(NodeDatabyMSWellLocal_Ncpus(1), NodeDatabyMSWellLocal)
      else

         call MeshSchema_csrdatawellrecv(NodeDatabyWellInjLocal, 0, 100)
         call MeshSchema_csrdatawellrecv(NodeDatabyWellProdLocal, 0, 200)
         call MeshSchema_csrdatawellrecv(NodeDatabyMSWellLocal, 0, 300)
      end if

      if (commRank == 0) then
         do i = 1, Ncpus
            call DefWell_deallocCSRDataWell(NodeDatabyWellInjLocal_Ncpus(i))
            call DefWell_deallocCSRDataWell(NodeDatabyWellProdLocal_Ncpus(i))
            call DefWell_deallocCSRDataWell(NodeDatabyMSWellLocal_Ncpus(i))
         end do
         deallocate (NodeDatabyWellInjLocal_Ncpus)
         deallocate (NodeDatabyWellProdLocal_Ncpus)
         deallocate (NodeDatabyMSWellLocal_Ncpus)
      end if

      ! Free TYPE MPI_DATANODEWELL
      call MPI_Type_free(MPI_DATANODEWELL, Ierr)

   end subroutine MeshSchema_sendrecv

   ! Send csr1 to "dest" with "tag"
   subroutine MeshSchema_csrsend(csr1, dest, tag, valsize)

      type(CSR) :: csr1
      integer, intent(in) :: tag, dest, valsize
      integer :: Nb, Nnz, Ierr

      call MPI_Send(csr1%Nb, 1, MPI_INTEGER, dest, tag + 1, ComPASS_COMM_WORLD, Ierr)
      Nb = csr1%Nb
      call MPI_Send(csr1%Pt, Nb + 1, MPI_INTEGER, dest, tag + 2, ComPASS_COMM_WORLD, Ierr)
      Nnz = csr1%Pt(Nb + 1)
      call MPI_Send(csr1%Num, Nnz, MPI_INTEGER, dest, tag + 3, ComPASS_COMM_WORLD, Ierr)
      if (valsize == VALSIZE_NB) then
         call MPI_Send(csr1%Val, Nb, MPI_INTEGER, dest, tag + 4, ComPASS_COMM_WORLD, Ierr)
      else if (valsize == VALSIZE_NNZ) then
         call MPI_Send(csr1%Val, Nnz, MPI_INTEGER, dest, tag + 4, ComPASS_COMM_WORLD, Ierr)
      end if

   end subroutine MeshSchema_csrsend

   ! Recv csr2 from "source" with "tag"
   subroutine MeshSchema_csrrecv(csr2, source, tag, valsize)

      type(CSR), intent(inout) :: csr2
      integer, intent(in) :: tag, source, valsize

      integer :: Nb, Nnz, Ierr
      integer stat(MPI_STATUS_SIZE)

      call MPI_Recv(csr2%Nb, 1, MPI_INTEGER, source, tag + 1, ComPASS_COMM_WORLD, stat, Ierr)

      Nb = csr2%Nb
      allocate (csr2%Pt(Nb + 1))
      call MPI_Recv(csr2%Pt, Nb + 1, MPI_INTEGER, source, tag + 2, ComPASS_COMM_WORLD, stat, Ierr)

      Nnz = csr2%Pt(Nb + 1)
      allocate (csr2%Num(Nnz))
      call MPI_Recv(csr2%Num, Nnz, MPI_INTEGER, source, tag + 3, ComPASS_COMM_WORLD, stat, Ierr)

      if (valsize == VALSIZE_NB) then
         allocate (csr2%Val(Nb))
         call MPI_Recv(csr2%Val, Nb, MPI_INTEGER, source, tag + 4, ComPASS_COMM_WORLD, stat, Ierr)
      else if (valsize == VALSIZE_NNZ) then
         allocate (csr2%Val(Nnz))
         call MPI_Recv(csr2%Val, Nnz, MPI_INTEGER, source, tag + 4, ComPASS_COMM_WORLD, stat, Ierr)
      end if

   end subroutine MeshSchema_csrrecv

   ! Send csr1 to "dest" with "tag"
   ! %val type is DataNodeWell
   subroutine MeshSchema_csrdatawellsend(csr1, dest, tag)

      type(TYPE_CSRDataNodeWell) :: csr1
      integer, intent(in) :: tag, dest
      integer :: Nb, Nnz, Ierr

      call MPI_Send(csr1%Nb, 1, MPI_INTEGER, dest, tag + 1, ComPASS_COMM_WORLD, Ierr)
      Nb = csr1%Nb
      call MPI_Send(csr1%Pt, Nb + 1, MPI_INTEGER, dest, tag + 2, ComPASS_COMM_WORLD, Ierr)
      Nnz = csr1%Pt(Nb + 1)
      call MPI_Send(csr1%Num, Nnz, MPI_INTEGER, dest, tag + 3, ComPASS_COMM_WORLD, Ierr)

      call MPI_Send(csr1%Val, Nnz, MPI_DATANODEWELL, dest, tag + 4, ComPASS_COMM_WORLD, Ierr)

   end subroutine MeshSchema_csrdatawellsend

   ! Recv csr2 from "source" with "tag"
   ! %val type is DataNodeWell
   subroutine MeshSchema_csrdatawellrecv(csr2, source, tag)

      type(TYPE_CSRDataNodeWell), intent(inout) :: csr2
      integer, intent(in) :: tag, source

      integer :: Nb, Nnz, Ierr
      integer stat(MPI_STATUS_SIZE)

      call MPI_Recv(csr2%Nb, 1, MPI_INTEGER, source, tag + 1, ComPASS_COMM_WORLD, stat, Ierr)

      Nb = csr2%Nb
      allocate (csr2%Pt(Nb + 1))
      call MPI_Recv(csr2%Pt, Nb + 1, MPI_INTEGER, source, tag + 2, ComPASS_COMM_WORLD, stat, Ierr)

      Nnz = csr2%Pt(Nb + 1)
      allocate (csr2%Num(Nnz))
      call MPI_Recv(csr2%Num, Nnz, MPI_INTEGER, source, tag + 3, ComPASS_COMM_WORLD, stat, Ierr)

      allocate (csr2%Val(Nnz))
      call MPI_Recv(csr2%Val, Nnz, MPI_DATANODEWELL, source, tag + 4, ComPASS_COMM_WORLD, stat, Ierr)

   end subroutine MeshSchema_csrdatawellrecv

   subroutine MeshSchema_send_FamilyDOFIdCOC(fidcoc, mpi_struct_id, dest, tag)

      type(FamilyDOFIdCOC), intent(in) :: fidcoc
      integer, intent(in) :: mpi_struct_id, dest, tag

      integer n, Ierr

      n = size(fidcoc%offsets)
      call MPI_Send(n, 1, MPI_INTEGER, dest, tag + 1, ComPASS_COMM_WORLD, Ierr)
      call MPI_Send(fidcoc%offsets, n, MPI_INTEGER, dest, tag + 2, ComPASS_COMM_WORLD, Ierr)
#ifndef NDEBUG
      if (.not. check_FamilyDOFIdCOC(fidcoc)) &
         call CommonMPI_abort("Inconsistent input COC")
#endif
      n = size(fidcoc%ids)
      call MPI_Send(fidcoc%ids, n, mpi_struct_id, dest, tag + 3, ComPASS_COMM_WORLD, Ierr)

   end subroutine MeshSchema_send_FamilyDOFIdCOC

   subroutine MeshSchema_recv_FamilyDOFIdCOC(fidcoc, mpi_struct_id, source, tag)

      type(FamilyDOFIdCOC), intent(inout) :: fidcoc
      integer, intent(in) :: mpi_struct_id, source, tag

      integer n, Ierr
      integer stat(MPI_STATUS_SIZE)

      call MPI_Recv(n, 1, MPI_INTEGER, source, tag + 1, ComPASS_COMM_WORLD, stat, Ierr)
      allocate (fidcoc%offsets(n))
      call MPI_Recv(fidcoc%offsets, n, MPI_INTEGER, source, tag + 2, ComPASS_COMM_WORLD, stat, Ierr)
      n = fidcoc%offsets(n)
      allocate (fidcoc%ids(n))
      call MPI_Recv(fidcoc%ids, n, mpi_struct_id, source, tag + 3, ComPASS_COMM_WORLD, stat, Ierr)

   end subroutine MeshSchema_recv_FamilyDOFIdCOC

   subroutine MeshSchema_send_recv_FamilyDOFIdCOC()

      integer :: i, Ierr
      integer :: mpi_FamilyDOFIdCOC_id ! FIXME: use MPI_DATATYPE

      mpi_FamilyDOFIdCOC_id = create_FamilyDOFId_MPI_struct()
      call MPI_Type_commit(mpi_FamilyDOFIdCOC_id, Ierr)

      if (commRank == 0) then
         do i = 1, Ncpus - 1
            call MeshSchema_send_FamilyDOFIdCOC(NumNodebyProc_Ncpus(i + 1), mpi_FamilyDOFIdCOC_id, i, 100)
            call MeshSchema_send_FamilyDOFIdCOC(NumFracbyProc_Ncpus(i + 1), mpi_FamilyDOFIdCOC_id, i, 300)
            call MeshSchema_send_FamilyDOFIdCOC(NumWellInjbyProc_Ncpus(i + 1), mpi_FamilyDOFIdCOC_id, i, 500)
            call MeshSchema_send_FamilyDOFIdCOC(NumWellProdbyProc_Ncpus(i + 1), mpi_FamilyDOFIdCOC_id, i, 700)
         end do
         NumNodebyProc = NumNodebyProc_Ncpus(1)
         NumFracbyProc = NumFracbyProc_Ncpus(1)
         NumWellInjbyProc = NumWellInjbyProc_Ncpus(1)
         NumWellProdbyProc = NumWellProdbyProc_Ncpus(1)
      else
         call MeshSchema_recv_FamilyDOFIdCOC(NumNodebyProc, mpi_FamilyDOFIdCOC_id, 0, 100)
         call MeshSchema_recv_FamilyDOFIdCOC(NumFracbyProc, mpi_FamilyDOFIdCOC_id, 0, 300)
         call MeshSchema_recv_FamilyDOFIdCOC(NumWellInjbyProc, mpi_FamilyDOFIdCOC_id, 0, 500)
         call MeshSchema_recv_FamilyDOFIdCOC(NumWellProdbyProc, mpi_FamilyDOFIdCOC_id, 0, 700)
      end if

      ! CHECKME: from Fortran 2003 nested allocatable are automatically deallocated
      if (commRank == 0) then
         do i = 1, Ncpus
            call free_ComPASS_struct(NumNodebyProc_Ncpus(i))
            call free_ComPASS_struct(NumFracbyProc_Ncpus(i))
            call free_ComPASS_struct(NumWellInjbyProc_Ncpus(i))
            call free_ComPASS_struct(NumWellProdbyProc_Ncpus(i))
         end do
         deallocate (NumNodebyProc_Ncpus)
         deallocate (NumFracbyProc_Ncpus)
         deallocate (NumWellInjbyProc_Ncpus)
         deallocate (NumWellProdbyProc_Ncpus)
      end if

      call MPI_Type_free(mpi_FamilyDOFIdCOC_id, Ierr)

   end subroutine MeshSchema_send_recv_FamilyDOFIdCOC

   ! Output:
   !  NumNodebyEdgebyWellLocal
   ! Input:
   !  NodeDatabyWellLocal
   ! Local list of Edge by Well (own+ghost)
   ! List is oriented (Id_parent, Id_son)
#ifdef NDEBUG
   pure &
#endif
      subroutine MeshSchema_NumNodebyEdgebyWellLocal(NodeDatabyWellLocal, NbEdgebyWellLocal, NumNodebyEdgebyWellLocal)

      type(TYPE_CSRDataNodeWell), intent(in) :: NodeDatabyWellLocal
      integer, allocatable, dimension(:), intent(inout)     :: NbEdgebyWellLocal
      integer, allocatable, dimension(:, :, :), intent(inout) :: NumNodebyEdgebyWellLocal

      integer :: NbWellLocal, i, j, NbEdgemax, comptNode

      NbWellLocal = NodeDatabyWellLocal%Nb ! Number of well inj
      allocate (NbEdgebyWellLocal(NbWellLocal))
      NbEdgemax = 0

      do i = 1, NbWellLocal
         ! Number of edges = (Number of nodes in local well i) - 1
         NbEdgebyWellLocal(i) = NodeDatabyWellLocal%Pt(i + 1) - NodeDatabyWellLocal%Pt(i) - 1
         NbEdgemax = max(NbEdgemax, NbEdgebyWellLocal(i))
      enddo

      allocate (NumNodebyEdgebyWellLocal(2, NbEdgemax, NbWellLocal))
      NumNodebyEdgebyWellLocal(:, :, :) = -1

      do i = 1, NbWellLocal

         ! comptNode = 0
         ! ! loop over every nodes of local well, minus the head node of well
         ! do j=NodeDatabyWellLocal%Pt(i)+1,NodeDatabyWellLocal%Pt(i+1)-1
         !    comptNode = comptNode + 1
         !    NumNodebyEdgebyWellLocal(1,comptNode,i) = NodeDatabyWellLocal%Val(j)%Parent
         !    NumNodebyEdgebyWellLocal(2,comptNode,i) = NodeDatabyWellLocal%Num(j)
         ! enddo

         comptNode = 0

         ! loop over every nodes of local well, minus the head node of well
         do j = 1, NodeDatabyWellLocal%Pt(i + 1) - NodeDatabyWellLocal%Pt(i) - 1
            NumNodebyEdgebyWellLocal(1, comptNode + j, i) = NodeDatabyWellLocal%Val(j + NodeDatabyWellLocal%Pt(i))%Parent
            NumNodebyEdgebyWellLocal(2, comptNode + j, i) = NodeDatabyWellLocal%Num(j + NodeDatabyWellLocal%Pt(i))
         enddo
         comptNode = comptNode + NodeDatabyWellLocal%Pt(i + 1) - NodeDatabyWellLocal%Pt(i) - 1
      enddo

   end subroutine MeshSchema_NumNodebyEdgebyWellLocal

   ! center of Cell
   subroutine MeshSchema_XCellLocal

      integer :: k, m, nbCellLocal
      double precision, dimension(3) :: xk(3)

      nbCellLocal = NbCellLocal_Ncpus(commRank + 1)

      allocate (XCellLocal(3, nbCellLocal))

      ! boucle sur les mailles
      do k = 1, nbCellLocal

         ! center of cell
         xk(:) = 0.d0
         do m = NodebyCellLocal%Pt(k) + 1, NodebyCellLocal%Pt(k + 1)
            xk(:) = xk(:) + XNodeLocal(:, NodebyCellLocal%Num(m))
         enddo
         xk(:) = xk(:)/dble(NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k))

         XCellLocal(:, k) = xk(:)
      enddo

   end subroutine MeshSchema_XCellLocal

   ! Vol of Cell
   subroutine MeshSchema_VolCellLocal

      ! calcul du volume et du centre de gravite
      ! maille polyèdrique quelconque non necessairement convexe
      ! faces non planes (decoupée en triangles avec un point au centre)

      integer :: i, j, k, m, n1, n2
      double precision :: volk, volT
      double precision, dimension(3) :: yk, xT, x1, x2, xs, e0, e1, e2, e3

      integer :: Ierr, errcode ! used for MPI_Abort

      ! check if XCellLocal is computed
      if (allocated(XCellLocal) .eqv. .false.) then
         if (commRank == 0) then
            print *, "XCellLocal: center of cell not computed"
         end if

         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      allocate (VolCellLocal(NbCellLocal_Ncpus(commRank + 1)))
      VolCellLocal(:) = 0.d0

      ! boucle sur les mailles
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         volk = 0.d0

         ! center of cell
         yk(:) = XCellLocal(:, k)

         ! boucle sur les faces i de la maille k
         do j = FacebyCellLocal%Pt(k) + 1, FacebyCellLocal%Pt(k + 1)
            i = FacebyCellLocal%Num(j)

            ! isobarycentre de la face
            xs(:) = XFaceLocal(:, i)

            ! boucle sur les nodes n1 de la face i
            do m = NodebyFaceLocal%Pt(i) + 1, NodebyFaceLocal%Pt(i + 1)
               n1 = NodebyFaceLocal%Num(m)
               x1(:) = XNodeLocal(:, n1)

               if (m == NodebyFaceLocal%Pt(i + 1)) then
                  n2 = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(i) + 1)
               else
                  n2 = NodebyFaceLocal%Num(m + 1)
               endif

               x2(:) = XNodeLocal(:, n2)
               xT(:) = (x1(:) + x2(:) + xs(:) + yk(:))/4.d0

               e0(:) = xT(:) - yk(:)
               e1(:) = x1(:) - yk(:)
               e2(:) = x2(:) - yk(:)
               e3(:) = xs(:) - yk(:)

               volT = e1(1)*e2(2)*e3(3) + e2(1)*e3(2)*e1(3) + e3(1)*e1(2)*e2(3) &
                      - e1(1)*e2(3)*e3(2) - e2(1)*e3(3)*e1(2) - e3(1)*e1(3)*e2(2)

               volT = abs(volT)/6.d0
               volk = volk + volT
            enddo
         enddo

#ifndef NDEBUG
         if (volk < 1E-10) then
            print *, "DEBUG - Small cell volume for cell", k
            do m = NodebyCellLocal%Pt(k) + 1, NodebyCellLocal%Pt(k + 1)
               n1 = NodebyCellLocal%Num(m)
               print *, "        Node", n1, ":", XNodeLocal(:, n1)
            end do
         endif
#endif

         VolCellLocal(k) = volk

      enddo

   end subroutine MeshSchema_VolCellLocal

   ! center of frac
   subroutine MeshSchema_XFaceLocal

      integer :: i, m
      integer :: nbFaceLocal
      double precision :: xf(3)

      nbFaceLocal = NbFaceLocal_Ncpus(commRank + 1)

      allocate (XFaceLocal(3, nbFaceLocal))    ! center of face

      ! boucle sur les face frac
      do i = 1, nbFaceLocal

         ! isobarycentre de la face
         xf(:) = 0.d0
         do m = NodebyFaceLocal%Pt(i) + 1, NodebyFaceLocal%Pt(i + 1)
            xf(:) = xf(:) + XNodeLocal(:, NodebyFaceLocal%Num(m))
         enddo
         XFaceLocal(:, i) = xf(:)/dble(NodebyFaceLocal%Pt(i + 1) - NodebyFaceLocal%Pt(i))

      end do ! end of loop frac

   end subroutine MeshSchema_XFaceLocal

   function MeshSchema_local_face_surface_from_nodes(barycenter, nodes) result(surface)

      double precision, dimension(3), intent(in) :: barycenter
      integer, intent(in) :: nodes(:)
      integer :: edges(size(nodes) + 1)
      real(c_double) :: surface

      integer :: i, nbnodes
      double precision, dimension(3) :: x1, x2 !, xt ! coordinates
      !double precision :: contribution12f

      surface = 0.d0
      ! loop on face edges
      nbnodes = size(nodes)
      edges(1:nbnodes) = nodes
      edges(nbnodes + 1) = nodes(1)
      do i = 1, nbnodes
         x1(:) = XNodeLocal(:, edges(i))
         x2(:) = XNodeLocal(:, edges(i + 1))
         surface = surface + MeshSchema_triangle_area(x1, x2, barycenter)
      end do

   end function MeshSchema_local_face_surface_from_nodes

   function MeshSchema_local_face_surface(fk) result(surface) &
      bind(C, name="face_surface")

      integer(c_int), intent(in) :: fk
      real(c_double) :: surface

      surface = MeshSchema_local_face_surface_from_nodes( &
                XFaceLocal(:, fk), &
                NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + 1:NodebyFaceLocal%Pt(fk + 1)) &
                )

   end function MeshSchema_local_face_surface

   subroutine MeshSchema_SurfFracLocal

      integer :: ifrac, errcode, Ierr

      ! check if XFaceLocal is computed
      if (allocated(XFaceLocal) .eqv. .false.) then
         if (commRank == 0) then
            print *, "XFaceLocal: center of cell not computed"
         end if

         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      allocate (SurfFracLocal(NbFracLocal_Ncpus(commRank + 1))) ! surf of face

      SurfFracLocal(:) = 0.d0

      do ifrac = 1, NbFracLocal_Ncpus(commRank + 1)
         SurfFracLocal(ifrac) = MeshSchema_local_face_surface(FracToFaceLocal(ifrac))
      end do

   end subroutine MeshSchema_SurfFracLocal

   subroutine MeshSchema_collect_fracture_nodes

      integer :: k, frac, face, pnode, nbfractures, nbnodes
      integer :: Ierr, errcode ! FIXME: used for MPI_Abort but not assigned

      nbfractures = NbFracLocal_Ncpus(commRank + 1)

      !print *, "DEBUG - Collecting nodes of", nbfractures, "fractures"

      if (allocated(NodebyFractureLocal%Pt)) then
         deallocate (NodebyFractureLocal%Pt)
      end if
      NodebyFractureLocal%Nb = nbfractures
      allocate (NodebyFractureLocal%Pt(nbfractures + 1))

      !print *, "DEBUG - Allocated CSR pointer"

      NodebyFractureLocal%Pt(1) = 0
      do frac = 1, nbfractures
         face = FracToFaceLocal(frac)
         nbnodes = NodebyFaceLocal%Pt(face + 1) - NodebyFaceLocal%Pt(face)
         !print *, "DEBUG - Fracture", frac, "has", nbnodes, "nodes"
         NodebyFractureLocal%Pt(frac + 1) = NodebyFractureLocal%Pt(frac) + nbnodes
      enddo

      if (allocated(NodebyFractureLocal%Num)) then
         deallocate (NodebyFractureLocal%Num)
      end if
      allocate (NodebyFractureLocal%Num(NodebyFractureLocal%Pt(nbfractures + 1)))
      k = 0
      do frac = 1, nbfractures
         face = FracToFaceLocal(frac)
         do pnode = NodebyFaceLocal%Pt(face) + 1, NodebyFaceLocal%Pt(face + 1)
            k = k + 1
            NodebyFractureLocal%Num(k) = NodebyFaceLocal%Num(pnode)
         end do
      end do

      if (k /= NodebyFractureLocal%Pt(nbfractures + 1)) then
         if (commRank == 0) then
            print *, "MeshSchema_collect_fracture_nodes: something went wrong"
         end if
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      !print *, "DEBUG - Collecting fracture nodes. DONE!"

   end subroutine MeshSchema_collect_fracture_nodes

   function MeshSchema_triangle_area(A, B, C) result(area)

      double precision, dimension(3), intent(in) :: A, B, C
      double precision, dimension(3) :: AB, AC
      double precision :: area

      AB = B - A
      AC = C - A
      area = dsqrt( &
             (AB(2)*AC(3) - AB(3)*AC(2))**2 + &
             (AB(3)*AC(1) - AB(1)*AC(3))**2 + &
             (AB(1)*AC(2) - AB(2)*AC(1))**2)/2.d0

   end function MeshSchema_triangle_area

   ! max number of nodes in a cell
   subroutine MeshSchema_NbNodeCellMax

      integer :: k, Nb

      NbNodeCellMax = 0
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         Nb = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         if (NbNodeCellMax < Nb) then
            NbNodeCellMax = Nb
         end if
      end do

   end subroutine MeshSchema_NbNodeCellMax

   ! max number of frac in a cell
   subroutine MeshSchema_NbFracCellMax

      integer :: k, Nb

      NbFracCellMax = 0
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         Nb = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)
         if (NbFracCellMax < Nb) then
            NbFracCellMax = Nb
         end if
      end do

   end subroutine MeshSchema_NbFracCellMax

   ! max number of nodes in a cell
   subroutine MeshSchema_NbNodeFaceMax

      integer :: i, Nb

      NbNodeFaceMax = 0
      do i = 1, NbFaceLocal_Ncpus(commRank + 1)

         Nb = NodebyFaceLocal%Pt(i + 1) - NodebyFaceLocal%Pt(i)
         if (NbNodeFaceMax < Nb) then
            NbNodeFaceMax = Nb
         end if
      end do

   end subroutine MeshSchema_NbNodeFaceMax

   ! Free mesh and some connecivities
   subroutine MeshSchema_Free

      deallocate (NbCellLocal_Ncpus)
      deallocate (NbCellOwn_Ncpus)
      deallocate (NbFaceLocal_Ncpus)
      deallocate (NbFaceOwn_Ncpus)
      deallocate (NbNodeLocal_Ncpus)
      deallocate (NbNodeOwn_Ncpus)
      deallocate (NbFracLocal_Ncpus)
      deallocate (NbFracOwn_Ncpus)

      call CommonType_deallocCSR(FacebyCellLocal)
      call CommonType_deallocCSR(FracbyCellLocal)
      call CommonType_deallocCSR(NodebyCellLocal)
      call CommonType_deallocCSR(NodebyFaceLocal)

      deallocate (XNodeLocal)
      deallocate (NodeFlagsLocal)
      deallocate (CellFlagsLocal)
      deallocate (FaceFlagsLocal)
      deallocate (CellTypesLocal)
      deallocate (FaceTypesLocal)

      nullify (NodeDarcyRocktypesLocal)
      nullify (FracDarcyRocktypesLocal)
      nullify (CellDarcyRocktypesLocal)
      deallocate (AllDarcyRocktypesLocal)

#ifdef _THERMIQUE_
      nullify (NodeFourierRocktypesLocal)
      nullify (FracFourierRocktypesLocal)
      nullify (CellFourierRocktypesLocal)
      deallocate (AllFourierRocktypesLocal)
#endif

      deallocate (IdCellLocal)
      deallocate (IdFaceLocal)
      deallocate (IdNodeLocal)
#ifdef _WITH_FREEFLOW_STRUCTURES_
      deallocate (IsFreeflowNode)
      deallocate (AtmState)
#endif

      deallocate (FracToFaceLocal)
      deallocate (FaceToFracLocal)

      call free_ComPASS_struct(NumNodebyProc)
      call free_ComPASS_struct(NumFracbyProc)
      call free_ComPASS_struct(NumWellInjbyProc)
      call free_ComPASS_struct(NumWellProdbyProc)

      deallocate (NbEdgebyWellInjLocal)
      deallocate (NbEdgebyWellProdLocal)
      deallocate (NumNodebyEdgebyWellInjLocal)
      deallocate (NumNodebyEdgebyWellProdLocal)

      deallocate (XCellLocal)
      deallocate (XFaceLocal)

      deallocate (VolCellLocal)
#ifdef _WITH_FREEFLOW_STRUCTURES_
      deallocate (SurfFreeFlowLocal)
#endif
      deallocate (SurfFracLocal)

      ! the follwoing free could be put after VAGFrac?
      ! but we keep them just if user wants to use them
      deallocate (PorositeCellLocal)
      deallocate (PorositeFracLocal)
      deallocate (PermCellLocal)
      deallocate (PermFracLocal)
#ifdef _THERMIQUE_
      deallocate (CondThermalCellLocal)
      deallocate (CondThermalFracLocal)
#endif

#ifdef _THERMIQUE_
      deallocate (CellThermalSourceLocal)
      deallocate (FracThermalSourceLocal)
#endif

      ! the folllowing free could be put after Assembly
      call CommonType_deallocCSR(NodebyNodeOwn)
      call CommonType_deallocCSR(FracbyNodeOwn)
      call CommonType_deallocCSR(CellbyNodeOwn)

      call CommonType_deallocCSR(NodebyFracOwn)
      call CommonType_deallocCSR(CellbyFracOwn)
      call CommonType_deallocCSR(FracbyFracOwn)

   end subroutine MeshSchema_Free

end module MeshSchema
