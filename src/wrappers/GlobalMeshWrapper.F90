!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

    module GlobalMeshWrapper

       use, intrinsic :: iso_c_binding
       use mpi, only: MPI_abort

       use CommonTypesWrapper, only: cpp_COC, cpp_array_wrapper, &
                                     retrieve_coc
       use CommonType, only: CSR

       use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

       use DefModel, only: NbComp

       use GlobalMesh, only: &
          NbNode, NbCell, &
          IdNode, IdFace, &
          NodebyCell, NodebyFace, CellbyNode, FracbyNode, &
          CellbyCell, CellbyFace, FacebyCell, FacebyNode, &
          XNode, CellTypes, FaceTypes, NodeFlags, CellFlags, FaceFlags, &
          NbDirNodeP, NbDirNodeT, &
          PermCell, PermFrac, &
          NodeRocktype, CellRocktype, FracRocktype, &
          CondThermalCell, CondThermalFrac, &
          PorositeCell, PorositeFrac, &
          CellThermalSource, CellComponentSource

       use MeshSchema, only: &
          NodebyFractureLocal, FacebyCellLocal, NodebyCellLocal, NodebyFaceLocal

       use GlobalMesh, only: &
          GlobalMesh_create_mesh, GlobalMesh_allocate_rocktype

       use DefWell, only: DefWell_Make_ComputeWellIndex
       use DefMSWell, only: DefMSWell_Make_ComputeWellIndex

       implicit none

       integer :: Ierr, errcode

       type, bind(C) :: cpp_MeshConnectivity
          type(cpp_COC) :: NodebyCell
          type(cpp_COC) :: NodebyFace
          type(cpp_COC) :: FacebyNode
          type(cpp_COC) :: FacebyCell
          type(cpp_COC) :: CellbyNode
          type(cpp_COC) :: CellbyFace
          type(cpp_COC) :: CellbyCell
       end type cpp_MeshConnectivity

       public :: &
          GlobalMesh_count_dirichlet_nodes, &
          get_global_number_of_nodes, &
          get_global_number_of_cells, &
          check_mesh_allocation, &
          retrieve_global_vertices, &
          retrieve_global_nodeflags, &
          retrieve_global_cellflags, &
          retrieve_global_faceflags, &
          GlobalMesh_allocate_rocktype_from_C, &
          retrieve_global_cell_rocktypes, &
          retrieve_global_node_rocktypes, &
          retrieve_global_fracture_rocktypes, &
          retrieve_global_celltypes, &
          retrieve_global_facetypes, &
          retrieve_global_id_faces, &
          retrieve_global_mesh_connectivity, &
          retrieve_mesh_connectivity, &
          retrieve_nodes_by_fractures, &
          retrieve_global_cell_porosity, &
          retrieve_global_fracture_porosity, &
          retrieve_global_cell_permeability, &
          retrieve_global_fracture_permeability, &
#ifdef _THERMIQUE_
          retrieve_global_cell_thermal_conductivity, &
          retrieve_global_fracture_thermal_conductivity, &
#endif
          retrieve_global_id_node, &
          DefWell_make_compute_well_index_from_C, &
          GlobalMesh_create_mesh_from_C

    contains

       subroutine GlobalMesh_count_dirichlet_nodes() &
          bind(C, name="GlobalMesh_count_dirichlet_nodes")

          integer :: i

          NbDirNodeP = 0
#ifdef _THERMIQUE_
          NbDirNodeT = 0
#endif
          do i = 1, NbNode
             if (IdNode(i)%P .eq. "d") NbDirNodeP = NbDirNodeP + 1
#ifdef _THERMIQUE_
             if (IdNode(i)%T .eq. "d") NbDirNodeT = NbDirNodeT + 1
#endif
          end do

       end subroutine GlobalMesh_count_dirichlet_nodes

       subroutine retrieve_global_vertices(cpp_array) &
          bind(C, name="retrieve_global_vertices")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Mesh is supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(XNode)) then
             print *, "Mesh vertices are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(XNode(1, 1))
          cpp_array%n = size(XNode, 2)

       end subroutine retrieve_global_vertices

       subroutine retrieve_global_nodeflags(cpp_array) &
          bind(C, name="retrieve_global_nodeflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(NodeFlags)) then
             print *, "Node flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(NodeFlags(1))
          cpp_array%n = size(NodeFlags)

       end subroutine retrieve_global_nodeflags

       subroutine retrieve_global_cellflags(cpp_array) &
          bind(C, name="retrieve_global_cellflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(CellFlags)) then
             print *, "cell flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellFlags(1))
          cpp_array%n = size(CellFlags)

       end subroutine retrieve_global_cellflags

       subroutine retrieve_global_faceflags(cpp_array) &
          bind(C, name="retrieve_global_faceflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(FaceFlags)) then
             print *, "Face flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FaceFlags(1))
          cpp_array%n = size(FaceFlags)

       end subroutine retrieve_global_faceflags

       subroutine retrieve_global_celltypes(cpp_array) &
          bind(C, name="retrieve_global_celltypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(CellTypes)) then
             print *, "Cell types are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellTypes(1))
          cpp_array%n = size(CellTypes)

       end subroutine retrieve_global_celltypes

       subroutine retrieve_global_facetypes(cpp_array) &
          bind(C, name="retrieve_global_facetypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(FaceTypes)) then
             print *, "Face types are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FaceTypes(1))
          cpp_array%n = size(FaceTypes)

       end subroutine retrieve_global_facetypes

       subroutine GlobalMesh_allocate_rocktype_from_C() &
          bind(C, name="GlobalMesh_allocate_rocktype")

          call GlobalMesh_allocate_rocktype

       end subroutine GlobalMesh_allocate_rocktype_from_C

       subroutine retrieve_global_cell_rocktypes(cpp_array) &
          bind(C, name="retrieve_global_cell_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(CellRocktype)) then
             print *, "cell rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellRocktype)
          cpp_array%n = size(CellRocktype, 2)

       end subroutine retrieve_global_cell_rocktypes

       subroutine retrieve_global_node_rocktypes(cpp_array) &
          bind(C, name="retrieve_global_node_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(NodeRocktype)) then
             print *, "node rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(NodeRocktype)
          cpp_array%n = size(NodeRocktype, 2)

       end subroutine retrieve_global_node_rocktypes

       subroutine retrieve_global_fracture_rocktypes(cpp_array) &
          bind(C, name="retrieve_global_fracture_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             !CHECKME: Maybe MPI_abort would be better here
             !buffer%p = c_null_ptr
             !buffer%n = 0
             print *, "Global values are supposed to be read by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(FracRocktype)) then
             print *, "fracture rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FracRocktype)
          cpp_array%n = size(FracRocktype, 2)

       end subroutine retrieve_global_fracture_rocktypes

       subroutine retrieve_global_cell_molar_sources(cpp_array) &
          bind(C, name="retrieve_global_cell_molar_sources")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) &
             call CommonMPI_abort("global values must be read by master process")
          if (.not. allocated(CellComponentSource)) &
             call CommonMPI_abort("cell component sources are not allocated")
          if (size(CellComponentSource) /= NbCell*NbComp) & !FIXME tester chaque dimension
             call CommonMPI_abort("cell component sources have inconsistent size")

          cpp_array%p = c_loc(CellComponentSource(1, 1))
          cpp_array%n = size(CellComponentSource, 2)

       end subroutine retrieve_global_cell_molar_sources

       subroutine retrieve_cell_heat_source(cpp_array) &
          bind(C, name="retrieve_cell_heat_source")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) &
             call CommonMPI_abort("global values must be read by master process")
          if (.not. allocated(CellThermalSource)) &
             call CommonMPI_abort("cell heat sources are not allocated")
          if (size(CellThermalSource) /= NbCell) &
             call CommonMPI_abort("cell heat sources have inconsistent size")

          cpp_array%p = c_loc(CellThermalSource(1))
          cpp_array%n = size(CellThermalSource)

       end subroutine retrieve_cell_heat_source

       subroutine retrieve_global_id_faces(cpp_array) &
          bind(C, name="retrieve_global_id_faces")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(IdFace)) then
             print *, "Faces ids are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(IdFace(1))
          cpp_array%n = size(IdFace)

       end subroutine retrieve_global_id_faces

       subroutine retrieve_global_cell_porosity(cpp_array) &
          bind(C, name="retrieve_global_cell_porosity")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(PorositeCell)) then
             print *, "Cell porosity array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PorositeCell(1))
          cpp_array%n = size(PorositeCell)

       end subroutine retrieve_global_cell_porosity

       subroutine retrieve_global_fracture_porosity(cpp_array) &
          bind(C, name="retrieve_global_fracture_porosity")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(PorositeFrac)) then
             print *, "face porosity array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PorositeFrac(1))
          cpp_array%n = size(PorositeFrac)

       end subroutine retrieve_global_fracture_porosity

       subroutine get_global_number_of_nodes(n) &
          bind(C, name="get_global_number_of_nodes")

          integer(c_long), intent(out) :: n
          n = NbNode

       end subroutine get_global_number_of_nodes

       subroutine get_global_number_of_cells(n) &
          bind(C, name="get_global_number_of_cells")

          integer(c_long), intent(out) :: n
          n = NbCell

       end subroutine get_global_number_of_cells

       subroutine retrieve_global_cell_permeability(cpp_array) &
          bind(C, name="retrieve_global_cell_permeability")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(PermCell)) then
             print *, "Cell permeability array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PermCell(1, 1, 1))
          cpp_array%n = size(PermCell, 3)

       end subroutine retrieve_global_cell_permeability

       subroutine retrieve_global_fracture_permeability(cpp_array) &
          bind(C, name="retrieve_global_fracture_permeability")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(PermFrac)) then
             print *, "face permeability array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PermFrac(1))
          cpp_array%n = size(PermFrac)

       end subroutine retrieve_global_fracture_permeability

#ifdef _THERMIQUE_

       subroutine retrieve_global_cell_thermal_conductivity(cpp_array) &
          bind(C, name="retrieve_global_cell_thermal_conductivity")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(CondThermalCell)) then
             print *, "Cell thermal conductivity array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CondThermalCell(1, 1, 1))
          cpp_array%n = size(CondThermalCell, 3)

       end subroutine retrieve_global_cell_thermal_conductivity

       subroutine retrieve_global_fracture_thermal_conductivity(cpp_array) &
          bind(C, name="retrieve_global_fracture_thermal_conductivity")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Global mesh is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(CondThermalFrac)) then
             print *, "face thermal conductivity array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CondThermalFrac(1))
          cpp_array%n = size(CondThermalFrac)

       end subroutine retrieve_global_fracture_thermal_conductivity

#endif

       subroutine retrieve_global_id_node(cpp_array) &
          bind(C, name="retrieve_global_id_node")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (commRank /= 0) then
             print *, "Id node is a global array and is supposed to be handled by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (.not. allocated(IdNode)) then
             print *, "Global id node is not allocated."
             print *, "It is possible that the mesh is already distributed."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(IdNode(1))
          cpp_array%n = size(IdNode)

       end subroutine retrieve_global_id_node

       subroutine retrieve_global_mesh_connectivity(connectivity) &
          bind(C, name="retrieve_global_mesh_connectivity")

          type(cpp_MeshConnectivity), intent(inout) :: connectivity

          call retrieve_coc(NodebyCell, connectivity%NodebyCell)
          call retrieve_coc(NodebyFace, connectivity%NodebyFace)
          call retrieve_coc(FacebyNode, connectivity%FacebyNode)
          call retrieve_coc(FacebyCell, connectivity%FacebyCell)
          call retrieve_coc(CellbyNode, connectivity%CellbyNode)
          call retrieve_coc(CellbyFace, connectivity%CellbyFace)
          call retrieve_coc(CellbyCell, connectivity%CellbyCell)

       end subroutine retrieve_global_mesh_connectivity

       subroutine retrieve_mesh_connectivity(connectivity) &
          bind(C, name="retrieve_mesh_connectivity")

          type(cpp_MeshConnectivity), intent(inout) :: connectivity
          type(CSR) :: empty_CSR
          empty_CSR%Nb = 0

          call retrieve_coc(NodebyCellLocal, connectivity%NodebyCell)
          call retrieve_coc(NodebyFaceLocal, connectivity%NodebyFace)
          call retrieve_coc(FacebyCellLocal, connectivity%FacebyCell)
          ! FIXME: Use a local connectivity structure
          ! The following connectivity elements are not defined locally
          call retrieve_coc(empty_CSR, connectivity%FacebyNode)
          call retrieve_coc(empty_CSR, connectivity%CellbyNode)
          call retrieve_coc(empty_CSR, connectivity%CellbyFace)
          call retrieve_coc(empty_CSR, connectivity%CellbyCell)

       end subroutine retrieve_mesh_connectivity

       subroutine retrieve_nodes_by_fractures(coc) &
          bind(C, name="retrieve_nodes_by_fractures")

          type(cpp_COC), intent(inout) :: coc

          call retrieve_coc(NodebyFractureLocal, coc)

       end subroutine retrieve_nodes_by_fractures

       function check_mesh_allocation() result(status)

          logical :: status

          status = .false.

          if (commRank /= 0) then
             print *, "check_mesh_allocation is supposed to be run by master process"
             return
          end if

          if (.not. allocated(XNode)) then
             print *, "Node array is not allocated."
             return
          end if

          if (size(XNode, 2) /= NbNode) then
             print *, "Node array is not consistent with number of nodes."
             return
          end if

          status = .true.

       end function check_mesh_allocation

       subroutine DefWell_make_compute_well_index_from_C() &
          bind(C, name="DefWell_make_compute_well_index")

          if (.not. check_mesh_allocation()) then
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          call DefWell_Make_ComputeWellIndex( &
             NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
             PermCell, PermFrac)

          call DefMSWell_Make_ComputeWellIndex( &
             NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
             PermCell, PermFrac)

       end subroutine DefWell_make_compute_well_index_from_C

       subroutine GlobalMesh_create_mesh_from_C(nbnodes, nbcells, nbfaces, &
                                                nodes, &
                                                cell_faces_ptr, cell_faces_val, &
                                                cell_nodes_ptr, cell_nodes_val, &
                                                face_nodes_ptr, face_nodes_val, &
                                                cell_id, face_id) &
          bind(C, name="GlobalMesh_create_mesh")

          integer(c_int), value, intent(in) :: nbnodes
          integer(c_int), value, intent(in) :: nbcells
          integer(c_int), value, intent(in) :: nbfaces
          real(c_double), dimension(3, nbnodes), intent(in) :: nodes
          integer(c_int), dimension(nbcells + 1), intent(in) :: cell_faces_ptr
          integer(c_int), dimension(cell_faces_ptr(nbcells + 1)), intent(in) :: cell_faces_val
          integer(c_int), dimension(nbcells + 1), intent(in) :: cell_nodes_ptr
          integer(c_int), dimension(cell_nodes_ptr(nbcells + 1)), intent(in) :: cell_nodes_val
          integer(c_int), dimension(nbfaces + 1), intent(in) :: face_nodes_ptr
          integer(c_int), dimension(face_nodes_ptr(nbfaces + 1)), intent(in) :: face_nodes_val
          integer(c_int), dimension(nbcells), intent(in) :: cell_id
          integer(c_int), dimension(nbfaces), intent(in) :: face_id

          call GlobalMesh_create_mesh(nodes, &
                                      cell_faces_ptr, cell_faces_val, &
                                      cell_nodes_ptr, cell_nodes_val, &
                                      face_nodes_ptr, face_nodes_val, &
                                      cell_id, face_id, &
                                      .true.)

       end subroutine GlobalMesh_create_mesh_from_C

    end module GlobalMeshWrapper
