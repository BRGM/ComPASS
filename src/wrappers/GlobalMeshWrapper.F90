!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module GlobalMeshWrapper

       use, intrinsic :: iso_c_binding

       use CommonType
       use CommonTypesWrapper
       use CommonMPI
       use GlobalMesh
       use DefWell
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       type, bind(C) :: cpp_MeshConnectivity
          type(cpp_COC) :: NodebyCell
          type(cpp_COC) :: NodebyFace;
          type(cpp_COC) :: FacebyNode;
          type(cpp_COC) :: FacebyCell;
          type(cpp_COC) :: CellbyNode;
          type(cpp_COC) :: CellbyFace;
          type(cpp_COC) :: CellbyCell;
       end type cpp_MeshConnectivity

       public :: &
          GlobalMesh_allocate_id_nodes, &
          GlobalMesh_count_dirichlet_nodes, &
          check_mesh_allocation, &
          retrieve_global_vertices, &
          retrieve_global_nodeflags, &
          retrieve_global_cellflags, &
          retrieve_global_faceflags, &
          GlobalMesh_allocate_rocktype_from_C, &
          retrieve_global_cellrocktype, &
          retrieve_global_fracrocktype, &
          retrieve_global_celltypes, &
          retrieve_global_facetypes, &
          retrieve_global_id_faces, &
          retrieve_global_mesh_connectivity, &
          retrieve_mesh_connectivity, &
          retrieve_cell_porosity, &
          retrieve_fracture_porosity, &
          retrieve_cell_permeability, &
          retrieve_fracture_permeability, &
          retrieve_global_id_node, &
          GlobalMesh_build_cartesian_grid_from_C, &
          GlobalMesh_make_post_read_from_C, &
          GlobalMesh_Make_post_read_fracture_and_dirBC_from_C, &
          GlobalMesh_Make_post_read_set_poroperm_from_C, &
          GlobalMesh_Make_post_read_well_connectivity_and_ip_from_C, &
          GlobalMesh_MeshBoundingBox_from_C, &
          GlobalMesh_Compute_all_connectivies_from_C, &
          GlobalMesh_SetFrac_from_C, &
          GlobalMesh_NodeOfFrac_from_C, &
          GlobalMesh_SetDirBC_from_C, &
          GlobalMesh_FracbyNode_from_C, &
          DefWell_make_compute_well_index_from_C, &
          GlobalMesh_create_mesh_from_C, &
          GlobalMesh_set_cartesian_mesh, &
          GlobalMesh_set_hexahedron_mesh, &
          GlobalMesh_set_tetrahedron_mesh, &
          GlobalMesh_set_wedge_mesh

    contains

    subroutine GlobalMesh_allocate_id_nodes() &
        bind(C, name="GlobalMesh_allocate_id_nodes")

    if(allocated(Idnode)) then
        deallocate(IdNode)
    end if
    allocate(IdNode(NbNode))

        end subroutine GlobalMesh_allocate_id_nodes

    subroutine GlobalMesh_count_dirichlet_nodes() &
        bind(C, name="GlobalMesh_count_dirichlet_nodes")

    integer :: i

    NbDirNodeP = 0
#ifdef _THERMIQUE_
    NbDirNodeT = 0
#endif
    do i=1, NbNode
        if(IdNode(i)%P.eq."d") NbDirNodeP = NbDirNodeP + 1
#ifdef _THERMIQUE_
        if(IdNode(i)%T.eq."d") NbDirNodeT = NbDirNodeT + 1
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


    subroutine retrieve_global_cellrocktype(cpp_array) &
          bind(C, name="retrieve_global_cellrocktype")

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
             print *, "cell rocktype are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellRocktype)
          cpp_array%n = (IndThermique+1)*size(CellRocktype, 2)

       end subroutine retrieve_global_cellrocktype


    subroutine retrieve_global_fracrocktype(cpp_array) &
          bind(C, name="retrieve_global_fracrocktype")

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
             print *, "frac rocktype are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FracRocktype)
          cpp_array%n = (IndThermique+1)*size(FracRocktype, 2)

       end subroutine retrieve_global_fracrocktype


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

       subroutine retrieve_cell_porosity(cpp_array) &
          bind(C, name="retrieve_cell_porosity")

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

       end subroutine retrieve_cell_porosity

       subroutine retrieve_fracture_porosity(cpp_array) &
          bind(C, name="retrieve_fracture_porosity")

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

       end subroutine retrieve_fracture_porosity

       subroutine retrieve_cell_permeability(cpp_array) &
          bind(C, name="retrieve_cell_permeability")

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

       end subroutine retrieve_cell_permeability

       subroutine retrieve_fracture_permeability(cpp_array) &
          bind(C, name="retrieve_fracture_permeability")

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

       end subroutine retrieve_fracture_permeability

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

       subroutine GlobalMesh_build_cartesian_grid_from_C(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz) &
          bind(C, name="GlobalMesh_build_cartesian_grid")

          real(kind=c_double), value, intent(in)  :: Ox, Oy, Oz
          real(kind=c_double), value, intent(in)  :: lx, ly, lz
          integer(kind=c_int), value, intent(in)  :: nx, ny, nz
          ! FIXME: set consistent values to error codes
          integer :: errcode, Ierr

          if (commRank /= 0) then
             print *, "Mesh is supposed to be created by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          call GlobalMesh_Build_cartesian_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)

       end subroutine GlobalMesh_build_cartesian_grid_from_C

       subroutine GlobalMesh_make_post_read_from_C() &
          bind(C, name="GlobalMesh_make_post_read")

          if (commRank /= 0) then
             print *, "GlobalMesh_make_post_read_from_C is supposed to be run by master process."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          call GlobalMesh_Make_post_read()

       end subroutine GlobalMesh_make_post_read_from_C

       subroutine GlobalMesh_Make_post_read_fracture_and_dirBC_from_C() &
          bind(C, name="GlobalMesh_make_post_read_fracture_and_dirBC")

          call GlobalMesh_Make_post_read_fracture_and_dirBC

       end subroutine GlobalMesh_Make_post_read_fracture_and_dirBC_from_C

       subroutine GlobalMesh_Make_post_read_set_poroperm_from_C() &
          bind(C, name="GlobalMesh_make_post_read_set_poroperm")

          call GlobalMesh_Make_post_read_set_poroperm

       end subroutine GlobalMesh_Make_post_read_set_poroperm_from_C

       subroutine GlobalMesh_Make_post_read_well_connectivity_and_ip_from_C() &
          bind(C, name="GlobalMesh_make_post_read_well_connectivity_and_ip")

          call GlobalMesh_Make_post_read_well_connectivity_and_ip

       end subroutine GlobalMesh_Make_post_read_well_connectivity_and_ip_from_C

       subroutine GlobalMesh_MeshBoundingBox_from_C() &
          bind(C, name="GlobalMesh_mesh_bounding_box")
          call GlobalMesh_MeshBoundingBox
       end subroutine GlobalMesh_MeshBoundingBox_from_C

       subroutine GlobalMesh_Compute_all_connectivies_from_C() &
          bind(C, name="GlobalMesh_compute_all_connectivies")
          call GlobalMesh_Compute_all_connectivies
       end subroutine GlobalMesh_Compute_all_connectivies_from_C

       subroutine GlobalMesh_SetFrac_from_C() &
          bind(C, name="GlobalMesh_set_frac")
          call GlobalMesh_SetFrac
       end subroutine GlobalMesh_SetFrac_from_C

       subroutine GlobalMesh_NodeOfFrac_from_C() &
          bind(C, name="GlobalMesh_node_of_frac")
          call GlobalMesh_NodeOfFrac
       end subroutine GlobalMesh_NodeOfFrac_from_C

       subroutine GlobalMesh_SetDirBC_from_C() &
          bind(C, name="GlobalMesh_set_dir_BC")
          call GlobalMesh_SetDirBC
       end subroutine GlobalMesh_SetDirBC_from_C

       subroutine GlobalMesh_FracbyNode_from_C() &
          bind(C, name="GlobalMesh_frac_by_node")
          call GlobalMesh_FracbyNode
       end subroutine GlobalMesh_FracbyNode_from_C

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

          subroutine GlobalMesh_set_cartesian_mesh() &
              bind(C, name="GlobalMesh_set_cartesian_mesh")
          MESH_TYPE = "cartesian-quad"
          end subroutine GlobalMesh_set_cartesian_mesh

          subroutine GlobalMesh_set_hexahedron_mesh() &
              bind(C, name="GlobalMesh_set_hexahedron_mesh")
          MESH_TYPE = "hexahedron-quad"
          end subroutine GlobalMesh_set_hexahedron_mesh

          subroutine GlobalMesh_set_tetrahedron_mesh() &
              bind(C, name="GlobalMesh_set_tetrahedron_mesh")
          MESH_TYPE = "tetrahedron-triangle"
          end subroutine GlobalMesh_set_tetrahedron_mesh

          subroutine GlobalMesh_set_wedge_mesh() &
              bind(C, name="GlobalMesh_set_wedge_mesh")
          MESH_TYPE = "wedge"
          end subroutine GlobalMesh_set_wedge_mesh

    end module GlobalMeshWrapper

