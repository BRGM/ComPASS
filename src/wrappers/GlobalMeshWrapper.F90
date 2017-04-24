
    module GlobalMeshWrapper

       use, intrinsic :: iso_c_binding

       use CommonTypesWrapper
       use CommonMPI
       use GlobalMesh
       use DefWell

       implicit none

       integer :: Ierr, errcode

       type, bind(C) :: cpp_MeshConnectivity
          type(cpp_COC) :: NodebyCell
          type(cpp_COC) :: NodebyFace; 
          type(cpp_COC) :: FacebyCell; 
          type(cpp_COC) :: CellbyNode; 
          type(cpp_COC) :: CellbyFace; 
          type(cpp_COC) :: CellbyCell; 
       end type cpp_MeshConnectivity

       public :: &
          check_mesh_allocation, &
          retrieve_vertices, &
          retrieve_id_faces, &
          retrieve_mesh_connectivity, &
          retrieve_cell_porosity, &
          retrieve_fracture_porosity, &
          retrieve_cell_permeability, &
          retrieve_fracture_permeability, &
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
          DefWell_make_compute_well_index_from_C

    contains

       subroutine retrieve_vertices(cpp_array) &
          bind(C, name="retrieve_vertices")

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

       end subroutine retrieve_vertices

       subroutine retrieve_id_faces(cpp_array) &
          bind(C, name="retrieve_id_faces")

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

       end subroutine retrieve_id_faces

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

          if (.not. allocated(PorositeFace)) then
             print *, "Fracture porosity array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PorositeFace(1))
          cpp_array%n = size(PorositeFace)

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

          cpp_array%p = c_loc(PermCell(1,1,1))
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
             print *, "Fracture permeability array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(PermFrac(1))
          cpp_array%n = size(PermFrac)

       end subroutine retrieve_fracture_permeability

       subroutine retrieve_mesh_connectivity(connectivity) bind(C, name="retrieve_mesh_connectivity")

          type(cpp_MeshConnectivity), intent(inout) :: connectivity

          call retrieve_coc(NodebyCell, connectivity%NodebyCell)
          call retrieve_coc(NodebyFace, connectivity%NodebyFace)
          call retrieve_coc(FacebyCell, connectivity%FacebyCell)
          call retrieve_coc(CellbyNode, connectivity%CellbyNode)
          call retrieve_coc(CellbyFace, connectivity%CellbyFace)
          call retrieve_coc(CellbyCell, connectivity%CellbyCell)

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

    end module GlobalMeshWrapper

