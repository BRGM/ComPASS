
    module GlobalMeshWrapper

       use, intrinsic :: iso_c_binding

       use CommonType
       use GlobalMesh

       implicit none

       integer :: Ierr, errcode
       
       type, bind(C) :: cpp_array_wrapper
          type(c_ptr)       :: p
          integer(c_size_t) :: n
       end type cpp_array_wrapper

       ! Container of container
       type, bind(C) :: cpp_COC
          integer(c_int)    :: nb_containers
          type(c_ptr)       :: container_offset
          type(c_ptr)       :: container_content
       end type cpp_COC

       type, bind(C) :: cpp_CSR
          integer(c_int)    :: nb_rows
          type(c_ptr)       :: row_begin
          type(c_ptr)       :: row_columns
          type(c_ptr)       :: row_nonzeros
       end type cpp_CSR

       type, bind(C) :: cpp_MeshConnectivity
	    type(cpp_COC) :: NodebyCell
	    type(cpp_COC) :: NodebyFace;
	    type(cpp_COC) :: FacebyCell;
	    type(cpp_COC) :: CellbyNode;
	    type(cpp_COC) :: CellbyFace;
	    type(cpp_COC) :: CellbyCell;
       end type cpp_MeshConnectivity

    protected :: &
            f2c_integer_array_to_pointer, &
            f2c_double_array_to_pointer, &
            f2c_double_array, &
            retrieve_coc

    public :: &
          retrieve_vertices, &
            retrieve_mesh_connectivity


    contains

       subroutine f2c_integer_array_to_pointer(array, p)

         integer(kind=c_int), allocatable, dimension(:), target, intent(in) :: array
         type(c_ptr), intent(inout) :: p

          p = c_loc(array(1))
         
end subroutine f2c_integer_array_to_pointer

       subroutine f2c_double_array_to_pointer(array, p)

         real(kind=c_double), allocatable, dimension(:), target, intent(in) :: array
         type(c_ptr), intent(inout) :: p

          p = c_loc(array(1))
         
end subroutine f2c_double_array_to_pointer

subroutine f2c_double_array(fortran_array, cpp_array)

          real(kind=c_double), allocatable, dimension(:), target, intent(in) :: fortran_array
          type(cpp_array_wrapper), intent(out) :: cpp_array

          if (.NOT. allocated(fortran_array)) then
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          call f2c_double_array_to_pointer(fortran_array, cpp_array%p)
          cpp_array%n = size(fortran_array)

       end subroutine f2c_double_array

       subroutine retrieve_vertices(cpp_array) bind(C, name="retrieve_vertices")

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

         subroutine retrieve_coc(fortran_csr, retrieved_coc)
         
          type(CSR), intent(in) :: fortran_csr
          type(cpp_COC), intent(inout) :: retrieved_coc

          if ((.not. allocated(fortran_csr%Pt)) .or. &
              (.not. allocated(fortran_csr%Num)) .or. &
              allocated(fortran_csr%Val)) &
             then
             print *, "Trying to retrieve as COC a CSR which is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          retrieved_coc%nb_containers = fortran_csr%Nb
          call f2c_integer_array_to_pointer(fortran_csr%Pt, retrieved_coc%container_offset)
          call f2c_integer_array_to_pointer(fortran_csr%Num, retrieved_coc%container_content)

       end subroutine retrieve_coc

       subroutine retrieve_mesh_connectivity(connectivity) bind(C, name="retrieve_mesh_connectivity")

         type(cpp_MeshConnectivity), intent(inout) :: connectivity

         call retrieve_coc(NodebyCell, connectivity%NodebyCell)
	 call retrieve_coc(NodebyFace, connectivity%NodebyFace)
	 call retrieve_coc(FacebyCell, connectivity%FacebyCell)
	 call retrieve_coc(CellbyNode, connectivity%CellbyNode)
	 call retrieve_coc(CellbyFace, connectivity%CellbyFace)
  call retrieve_coc(CellbyCell, connectivity%CellbyCell)
  
        end subroutine retrieve_mesh_connectivity

    end module GlobalMeshWrapper
         
