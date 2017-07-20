
    module LocalMeshWrapper

       use, intrinsic :: iso_c_binding

       use CommonTypesWrapper
       use CommonMPI
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_vertices, &
          retrieve_nodeflags

    contains

       subroutine retrieve_vertices(cpp_array) &
          bind(C, name="retrieve_vertices")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(XNodeLocal)) then
             print *, "Local mesh vertices are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(XNodeLocal(1, 1))
          cpp_array%n = size(XNodeLocal, 2)

       end subroutine retrieve_vertices

       subroutine retrieve_nodeflags(cpp_array) &
          bind(C, name="retrieve_nodeflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(NodeFlagsLocal)) then
             print *, "Local node flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(NodeFlagsLocal(1))
          cpp_array%n = size(NodeFlagsLocal)

       end subroutine retrieve_nodeflags

    end module LocalMeshWrapper

