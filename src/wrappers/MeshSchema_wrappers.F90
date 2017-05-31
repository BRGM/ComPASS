
    module MeshSchema_wrappers

       use, intrinsic :: iso_c_binding

       use CommonTypesWrapper
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_id_node

    contains

       subroutine retrieve_id_node(cpp_array) &
          bind(C, name="retrieve_id_node")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(IdNodeLocal)) then
             print *, "Local id node is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(IdNodeLocal(1))
          cpp_array%n = size(IdNodeLocal)

       end subroutine retrieve_id_node

    end module MeshSchema_wrappers

