
    module GlobalMeshWrapper

       use, intrinsic :: iso_c_binding

       use GlobalMesh

       implicit none

       type, bind(C) :: cpp_vertices
          type(c_ptr)       :: p
          integer(c_size_t) :: n
       end type cpp_vertices

       public :: &
          vertices_buffer

    contains

       subroutine vertices_buffer(buffer) bind(C, name="vertices_buffer")

          type(cpp_vertices), intent(inout) :: buffer

          if (commRank == 0) then
             buffer%p = c_loc(XNode(1, 1))
             buffer%n = size(XNode, 2)
          else
             !CHECKME: Maybe MPI_abort would be better here
             print *, "Mesh is supposed to be read by master process."
             buffer%p = c_null_ptr
             buffer%n = 0
          end if

       end subroutine vertices_buffer

    end module GlobalMeshWrapper
