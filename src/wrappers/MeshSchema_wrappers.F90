
    module MeshSchema_wrappers

       use, intrinsic :: iso_c_binding

       use CommonTypesWrapper
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_id_node, &
          retrieve_frac_face_id, &
          retrieve_face_frac_id

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

       subroutine retrieve_frac_face_id(cpp_array) &
          bind(C, name="retrieve_frac_face_id")

          type(cpp_array_wrapper), intent(inout) :: cpp_array
          integer :: n
          
          n = size(FracToFaceLocal)
          cpp_array%n = n

          if (n>0 .and. (.not. allocated(FracToFaceLocal))) then
             print *, "Local frac face id is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if(n==0) then
              cpp_array%p = C_NULL_PTR
          else 
              cpp_array%p = c_loc(FracToFaceLocal(1))
          end if
          
       end subroutine retrieve_frac_face_id

       subroutine retrieve_face_frac_id(cpp_array) &
          bind(C, name="retrieve_face_frac_id")
       
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          integer :: n
          
          n = size(FaceToFracLocal)
          cpp_array%n = n
       
          if (n==0 .or. (.not. allocated(FaceToFracLocal))) then
             print *, "Local frac face id is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if
       
          if(n==0) then 
              cpp_array%p = C_NULL_PTR
          else 
              cpp_array%p = c_loc(FaceToFracLocal(1))
          end if
          
       end subroutine retrieve_face_frac_id

    end module MeshSchema_wrappers

