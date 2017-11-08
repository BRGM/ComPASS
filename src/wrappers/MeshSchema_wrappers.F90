
    module MeshSchema_wrappers

       use, intrinsic :: iso_c_binding

       use CommonMPI
       use CommonTypesWrapper
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_id_node, &
          retrieve_frac_face_id, &
          retrieve_face_frac_id, &
          retrieve_nb_cells_own, &
          retrieve_nb_faces_own, &
          retrieve_nb_nodes_own, &
          retrieve_nb_fractures_own

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

          subroutine check_nb_own_array(array)

          integer(c_int), allocatable, dimension(:), intent(in) :: array

          if( .not. allocated(array) ) then
              print *, "Partition info array is not allocated."
              !CHECKME: MPI_Abort is supposed to end all MPI processes
              call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if( size(array) /= Ncpus ) then
              print *, "Inconsistent partition info array size."
              !CHECKME: MPI_Abort is supposed to end all MPI processes
              call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          end subroutine check_nb_own_array

          
          subroutine retrieve_nb_cells_own(cpp_array) &
              bind(C, name="retrieve_nb_cells_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbCellOwn_Ncpus)
          cpp_array%n = Ncpus
          if(Ncpus==0) then
              cpp_array%p = C_NULL_PTR
          else
              cpp_array%p = c_loc(NbCellOwn_Ncpus(1))
          end if
          end subroutine retrieve_nb_cells_own

          subroutine retrieve_nb_faces_own(cpp_array) &
              bind(C, name="retrieve_nb_faces_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbFaceOwn_Ncpus)
          cpp_array%n = Ncpus
          if(Ncpus==0) then
              cpp_array%p = C_NULL_PTR
          else
              cpp_array%p = c_loc(NbFaceOwn_Ncpus(1))
          end if
          end subroutine retrieve_nb_faces_own

          subroutine retrieve_nb_nodes_own(cpp_array) &
              bind(C, name="retrieve_nb_nodes_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbNodeOwn_Ncpus)
          cpp_array%n = Ncpus
          if(Ncpus==0) then
              cpp_array%p = C_NULL_PTR
          else
              cpp_array%p = c_loc(NbNodeOwn_Ncpus(1))
          end if
          end subroutine retrieve_nb_nodes_own

          subroutine retrieve_nb_fractures_own(cpp_array) &
              bind(C, name="retrieve_nb_fractures_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbFracOwn_Ncpus)
          cpp_array%n = Ncpus
          if(Ncpus==0) then
              cpp_array%p = C_NULL_PTR
          else
              cpp_array%p = c_loc(NbFracOwn_Ncpus(1))
          end if
          end subroutine retrieve_nb_fractures_own

    end module MeshSchema_wrappers

