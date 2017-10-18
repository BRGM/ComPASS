!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module LocalMeshWrapper

       use, intrinsic :: iso_c_binding

       use CommonTypesWrapper
       use CommonMPI
       use MeshSchema

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_vertices, &
          retrieve_nodeflags, &
          retrieve_cellflags, &
          retrieve_faceflags

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

       subroutine retrieve_cellflags(cpp_array) &
          bind(C, name="retrieve_cellflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(CellFlagsLocal)) then
             print *, "Local cell flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellFlagsLocal(1))
          cpp_array%n = size(CellFlagsLocal)

       end subroutine retrieve_cellflags

       subroutine retrieve_faceflags(cpp_array) &
          bind(C, name="retrieve_faceflags")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(FaceFlagsLocal)) then
             print *, "Local face flags are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FaceFlagsLocal(1))
          cpp_array%n = size(FaceFlagsLocal)

       end subroutine retrieve_faceflags

    end module LocalMeshWrapper

