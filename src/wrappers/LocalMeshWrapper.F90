!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LocalMeshWrapper

   use, intrinsic :: iso_c_binding
   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort
   use InteroperabilityStructures, only: cpp_array_wrapper
   use MeshSchema, only: &
      NodeFlagsLocal, CellFlagsLocal, FaceFlagsLocal, &
      CellTypesLocal, FaceTypesLocal, XNodeLocal, &
      NbNodeLocal_Ncpus, NbCellLocal_Ncpus, NbFracLocal_Ncpus

   implicit none

   integer :: Ierr, errcode

   public :: &
      retrieve_vertices, &
      retrieve_nodeflags, &
      retrieve_cellflags, &
      retrieve_faceflags, &
      retrieve_celltypes, &
      retrieve_facetypes

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

   subroutine retrieve_celltypes(cpp_array) &
      bind(C, name="retrieve_celltypes")

      type(cpp_array_wrapper), intent(inout) :: cpp_array

      if (.not. allocated(CellTypesLocal)) then
         print *, "Local cell types are not allocated."
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      cpp_array%p = c_loc(CellTypesLocal(1))
      cpp_array%n = size(CellTypesLocal)

   end subroutine retrieve_celltypes

   subroutine retrieve_facetypes(cpp_array) &
      bind(C, name="retrieve_facetypes")

      type(cpp_array_wrapper), intent(inout) :: cpp_array

      if (.not. allocated(FaceTypesLocal)) then
         print *, "Local face types are not allocated."
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      cpp_array%p = c_loc(FaceTypesLocal(1))
      cpp_array%n = size(FaceTypesLocal)

   end subroutine retrieve_facetypes

   subroutine retrieve_size_if_allocated(n, array, array_name)
      integer(c_size_t), intent(out) :: n
      integer(c_int), allocatable, dimension(:), intent(in) :: array
      character(len=*), intent(in) :: array_name
      if (.not. allocated(array)) then
         call CommonMPI_abort(array_name//" not allocated.")
      endif
      n = array(commRank + 1)
   end subroutine retrieve_size_if_allocated

   function nb_nodes() result(n) &
      bind(C, name="nb_nodes")
      integer(c_size_t) :: n
      call retrieve_size_if_allocated(n, NbNodeLocal_Ncpus, "NbNodeLocal_Ncpus")
   end function nb_nodes

   function nb_fractures() result(n) &
      bind(C, name="nb_fractures")
      integer(c_size_t) :: n
      call retrieve_size_if_allocated(n, NbFracLocal_Ncpus, "NbFracLocal_Ncpus")
   end function nb_fractures

   function nb_cells() result(n) &
      bind(C, name="nb_cells")
      integer(c_size_t) :: n
      call retrieve_size_if_allocated(n, NbCellLocal_Ncpus, "NbCellLocal_Ncpus")
   end function nb_cells

end module LocalMeshWrapper

