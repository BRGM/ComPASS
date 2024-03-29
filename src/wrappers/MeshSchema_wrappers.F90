!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

    module MeshSchema_wrappers

       use, intrinsic :: iso_c_binding
       use mpi, only: MPI_Abort
       use CommonMPI, only: Ncpus, ComPASS_COMM_WORLD
       use InteroperabilityStructures, only: cpp_array_wrapper
       use CommonTypesWrapper, only: cpp_COC, retrieve_dof_family
       use MeshSchema, only: &
          IdNodeLocal, &
          FaceToFracLocal, FracToFaceLocal, &
          NodeDarcyRocktypesLocal, CellDarcyRocktypesLocal, FracDarcyRocktypesLocal, &
          NbNodeOwn_Ncpus, NbCellOwn_Ncpus, NbFaceOwn_Ncpus, NbFracOwn_Ncpus, &
          NbWellInjOwn_Ncpus, NbWellProdOwn_Ncpus, &
          XCellLocal, XFaceLocal, &
          NumNodebyProc, NumFracbyProc, NumWellInjbyProc, NumWellProdbyProc, NumMSWellNodebyProc

       implicit none

       integer :: Ierr, errcode

       public :: &
          retrieve_id_node, &
          retrieve_frac_face_id, &
          retrieve_face_frac_id, &
          retrieve_nb_cells_own, &
          retrieve_nb_faces_own, &
          retrieve_nb_nodes_own, &
          retrieve_nb_fractures_own, &
          retrieve_cell_centers, &
          retrieve_face_centers, &
          retrieve_cell_rocktypes, &
          retrieve_node_rocktypes, &
          retrieve_fracture_rocktypes

    contains

       subroutine retrieve_cell_rocktypes(cpp_array) &
          bind(C, name="retrieve_cell_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. associated(CellDarcyRocktypesLocal)) then
             print *, "cell rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(CellDarcyRocktypesLocal(1))
          cpp_array%n = size(CellDarcyRocktypesLocal)

       end subroutine retrieve_cell_rocktypes

       subroutine retrieve_node_rocktypes(cpp_array) &
          bind(C, name="retrieve_node_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. associated(NodeDarcyRocktypesLocal)) then
             print *, "node rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(NodeDarcyRocktypesLocal(1))
          cpp_array%n = size(NodeDarcyRocktypesLocal)

       end subroutine retrieve_node_rocktypes

       subroutine retrieve_fracture_rocktypes(cpp_array) &
          bind(C, name="retrieve_fracture_rocktypes")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. associated(FracDarcyRocktypesLocal)) then
             print *, "frac rocktypes are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(FracDarcyRocktypesLocal(1))
          cpp_array%n = size(FracDarcyRocktypesLocal)

       end subroutine retrieve_fracture_rocktypes

       subroutine retrieve_cell_centers(cpp_array) &
          bind(C, name="retrieve_cell_centers")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(XCellLocal)) then
             print *, "Local cell centers are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(XCellLocal(1, 1))
          cpp_array%n = size(XCellLocal, 2)

       end subroutine retrieve_cell_centers

       subroutine retrieve_face_centers(cpp_array) &
          bind(C, name="retrieve_face_centers")

          type(cpp_array_wrapper), intent(inout) :: cpp_array

          if (.not. allocated(XFaceLocal)) then
             print *, "Local face centers are not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          cpp_array%p = c_loc(XFaceLocal(1, 1))
          cpp_array%n = size(XFaceLocal, 2)

       end subroutine retrieve_face_centers

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

          if (n > 0 .and. (.not. allocated(FracToFaceLocal))) then
             print *, "Local frac face id is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (n == 0) then
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

          if (n == 0 .or. (.not. allocated(FaceToFracLocal))) then
             print *, "Local frac face id is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (n == 0) then
             cpp_array%p = C_NULL_PTR
          else
             cpp_array%p = c_loc(FaceToFracLocal(1))
          end if

       end subroutine retrieve_face_frac_id

       subroutine check_nb_own_array(array)

          integer(c_int), allocatable, dimension(:), intent(in) :: array

          if (.not. allocated(array)) then
             print *, "Partition info array is not allocated."
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          if (size(array) /= Ncpus) then
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
          if (Ncpus == 0) then
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
          if (Ncpus == 0) then
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
          if (Ncpus == 0) then
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
          if (Ncpus == 0) then
             cpp_array%p = C_NULL_PTR
          else
             cpp_array%p = c_loc(NbFracOwn_Ncpus(1))
          end if
       end subroutine retrieve_nb_fractures_own

       subroutine retrieve_nb_wellinj_own(cpp_array) &
          bind(C, name="retrieve_nb_wellinj_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbWellInjOwn_Ncpus)
          cpp_array%n = Ncpus
          if (Ncpus == 0) then
             cpp_array%p = C_NULL_PTR
          else
             cpp_array%p = c_loc(NbWellInjOwn_Ncpus(1))
          end if
       end subroutine retrieve_nb_wellinj_own

       subroutine retrieve_nb_wellprod_own(cpp_array) &
          bind(C, name="retrieve_nb_wellprod_own")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          call check_nb_own_array(NbWellProdOwn_Ncpus)
          cpp_array%n = Ncpus
          if (Ncpus == 0) then
             cpp_array%p = C_NULL_PTR
          else
             cpp_array%p = c_loc(NbWellProdOwn_Ncpus(1))
          end if
       end subroutine retrieve_nb_wellprod_own
       subroutine retrieve_NumNodebyProc(coc) &
          bind(C, name="retrieve_NumNodebyProc")
          type(cpp_COC), intent(inout) :: coc

          call retrieve_dof_family(NumNodebyProc, coc)

       end subroutine retrieve_NumNodebyProc

       subroutine retrieve_NumFracbyProc(coc) &
          bind(C, name="retrieve_NumFracbyProc")
          type(cpp_COC), intent(inout) :: coc

          call retrieve_dof_family(NumFracbyProc, coc)

       end subroutine retrieve_NumFracbyProc

       subroutine retrieve_NumWellProdbyProc(coc) &
          bind(C, name="retrieve_NumWellProdbyProc")
          type(cpp_COC), intent(inout) :: coc

          call retrieve_dof_family(NumWellProdbyProc, coc)

       end subroutine retrieve_NumWellProdbyProc

       subroutine retrieve_NumWellInjbyProc(coc) &
          bind(C, name="retrieve_NumWellInjbyProc")
          type(cpp_COC), intent(inout) :: coc

          call retrieve_dof_family(NumWellInjbyProc, coc)

       end subroutine retrieve_NumWellInjbyProc

       subroutine retrieve_NumMSWellNodebyProc(coc) &
          bind(C, name="retrieve_NumMSWellNodebyProc")
          type(cpp_COC), intent(inout) :: coc

          call retrieve_dof_family(NumMSWellNodebyProc, coc)

       end subroutine retrieve_NumMSWellNodebyProc

    end module MeshSchema_wrappers
