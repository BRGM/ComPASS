!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

    module CommonTypesWrapper

       use, intrinsic :: iso_c_binding

       use CommonMPI
       use CommonType
       use InteroperabilityStructures
       
       implicit none

       !type, bind(C) :: cpp_narray_wrapper
       !   integer(c_int) :: dim
       !   type(c_ptr)    :: shape
       !   type(c_ptr)    :: data
       !end type cpp_narray_wrapper

       type, bind(C) :: cpp_narray_wrapper
          type(c_ptr)    :: shape
          type(c_ptr)    :: data
       end type cpp_narray_wrapper

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

       public :: &
          f2c_integer_array_to_pointer, &
          f2c_double_array_to_pointer, &
          f2c_double_array, &
          retrieve_double_array, &
          retrieve_coc, &
          bind_2array, &
          bind_3array

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
          ! FIXME: set consistent values to error codes
          integer :: errcode, Ierr

          if (.NOT. allocated(fortran_array)) then
             !CHECKME: MPI_Abort is supposed to end all MPI processes
             call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
          end if

          call f2c_double_array_to_pointer(fortran_array, cpp_array%p)
          cpp_array%n = size(fortran_array)

       end subroutine f2c_double_array

       subroutine retrieve_double_array(fortran_array, cpp_array)

          real(kind=c_double), allocatable, dimension(:), target, intent(in) :: fortran_array
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: n

          if (.not. allocated(fortran_array)) then
             cpp_array%p = C_NULL_PTR
             cpp_array%n = 0
          else
             n = size(fortran_array)
             cpp_array%n = n
             if (n == 0) then
                ! FIXME: Remove comment
#ifdef TRACK_ZERO_SIZE_ARRAY
                write (*, *) '!!!!!!!!!!!!!!!!!!!!!!! Zero size array'
#endif
                cpp_array%p = C_NULL_PTR
             else
                cpp_array%p = c_loc(fortran_array(1))
             end if
          end if

       end subroutine retrieve_double_array

       subroutine retrieve_coc(fortran_csr, retrieved_coc)

          type(CSR), intent(in) :: fortran_csr
          type(cpp_COC), intent(inout) :: retrieved_coc
          ! FIXME: set consistent values to error codes
          integer :: n, errcode, Ierr

          n = fortran_csr%Nb
          retrieved_coc%nb_containers = n

          if (n == 0) then
#ifdef TRACK_ZERO_SIZE_ARRAY
             write (*, *) 'WARNING - Retrieving zero size COC.'
#endif
             retrieved_coc%container_offset = C_NULL_PTR
             retrieved_coc%container_content = C_NULL_PTR
          else
             if ((.not. allocated(fortran_csr%Pt)) .or. &
                 (.not. allocated(fortran_csr%Num))) then
                print *, "Trying to retrieve as COC a CSR which is not allocated."
                !CHECKME: MPI_Abort is supposed to end all MPI processes
                call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
             end if
             if (allocated(fortran_csr%Val)) then
                print *, "Trying to retrieve as COC a CSR which has allocated values (i.e. it is not a COC but rather a true CSR)."
                !CHECKME: MPI_Abort is supposed to end all MPI processes
                call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
             end if
             call f2c_integer_array_to_pointer(fortran_csr%Pt, retrieved_coc%container_offset)
             call f2c_integer_array_to_pointer(fortran_csr%Num, retrieved_coc%container_content)
          end if

       end subroutine retrieve_coc

       subroutine bind_2array(a, p)

          type(cpp_narray_wrapper), intent(in) :: a
          integer(c_size_t), pointer :: ashape(:)
          real(c_double), pointer, intent(out) :: p(:, :)

          call c_f_pointer(a%shape, ashape, shape=[2])
          call c_f_pointer(a%data, p, shape=ashape(2:1:-1)) ! switch to Fortran column major mode

       end subroutine bind_2array

       subroutine bind_3array(a, p)

          type(cpp_narray_wrapper), intent(in) :: a
          integer(c_size_t), pointer :: ashape(:)
          real(c_double), pointer, intent(out) :: p(:, :, :)

          call c_f_pointer(a%shape, ashape, shape=[3])
          call c_f_pointer(a%data, p, shape=ashape(3:1:-1)) ! switch to Fortran column major mode

       end subroutine bind_3array

    end module CommonTypesWrapper

