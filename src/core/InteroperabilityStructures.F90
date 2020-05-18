!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module InteroperabilityStructures

   use, intrinsic :: iso_c_binding
   use CommonType, only: CSRArray2dble
   use CommonMPI, only: CommonMPI_abort

   implicit none

   type, bind(C) :: cpp_array_wrapper
      type(c_ptr)       :: p
      integer(c_size_t) :: n
   end type cpp_array_wrapper

   type, bind(C) :: cpp_array_wrapper_dim2
      type(c_ptr)       :: p
      integer(c_size_t) :: ni
      integer(c_size_t) :: nj
   end type cpp_array_wrapper_dim2

   type, bind(C) :: csr_block_matrix_wrapper
      type(c_ptr)    :: row_offset
      type(c_ptr)    :: column
      type(c_ptr)    :: block_data
      integer(c_int) :: nb_rows
      integer(c_int) :: block_size
   end type csr_block_matrix_wrapper

contains

   subroutine retrieve_csr_block_matrix(matrix, c_wrapper)

      type(CSRArray2dble), intent(in) :: matrix
      type(csr_block_matrix_wrapper), intent(inout) :: c_wrapper
      if (size(matrix%Val, 1) /= size(matrix%Val, 2)) &
         call CommonMPI_abort("Expecting square subblocks.")
      if (size(matrix%Val, 3) /= matrix%Pt(matrix%Nb + 1)) &
         call CommonMPI_abort("Inconsistent offset pointers.")

      c_wrapper%row_offset = c_loc(matrix%Pt)
      c_wrapper%column = c_loc(matrix%Num)
      c_wrapper%block_data = c_loc(matrix%Val)
      c_wrapper%nb_rows = matrix%Nb
      c_wrapper%block_size = size(matrix%Val, 1)

   end subroutine retrieve_csr_block_matrix

   subroutine retrieve_id_array(p, c_wrapper)

      integer(c_int), dimension(:), pointer, intent(in) :: p
      type(cpp_array_wrapper), intent(inout) :: c_wrapper
      if (.not. associated(p)) &
         call CommonMPI_abort("Pointer is not associated.")

      c_wrapper%p = c_loc(p(1))
      c_wrapper%n = size(p, 1)

   end subroutine retrieve_id_array

   subroutine retrieve_double_array_dim2(p, a)
      real(c_double), dimension(:, :), pointer, intent(in) :: p
      type(cpp_array_wrapper_dim2), intent(inout) :: a

      if (.not. associated(p)) &
         call CommonMPI_abort("Pointer is not associated.")

      a%p = c_loc(p(1, 1))
      a%ni = size(p, 2) ! F order to C order
      a%nj = size(p, 1)

   end subroutine retrieve_double_array_dim2

end module InteroperabilityStructures
