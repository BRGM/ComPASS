    module PathUtilities

      use, intrinsic :: iso_c_binding

      implicit none

      interface
        subroutine c_make_directory(s) bind(C)
          use, intrinsic :: iso_c_binding, only:c_ptr
          implicit none
          type(c_ptr), value :: s
        end subroutine C_make_directory
      end interface

      public :: &
        make_directory

    contains

      subroutine make_directory(fortran_path)

        character(len=*), intent(in) :: fortran_path
        character(len=1, kind=c_char), target :: path_as_C_string(len(fortran_path) + 1)
        integer :: i

        do i=1, len(fortran_path)
            path_as_C_string(i) = fortran_path(i:i)
        end do
        path_as_C_string(len(fortran_path)+1) = C_NULL_CHAR
        call c_make_directory( c_loc(path_as_C_string(1)) )

      end subroutine make_directory

    end module PathUtilities
