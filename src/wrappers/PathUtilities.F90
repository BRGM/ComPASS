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

      subroutine make_directory(path)

        character(len=*), intent(in)  :: path
        character(len=1, kind=c_char), target     :: path_as_C_string(len(path) + 1)

        write (path_as_C_string, '(A)') trim(path)//C_NULL_CHAR
        call c_make_directory( c_loc(path_as_C_string) )

      end subroutine make_directory

    end module PathUtilities
