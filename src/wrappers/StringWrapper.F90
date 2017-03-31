
    module StringWrapper

      use, intrinsic :: iso_c_binding

      implicit none

      type, bind(C) :: cpp_string_wrapper
        type(c_ptr)       :: p
        integer(c_size_t) :: length
      end type cpp_string_wrapper

      public :: &
        fortran_string

    contains

      function fortran_string(wrapper) result(s)

        type(cpp_string_wrapper), intent(in) :: wrapper
        character(kind=c_char, len=wrapper%length), pointer :: s

        call c_f_pointer(wrapper%p, s)

      end function fortran_string

    end module StringWrapper
