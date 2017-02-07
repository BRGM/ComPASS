
    module passingstring_module

    use, intrinsic :: iso_c_binding

    implicit none

    type, bind(C) :: cpp_string_wrapper
        type(c_ptr)       :: p
        integer(c_size_t) :: length
    end type cpp_string_wrapper

    contains

    function fortran_string(wrapper) result(s)

    type(cpp_string_wrapper), intent(in) :: wrapper
    character(c_char), pointer           :: s(:)

    call c_f_pointer(wrapper%p, s, (/wrapper%length/))
    
    end function fortran_string


    subroutine fortran_dump(cpp_string) bind(C, name="fortran_dump")

    type(cpp_string_wrapper), intent(in) :: cpp_string

    write(*,*) 'Fortran writes: ', fortran_string(cpp_string)

    end subroutine fortran_dump

    end module passingstring_module
