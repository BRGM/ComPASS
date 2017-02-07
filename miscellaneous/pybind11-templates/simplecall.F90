
    module example_module

    use, intrinsic :: iso_c_binding

    implicit none

    contains

    
    function my_fortran_add(a, b) result(S) bind(C, name="my_fortran_add")

    integer(c_int), value, intent(in)   :: a, b
    integer(c_int)                      :: S

    write(*, *) 'Hello from Fortran!'
    S = a + b
    write(*,*) 'The Fortran sum: ', a, ' + ', b, ' evaluates to ', S

    end function my_fortran_add

    
    function my_fortran_add_byref(a, b) result(S) bind(C, name="my_fortran_add_byref")

    integer(c_int), intent(in)          :: a, b
    integer(c_int)                      :: S

    write(*, *) 'Hello from Fortran!'
    S = a + b
    write(*,*) 'The Fortran sum: ', a, ' + ', b, ' evaluates to ', S

    end function my_fortran_add_byref

    end module example_module
