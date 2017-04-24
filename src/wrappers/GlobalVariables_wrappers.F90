
    module GlobalVariablesWrappers

       use, intrinsic :: iso_c_binding

       use DefModel
       use NN

       implicit none

       public :: &
          get_current_time, &
          get_delta_t, &
          get_final_time, &
          set_final_time

    contains

       function get_current_time() result(t) &
          bind(C, name="get_current_time")
       real(c_double) :: t
       t = TimeCurrent
       end function get_current_time

       function get_delta_t() result(t) &
          bind(C, name="get_delta_t")
       real(c_double) :: t
       t = Delta_t
       end function get_delta_t

       function get_final_time() result(t) &
          bind(C, name="get_final_time")
       real(c_double) :: t
       t = TimeFinal
       end function get_final_time

       subroutine set_final_time(t) &
          bind(C, name="set_final_time")
       real(c_double), value, intent(in) :: t
       TimeFinal = t
       end subroutine set_final_time

    end module GlobalVariablesWrappers

