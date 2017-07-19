
    module IncCVWrapper

       use, intrinsic :: iso_c_binding

       use DefModel
       use IncCV
       use CommonTypesWrapper

       implicit none

       public :: &
          get_NbComp, &
          get_NbPhase, &
          size_of_unknowns, &
          retrieve_dirichlet_node_states, &
          retrieve_node_states, &
          retrieve_fracture_states, &
          retrieve_cell_states, &
          retrieve_injection_whp, &
          retrieve_production_whp

    contains

       function get_NbComp() result(n) &
          bind(C, name="number_of_components")
          integer(c_int) :: n
          n = NbComp
       end function get_NbComp

       function get_NbPhase() result(n) &
          bind(C, name="number_of_phases")
          integer(c_int) :: n
          n = NbPhase
       end function get_NbPhase

       function size_of_unknowns() result(n) &
          bind(C, name="size_of_unknowns")
          integer(c_int) :: n
          type(TYPE_IncCV) :: dummy
          n = sizeof(dummy)
       end function size_of_unknowns

       subroutine retrieve_state_array(states, cpp_array)

          TYPE(TYPE_IncCV), allocatable, dimension(:), target, intent(inout) :: states
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: n

          if (.not. allocated(states)) then
              cpp_array%p = C_NULL_PTR
              cpp_array%n = 0
          else
              n = size(states)
              cpp_array%n = n
              if (n==0) then
#ifndef NDEBUG              
                  ! FIXME: Remove comment
                  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!! Zero size array'
#endif !NDEBUG              
                  cpp_array%p = C_NULL_PTR
              else
                  cpp_array%p = c_loc(states(1))
              end if
          end if

       end subroutine retrieve_state_array

       subroutine retrieve_dirichlet_node_states(cpp_array) &
          bind(C, name="retrieve_dirichlet_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_state_array(IncNodeDirBC, cpp_array)
       end subroutine retrieve_dirichlet_node_states

       subroutine retrieve_node_states(cpp_array) &
          bind(C, name="retrieve_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_state_array(IncNode, cpp_array)
       end subroutine retrieve_node_states

       subroutine retrieve_fracture_states(cpp_array) &
          bind(C, name="retrieve_fracture_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_state_array(IncFrac, cpp_array)
       end subroutine retrieve_fracture_states

       subroutine retrieve_cell_states(cpp_array) &
          bind(C, name="retrieve_cell_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_state_array(IncCell, cpp_array)
       end subroutine retrieve_cell_states

       subroutine retrieve_injection_whp(cpp_array) &
          bind(C, name="retrieve_injection_whp")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_double_array(IncPressionWellInj, cpp_array)
       end subroutine retrieve_injection_whp

       subroutine retrieve_production_whp(cpp_array) &
          bind(C, name="retrieve_production_whp")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_double_array(IncPressionWellProd, cpp_array)
       end subroutine retrieve_production_whp

    end module IncCVWrapper

