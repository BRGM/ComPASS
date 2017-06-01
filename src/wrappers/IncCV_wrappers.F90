
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
          retrieve_boundary_node_states, &
          retrieve_node_states, &
          retrieve_fracture_states, &
          retrieve_cell_states

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

       subroutine retrieve_boundary_node_states(cpp_array) &
          bind(C, name="retrieve_boundary_node_states")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          if (.not. allocated(IncNodeDirBC)) then
             cpp_array%p = C_NULL_PTR
             cpp_array%n = 0
          else
             cpp_array%p = c_loc(IncNodeDirBC(1))
             cpp_array%n = size(IncNodeDirBC)
          end if
       end subroutine retrieve_boundary_node_states

       subroutine retrieve_node_states(cpp_array) &
          bind(C, name="retrieve_node_states")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          if (.not. allocated(IncCell)) then
             cpp_array%p = C_NULL_PTR
             cpp_array%n = 0
          else
             cpp_array%p = c_loc(IncCell(1))
             cpp_array%n = size(IncCell)
          end if
       end subroutine retrieve_node_states

       subroutine retrieve_fracture_states(cpp_array) &
          bind(C, name="retrieve_fracture_states")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          if (.not. allocated(IncFrac)) then
             cpp_array%p = C_NULL_PTR
             cpp_array%n = 0
          else
             cpp_array%p = c_loc(IncFrac(1))
             cpp_array%n = size(IncFrac)
          end if
       end subroutine retrieve_fracture_states

       subroutine retrieve_cell_states(cpp_array) &
          bind(C, name="retrieve_cell_states")
          type(cpp_array_wrapper), intent(inout) :: cpp_array
          if (.not. allocated(IncCell)) then
             cpp_array%p = C_NULL_PTR
             cpp_array%n = 0
          else
             cpp_array%p = c_loc(IncCell(1))
             cpp_array%n = size(IncCell)
          end if
       end subroutine retrieve_cell_states

    end module IncCVWrapper

