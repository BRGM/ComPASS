!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module IncCVWrapper

       use, intrinsic :: iso_c_binding
       use CommonMPI, only: commRank, CommonMPI_abort
       use DefModel, only: NbComp, NbPhase
       use IncCVReservoir, only: TYPE_IncCVReservoir, IncNode, IncCell, IncFrac
       use MeshSchema, only: NbNodeOwn_Ncpus, NbCellOwn_Ncpus, NbFracOwn_Ncpus
       use IncCVWells, only: IncPressionWellProd, IncPressionWellInj
       use DirichletContribution, only: IncNodeDirBC
       use InteroperabilityStructures, only: cpp_array_wrapper

       use CommonTypesWrapper, only: retrieve_double_array

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
          type(TYPE_IncCVReservoir) :: dummy
          n = sizeof(dummy)
       end function size_of_unknowns

       subroutine retrieve_state_array(states, cpp_array)

          type(TYPE_IncCVReservoir), allocatable, dimension(:), target, intent(inout) :: states
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: n

          if (.not. allocated(states)) then
              cpp_array%p = C_NULL_PTR
              cpp_array%n = 0
          else
              n = size(states)
              cpp_array%n = n
              if (n==0) then
#ifdef TRACK_ZERO_SIZE_ARRAY
                  ! FIXME: Remove comment
                  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!! Zero size array'
#endif
                  cpp_array%p = C_NULL_PTR
              else
                  cpp_array%p = c_loc(states(1))
              end if
          end if

       end subroutine retrieve_state_array

       subroutine retrieve_pointed_state_array(states, cpp_array)

          type(TYPE_IncCVReservoir), dimension(:), pointer, intent(in) :: states
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: n

          if (.not. associated(states)) then
              cpp_array%p = C_NULL_PTR
              cpp_array%n = 0
          else
              n = size(states)
              cpp_array%n = n
              if (n==0) then
#ifdef TRACK_ZERO_SIZE_ARRAY
                  ! FIXME: Remove comment
                  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!! Zero size array'
#endif
                  cpp_array%p = C_NULL_PTR
              else
                  cpp_array%p = c_loc(states(1))
              end if
          end if

       end subroutine retrieve_pointed_state_array

       subroutine retrieve_dirichlet_node_states(cpp_array) &
          bind(C, name="retrieve_dirichlet_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_state_array(IncNodeDirBC, cpp_array)
       end subroutine retrieve_dirichlet_node_states

       subroutine retrieve_own_dirichlet_node_states(cpp_array) &
          bind(C, name="retrieve_own_dirichlet_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: nb_owns = 0
          call retrieve_state_array(IncNodeDirBC, cpp_array)
          nb_owns = NbNodeOwn_Ncpus(commRank+1)
          if(cpp_array%n<nb_owns) &
            call CommonMPI_abort("inconsistent node sizes")
          cpp_array%n = nb_owns
       end subroutine retrieve_own_dirichlet_node_states

       subroutine retrieve_node_states(cpp_array) &
          bind(C, name="retrieve_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_pointed_state_array(IncNode, cpp_array)
       end subroutine retrieve_node_states

       subroutine retrieve_own_node_states(cpp_array) &
          bind(C, name="retrieve_own_node_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: nb_owns = 0
          call retrieve_pointed_state_array(IncNode, cpp_array)
          nb_owns = NbNodeOwn_Ncpus(commRank+1)
          if(cpp_array%n<nb_owns) &
            call CommonMPI_abort("inconsistent node sizes")
          cpp_array%n = nb_owns
       end subroutine retrieve_own_node_states

       subroutine retrieve_fracture_states(cpp_array) &
          bind(C, name="retrieve_fracture_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_pointed_state_array(IncFrac, cpp_array)
       end subroutine retrieve_fracture_states

       subroutine retrieve_own_fracture_states(cpp_array) &
          bind(C, name="retrieve_own_fracture_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: nb_owns = 0
          call retrieve_pointed_state_array(IncFrac, cpp_array)
          nb_owns = NbFracOwn_Ncpus(commRank+1)
          if(cpp_array%n<nb_owns) &
            call CommonMPI_abort("inconsistent node sizes")
          cpp_array%n = nb_owns
       end subroutine retrieve_own_fracture_states

       subroutine retrieve_cell_states(cpp_array) &
          bind(C, name="retrieve_cell_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          call retrieve_pointed_state_array(IncCell, cpp_array)
       end subroutine retrieve_cell_states

       subroutine retrieve_own_cell_states(cpp_array) &
          bind(C, name="retrieve_own_cell_states")
          type(cpp_array_wrapper), intent(out) :: cpp_array
          integer(c_size_t) :: nb_owns = 0
          call retrieve_pointed_state_array(IncCell, cpp_array)
          nb_owns = NbCellOwn_Ncpus(commRank+1)
          if(cpp_array%n<nb_owns) &
            call CommonMPI_abort("inconsistent node sizes")
          cpp_array%n = nb_owns
       end subroutine retrieve_own_cell_states

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

