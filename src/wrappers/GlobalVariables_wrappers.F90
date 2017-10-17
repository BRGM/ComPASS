!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module GlobalVariablesWrappers

       use, intrinsic :: iso_c_binding

       use DefModel
       use NN

       implicit none

       public :: &
          get_current_time, &
          get_delta_t, &
          get_final_time, &
          set_final_time, &
          get_gravity, &
          get_initial_timestep, &
          get_maximum_timestep, &
          set_initial_timestep, &
          set_maximum_timestep

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

       function get_gravity() result(g) &
          bind(C, name="gravity")
          real(c_double) :: g
          g = Gravite
       end function get_gravity

       function get_initial_timestep() result(t) &
          bind(C, name="get_initial_timestep")
          real(c_double) :: t
          t = TimeStepInit
       end function get_initial_timestep

       function get_maximum_timestep() result(t) &
          bind(C, name="get_maximum_timestep")
          real(c_double) :: t
          t = TimeStepMax
       end function get_maximum_timestep

       subroutine set_initial_timestep(t) &
          bind(C, name="set_initial_timestep")
          real(c_double), value, intent(in) :: t
          TimeStepInit = t
       end subroutine set_initial_timestep

       subroutine set_maximum_timestep(t) &
          bind(C, name="set_maximum_timestep")
          real(c_double), value, intent(in) :: t
          TimeStepMax = t
       end subroutine set_maximum_timestep

    end module GlobalVariablesWrappers

