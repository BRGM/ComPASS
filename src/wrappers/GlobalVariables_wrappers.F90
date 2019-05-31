!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module GlobalVariablesWrappers

       use, intrinsic :: iso_c_binding

       use DefModel, only: NbIncTotalMax, NbIncTotalPrimMax
       use Physics, only: Thickness, gravity, CpRoche, atm_pressure
       use SchemeParameters, only: TimeStepMax, TimeStepInit, TimeFinal
       use NN, only: Delta_t, TimeCurrent

       implicit none

       public :: &
          nb_primary_variables, &
          get_current_time, &
          set_current_time, &
          get_delta_t, &
          get_final_time, &
          set_final_time, &
          get_gravity, &
          set_gravity, &
          get_atm_pressure, &
          set_atm_pressure, &
          get_volumetric_heat_capacity, &
          set_volumetric_heat_capacity, &
          get_fracture_thickness, &
          set_fracture_thickness, &
          get_initial_timestep, &
          get_maximum_timestep, &
          set_initial_timestep, &
          set_maximum_timestep

    contains

        function nb_primary_variables() result(n) &
        bind(C, name="nb_primary_variables")
        integer(c_size_t) :: n
        n = NbIncTotalPrimMax ! NbIncTotalMax
        end function nb_primary_variables
        
        function get_current_time() result(t) &
          bind(C, name="get_current_time")
          real(c_double) :: t
          t = TimeCurrent
       end function get_current_time

       subroutine set_current_time(t) &
          bind(C, name="set_current_time")
          real(c_double), intent(in), value :: t
          TimeCurrent = t
       end subroutine set_current_time

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
          bind(C, name="get_gravity")
          real(c_double) :: g
          g = gravity
       end function get_gravity

       subroutine set_gravity(g) &
          bind(C, name="set_gravity")
          real(c_double), value, intent(in) :: g
          gravity = g
       end subroutine set_gravity

       function get_atm_pressure() result(p) &
         bind(C, name="get_atm_pressure")
         real(c_double) :: p
         p = atm_pressure
      end function get_atm_pressure

       subroutine set_atm_pressure(p) &
          bind(C, name="set_atm_pressure")
          real(c_double), value, intent(in) :: p
          atm_pressure = p
       end subroutine set_atm_pressure

       function get_volumetric_heat_capacity() result(cp) &
          bind(C, name="get_rock_volumetric_heat_capacity")
          real(c_double) :: cp
          cp = CpRoche
       end function get_volumetric_heat_capacity

       subroutine set_volumetric_heat_capacity(cp) &
          bind(C, name="set_rock_volumetric_heat_capacity")
          real(c_double), value, intent(in) :: cp
          CpRoche = cp
       end subroutine set_volumetric_heat_capacity

       function get_fracture_thickness() result(d) &
          bind(C, name="get_fracture_thickness")
          real(c_double) :: d
          d = Thickness
       end function get_fracture_thickness

       subroutine set_fracture_thickness(d) &
          bind(C, name="set_fracture_thickness")
          real(c_double), value, intent(in) :: d
          Thickness = d
       end subroutine set_fracture_thickness

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

