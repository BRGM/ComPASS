!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module GlobalVariablesWrappers

       use, intrinsic :: iso_c_binding

       use DefModel, only: NbIncTotalMax, NbIncTotalPrimMax, LIQUID_PHASE
       use Physics, only: Thickness, gravity, CpRoche, atm_pressure, atm_temperature, &
                          atm_flux_radiation, soil_emissivity, rain_flux
       use SchemeParameters, only: TimeStepMax, TimeStepInit, TimeFinal

       implicit none

       public :: &
          nb_primary_variables, &
          get_gravity, &
          set_gravity, &
#ifdef _WIP_FREEFLOW_STRUCTURES_
          get_atm_pressure, &
          set_atm_pressure, &
          get_atm_temperature, &
          set_atm_temperature, &
          get_atm_flux_radiation, &
          set_atm_flux_radiation, &
          get_soil_emissivity, &
          set_soil_emissivity, &
          get_atm_rain_flux, &
          set_atm_rain_flux, &
#endif
          get_volumetric_heat_capacity, &
          set_volumetric_heat_capacity, &
          get_fracture_thickness, &
          set_fracture_thickness

    contains

        function nb_primary_variables() result(n) &
        bind(C, name="nb_primary_variables")
        integer(c_size_t) :: n
        n = NbIncTotalPrimMax ! NbIncTotalMax
        end function nb_primary_variables

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

#ifdef _WIP_FREEFLOW_STRUCTURES_
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

       function get_atm_temperature() result(T) &
         bind(C, name="get_atm_temperature")
         real(c_double) :: T
         T = atm_temperature
      end function get_atm_temperature

       subroutine set_atm_temperature(T) &
          bind(C, name="set_atm_temperature")
          real(c_double), value, intent(in) :: T
          atm_temperature = T
       end subroutine set_atm_temperature

       function get_atm_flux_radiation() result(q) &
         bind(C, name="get_atm_flux_radiation")
         real(c_double) :: q
         q = atm_flux_radiation
      end function get_atm_flux_radiation

       subroutine set_atm_flux_radiation(q) &
          bind(C, name="set_atm_flux_radiation")
          real(c_double), value, intent(in) :: q
          atm_flux_radiation = q
       end subroutine set_atm_flux_radiation

       function get_soil_emissivity() result(cst) &
         bind(C, name="get_soil_emissivity")
         real(c_double) :: cst
         cst = soil_emissivity
      end function get_soil_emissivity

       subroutine set_soil_emissivity(cst) &
          bind(C, name="set_soil_emissivity")
          real(c_double), value, intent(in) :: cst
          soil_emissivity = cst
       end subroutine set_soil_emissivity

       function get_atm_rain_flux() result(q_rain) &
         bind(C, name="get_atm_rain_flux")
         real(c_double) :: q_rain
         q_rain = rain_flux(LIQUID_PHASE)
      end function get_atm_rain_flux

       subroutine set_atm_rain_flux(q_rain) &
          bind(C, name="set_atm_rain_flux")
          real(c_double), value, intent(in) :: q_rain
          rain_flux(LIQUID_PHASE) = q_rain
       end subroutine set_atm_rain_flux
#endif

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

    end module GlobalVariablesWrappers

