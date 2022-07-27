!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Physics

   use iso_c_binding, only: c_double
   use DefModel, only: NbComp

   real(c_double) :: gravity = 9.81d0

   ! Far field atm constants
   ! must be compatible with reference pressure !!! (called in f_PartialMolarEnthalpy which does not depend on P)
   real(c_double) :: atm_pressure = 1.d5
   real(c_double) :: atm_temperature = 300.d0
   real(c_double) :: rain_temperature = 300.d0
   real(c_double) :: atm_comp(2, 2) = &
                     RESHAPE((/ &
                             0.999d0, 1.d-3, & ! gas
                             0.01d0, 0.99d0 & ! liquid
                             /), (/2, 2/))
   ! Convective-diffusive boundary layer constants
   real(c_double) :: Hm(2) = RESHAPE((/20.d0/29.d0, &  ! gas Convective-diffusive const of the FreeFlow boundary layer
                                       0.d0 &  ! everything which is not gas is nul
                                       /), (/2/))
   real(c_double) :: HT = 20.d0   ! Thermal Convective const of the FreeFlow boundary layer
   ! Net radiation constants
   real(c_double) :: atm_flux_radiation = 0.d0   ! 340.d0 W/m^2
   real(c_double) :: soil_emissivity = 0.d0   ! 0.97.d0
   real(c_double) :: Stephan_Boltzmann_cst = 5.67d-8  ! W/m^2/K^4
   ! Rain input flux
   real(c_double) :: rain_flux(2) = RESHAPE((/0.d0, &  ! gas rain = 0
                                              -0.d0 &  ! liquid rain source term
                                              /), (/2/))

   !FIXME: must be put elsewhere and must be an array
   real(c_double) :: Thickness = 1.d0 !< Thickness of the fractures

   !FIXME: must be put elsewhere and must be an array
   !       we wait for v.4.1.0 and parameter distribution from python
   ! Volumetric heat capacity
   real(c_double) :: CpRoche = 800.d0*2000.d0 !< J/m3

   type, bind(C) :: Xalpha
      real(c_double) :: pressure, temperature, molar_fractions(NbComp)
   end type Xalpha

contains

   !< P is phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   function make_Xalpha(p, T, C) result(X)
      real(c_double), intent(in) :: p, T, C(NbComp)
      type(Xalpha) :: X
      X%pressure = p
      X%temperature = T
      X%molar_fractions = C
   end function make_Xalpha

   !< P is phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine extract_from_Xalpha(X, p, T, C)
      type(Xalpha), intent(in) :: X
      real(c_double), intent(out) :: p, T, C(NbComp)
      p = X%pressure
      T = X%temperature
      C = X%molar_fractions
   end subroutine extract_from_Xalpha

end module Physics
