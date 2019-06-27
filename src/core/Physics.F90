!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Physics
    
    use iso_c_binding, only: c_double

    real(c_double) :: gravity = 9.81d0

    ! Far field atm constants
    real(c_double) :: atm_pressure = 1.d5
    real(c_double) :: atm_temperature = 300.d0
    real(c_double) :: atm_comp(2, 2) = &
                                    RESHAPE((/ &
                                        0.999d0, 1.d-3, & ! gas
                                        0.d0,    1.d0   & ! liquid
                                    /), (/2, 2/))
    ! Convective-diffusive boundary layer constants
    real(c_double) :: Hm(2) = RESHAPE((/ 20.d0/29.d0, &  ! gas Convective-diffusive const of the FreeFlow boundary layer
                                         0.d0         &  ! everything which is not gas is nul
                                         /), (/2/))
    real(c_double) :: HT =  20.d0   ! Thermal Convective const of the FreeFlow boundary layer
    ! Net radiation constants
    real(c_double) :: atm_flux_radiation = 340.d0   ! W/m^2
    real(c_double) :: soil_emissivity = 0.97d0
    real(c_double) :: Stephan_Boltzmann_cst = 5.67d-8  ! W/m^2/K^4
    ! Rain input flux
    real(c_double) :: rain_flux(2) = RESHAPE((/ 0.d0,   &  ! gas rain = 0
                                                -32.d-3 &  ! liquid rain source term
                                                /), (/2/))

    !FIXME: must be put elsewhere and must be an array
    real(c_double) :: Thickness = 1.d0 !< Thickness of the fractures

    !FIXME: must be put elsewhere and must be an array
    !       we wait for v.4.1.0 and parameter distribution from python
    ! Volumetric heat capacity
    real(c_double) :: CpRoche = 800.d0*2000.d0 !< J/m3

end module Physics
