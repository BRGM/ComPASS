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

    real(c_double) :: atm_pressure = 1.d5
    
    !FIXME: must be put elsewhere and must be an array
    real(c_double) :: Thickness = 1.d0 !< Thickness of the fractures

    !FIXME: must be put elsewhere and must be an array
    !       we wait for v.4.1.0 and parameter distribution from python
    ! Volumetric heat capacity
    real(c_double) :: CpRoche = 800.d0*2000.d0 !< J/m3

end module Physics
