!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module PhysicalConstants

   use iso_c_binding, only: c_double

   implicit none

   ! Molar masses
   real(c_double), parameter :: M_H2O = 18.d-3
   real(c_double), parameter :: M_air = 29.d-3

end module PhysicalConstants
