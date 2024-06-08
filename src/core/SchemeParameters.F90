!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module SchemeParameters

   use iso_c_binding, only: c_double

   implicit none

   ! defined in python
   ! ! ! ****** Obj values used to compute Newton increment ****** ! !

   ! real(c_double), parameter :: &
   !    NewtonIncreObj_P = 5.001d6, &
   !    NewtonIncreObj_T = 20.d0, &
   !    NewtonIncreObj_C = 1.d0, &
   !    NewtonIncreObj_S = 0.2d0

   ! eps
   real(c_double), parameter :: eps = 1.0d-10

end module SchemeParameters
