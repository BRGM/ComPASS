!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module SchemeParameters

   use iso_c_binding

   implicit none

   ! ****** Times ****** ! !

   ! One day
   real(c_double), parameter :: OneSecond = 1.d0
   real(c_double), parameter :: OneDay = 24.d0*3600.d0
   real(c_double), parameter :: OneMonth = 3.d1*OneDay
   real(c_double), parameter :: OneYear = 3.6525d2*OneDay

   ! time step init and max
   ! FIXME: parameter is removed to assign variable from python
   real(c_double) :: TimeFinal = 30*OneYear

   real(c_double) :: TimeStepInit = OneDay
   real(c_double) :: TimeStepMax = OneYear

   ! output_frequency for visu
   real(c_double), parameter :: output_frequency = OneYear

   ! ! ****** Newton iters max and stop condition ****** ! !
   integer, parameter :: NewtonNiterMax = 40
   real(c_double), parameter :: NewtonTol = 1.d-5

   ! ! ****** ksp linear solver iters max and stop condition ****** ! !
   integer, parameter :: KspNiterMax = 150 ! max nb of iterations
   real(c_double), parameter :: KspTol = 1.d-6 ! tolerance

   ! ! ****** Obj values used to compute Newton increment ****** ! !

   real(c_double), parameter :: &
      NewtonIncreObj_P = 5.d5, &
      NewtonIncreObj_T = 20.d0, &
      NewtonIncreObj_C = 1.d0, &
      NewtonIncreObj_S = 0.2d0

   ! ! ****** Obj values used to compute next time step ****** ! !

   real(c_double), parameter :: &
      TimeStepObj_P = 5.d5, &
      TimeStepObj_T = 20.d0, &
      TimeStepObj_C = 1.d0, &
      TimeStepObj_S = 0.6d0

   ! ! ****** Parameters of VAG schme (volume distribution) ****** ! !

   real(c_double), parameter :: &
      omegaDarcyCell = 0.075, & ! darcy cell/frac
      omegaDarcyFrac = 0.15

   real(c_double), parameter :: &
      omegaFourierCell = 0.075, & ! fourier cell/frac
      omegaFourierFrac = 0.15

   ! ! ****** Others ****** ! !

   ! eps
   real(c_double), parameter :: eps = 1.0d-10

end module SchemeParameters
