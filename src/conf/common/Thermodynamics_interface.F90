!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Thermodynamics_interface

   use, intrinsic :: iso_c_binding, only: &
      c_int, c_double, c_ptr, c_loc
   use DefModel, only: NbPhase, NbComp
   use Physics, only: Xalpha, make_Xalpha, extract_from_Xalpha

   implicit none

   abstract interface
      subroutine law_with_derivatives(X, f, dfdX)
         import Xalpha, c_ptr
         type(Xalpha), intent(in) :: X, dfdX
         type(c_ptr), value, intent(in) :: f
      end subroutine law_with_derivatives
   end interface

   type Law_with_derivatives_ptr
      procedure(law_with_derivatives), pointer, nopass :: f
   end type Law_with_derivatives_ptr

   type(Law_with_derivatives_ptr), dimension(NbPhase) :: viscosity_with_derivatives

   abstract interface
      function law(X) result(f)
         import Xalpha, c_double
         type(Xalpha), intent(in) :: X
         real(c_double) :: f
      end function law
   end interface

   type Law_ptr
      procedure(law), pointer, nopass :: f
   end type Law_ptr

   type(Law_ptr), dimension(NbPhase) :: viscosity

contains

   ! Phase viscosity
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   function f_Viscosity(iph, p, T, C) result(f) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double) :: f

      ! use viscosity defined in Python
      f = viscosity(iph)%f(make_Xalpha(p, T, C))

   end function f_Viscosity

   ! Phase viscosity and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_Viscosity_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f
      type(Xalpha) :: dfdX

      ! use viscosity defined in Python
      call viscosity_with_derivatives(iph)%f( &
         make_Xalpha(p, T, C), c_loc(f), dfdX)
      call extract_from_Xalpha(dfdX, dfdP, dfdT, dfdC)

   end subroutine f_Viscosity_with_derivatives

end module Thermodynamics_interface