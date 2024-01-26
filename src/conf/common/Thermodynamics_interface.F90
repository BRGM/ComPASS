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
   type(Law_with_derivatives_ptr), dimension(NbPhase) :: molar_density_with_derivatives
   type(Law_with_derivatives_ptr), dimension(NbPhase) :: volumetric_mass_density_with_derivatives
   type(Law_with_derivatives_ptr), dimension(NbPhase) :: molar_enthalpy_with_derivatives

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
   type(Law_ptr), dimension(NbPhase) :: molar_density
   type(Law_ptr), dimension(NbPhase) :: volumetric_mass_density
   type(Law_ptr), dimension(NbPhase) :: molar_enthalpy

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

   ! Phase molar density
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   function f_MolarDensity(iph, p, T, C) result(f) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double) :: f

      ! use molar density defined in Python
      f = molar_density(iph)%f(make_Xalpha(p, T, C))

   end function f_MolarDensity

   ! Phase molar density and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_MolarDensity_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_molar_density_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f
      type(Xalpha) :: dfdX

      ! use molar density defined in Python
      call molar_density_with_derivatives(iph)%f( &
         make_Xalpha(p, T, C), c_loc(f), dfdX)
      call extract_from_Xalpha(dfdX, dfdP, dfdT, dfdC)

   end subroutine f_MolarDensity_with_derivatives

   ! Phase volumetric mass density
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   function f_VolumetricMassDensity(iph, p, T, C) result(f) &
      bind(C, name="FluidThermodynamics_volumetric_mass_density")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double) :: f

      ! use volumetric mass density defined in Python
      f = volumetric_mass_density(iph)%f(make_Xalpha(p, T, C))

   end function f_VolumetricMassDensity

   ! Phase volumetric mass density and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_VolumetricMassDensity_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_volumetric_mass_density_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f
      type(Xalpha) :: dfdX

      ! use volumetric mass density defined in Python
      call volumetric_mass_density_with_derivatives(iph)%f( &
         make_Xalpha(p, T, C), c_loc(f), dfdX)
      call extract_from_Xalpha(dfdX, dfdP, dfdT, dfdC)

   end subroutine f_VolumetricMassDensity_with_derivatives

   ! Phase molar enthalpy
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   function f_MolarEnthalpy(iph, p, T, C) result(f) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double) :: f

      ! use molar enthalpy defined in Python
      f = molar_enthalpy(iph)%f(make_Xalpha(p, T, C))

   end function f_MolarEnthalpy

   ! Phase molar enthalpy and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_MolarEnthalpy_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_molar_enthalpy_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f
      type(Xalpha) :: dfdX

      ! use molar enthalpy defined in Python
      call molar_enthalpy_with_derivatives(iph)%f( &
         make_Xalpha(p, T, C), c_loc(f), dfdX)
      call extract_from_Xalpha(dfdX, dfdP, dfdT, dfdC)

   end subroutine f_MolarEnthalpy_with_derivatives

end module Thermodynamics_interface
