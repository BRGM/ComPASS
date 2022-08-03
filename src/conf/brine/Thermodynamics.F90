!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: 2 phase 1 comp thermal, MCP=(1,1)

! 1: Gas
! 2: Water

module Thermodynamics

   use, intrinsic :: iso_c_binding, only: c_double, c_int
   use DefModel, only: NbPhase, NbComp, IndThermique
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir
   use Thermodynamics_interface, only: &
      f_MolarDensity_with_derivatives, f_MolarDensity, &
      f_Viscosity_with_derivatives, f_Viscosity

#ifndef NDEBUG
   use CommonMPI, only: CommonMPI_abort
#endif

   implicit none

   public :: &
      f_Fugacity, & ! Fugacity (raises error because single phase)
      f_DensiteMassique, & ! \rho^alpha(P,T,C)
      f_EnergieInterne, &
      f_Enthalpie

contains

   ! Fugacity coefficients
   !< icp component identifier
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the phase pressure
   !< T is the temperature
   !< C is the phase molar frcations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Fugacity(icp, iph, p, T, C, f, dfdp, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_fugacity")
      integer(c_int), intent(in) :: icp, iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdp, dfdT, dfdC(NbComp)

#ifndef NDEBUG
      call CommonMPI_abort("Should never be called with a single phase.")
#endif

   end subroutine f_Fugacity

   ! FIXME #51 densite massique
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_DensiteMassique(iph, P, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_MolarDensity_with_derivatives(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_DensiteMassique

   ! EnergieInterne
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_EnergieInterne(iph, P, T, C, f, dPf, dTf, dCf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_Enthalpie(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_EnergieInterne

   ! FIXME: pure liquid water...
   pure subroutine f_proxy_enthalpy(P, T, f, dPf, dTf)
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(out) :: f, dPf, dTf

      real(c_double), parameter :: a = -14.4319d+3
      real(c_double), parameter :: b = 4.70915d+3
      real(c_double), parameter :: cc = -4.87534d0
      real(c_double), parameter :: d = 1.45008d-2
      real(c_double), parameter :: T0 = 273.d0

      f = a + b*(T - T0) + cc*(T - T0)**2 + d*(T - T0)**3
      dPf = 0.d0
      dTf = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2

   end subroutine f_proxy_enthalpy

   ! Enthalpie
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   ! If Enthalpide depends on the compositon C, change DefFlash.F90
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Enthalpie(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

! FIXME: this should depend on salt concentration and use specific enthalpy
#ifndef NDEBUG
      call CommonMPI_abort("f_Enthalpie: not implemented correctly.")
#endif

      call f_proxy_enthalpy(P, T, f, dPf, dTf)
      dCf = 0.d0

   end subroutine f_Enthalpie

end module Thermodynamics
