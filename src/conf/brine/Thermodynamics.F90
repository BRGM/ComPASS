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
      ! FIXME #51 f_VolumetricMassDensity = f_MolarDensity
      f_MolarDensity_with_derivatives, f_MolarDensity, &
      f_VolumetricMassDensity_with_derivatives, f_VolumetricMassDensity, &
      f_Viscosity_with_derivatives, f_Viscosity, &
      f_MolarEnthalpy_with_derivatives, f_MolarEnthalpy

#ifndef NDEBUG
   use CommonMPI, only: CommonMPI_abort
#endif

   implicit none

   public :: &
      f_Fugacity, & ! Fugacity (raises error because single phase)
      f_EnergieInterne

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

   ! EnergieInterne
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_EnergieInterne(iph, P, T, C, f, dPf, dTf, dCf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_MolarEnthalpy_with_derivatives(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_EnergieInterne

end module Thermodynamics
