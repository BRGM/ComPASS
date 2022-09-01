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
   use CommonMPI, only: CommonMPI_abort
   use DefModel, only: &
      NbPhase, NbComp, IndThermique, &
      GAS_PHASE, LIQUID_PHASE
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir
   use Thermodynamics_interface, only: &
      ! FIXME #51 f_VolumetricMassDensity = f_MolarDensity
      f_MolarDensity_with_derivatives, f_MolarDensity, &
      f_VolumetricMassDensity_with_derivatives, f_VolumetricMassDensity, &
      f_Viscosity_with_derivatives, f_Viscosity, &
      f_MolarEnthalpy_with_derivatives, f_MolarEnthalpy

   implicit none

   public :: &
      f_Fugacity, & ! Fugacity
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat, &
      f_EnergieInterne

contains

   ! Fugacity = fugacity coefficient *  component phase molar fraction
   !< icp component identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< P is the phase pressure
   !< T is the temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Fugacity(icp, iph, p, T, C, f, dfdp, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_fugacity")
      integer(c_int), intent(in) :: icp, iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdp, dfdT, dfdC(NbComp)

      dfdp = 0.d0
      dfdT = 0.d0
      dfdC = 0.d0

      if (iph == GAS_PHASE) then
         f = p*C(icp)
         dfdp = C(icp)
         dfdC(icp) = p
      else if (iph == LIQUID_PHASE) then
         ! FIXME: use Poynting correction?
         call FluidThermodynamics_Psat(T, f, dfdT)
         ! f is Psat here
         ! we compute fugacity = Psat * C(icp) in place
         ! so the order of operations below is important
         dfdT = dfdT*C(icp)
         dfdC(icp) = f
         f = f*C(icp)
#ifndef NDEBUG
      else
         call CommonMPI_abort('Unknow phase in f_Fugacity')
#endif
      end if

   end subroutine f_Fugacity

   ! EnergieInterne
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_EnergieInterne(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_MolarEnthalpy_with_derivatives(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_EnergieInterne

   !< T is the Temperature
   pure subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")
      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: Psat, dT_PSat

      Psat = (T - 273.d0)**4.d0/1.0d3
      dT_PSat = 4.d0*(T - 273.d0)**3.d0/1.0d3

   end subroutine FluidThermodynamics_Psat

   !< P is a Pressure
#ifdef NDEBUG
   pure &
#endif
      subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")
      real(c_double), intent(in) :: P
      real(c_double), intent(out) :: Tsat, dP_Tsat

#ifndef NDEBUG
      if (P < 0) then
         write (*, *) "Saturation temperature at pressure:", P
         call CommonMPI_abort("Negative reference pressure in FluidThermodynamics_Tsat!")
      end if
#endif

      Tsat = 100.d0*(P/1.d5)**0.25d0 + 273.d0
      dP_Tsat = 0.25d0*1.d-3*(P/1.d5)**(-0.75d0)

   end subroutine FluidThermodynamics_Tsat

end module Thermodynamics
