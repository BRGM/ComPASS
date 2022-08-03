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
      f_MolarDensity_with_derivatives, f_MolarDensity, &
      f_Viscosity_with_derivatives, f_Viscosity

   implicit none

   public :: &
      f_Fugacity, & ! Fugacity
      f_DensiteMassique, & ! \rho^alpha(P,T,C)
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat, &
      f_EnergieInterne, &
      f_Enthalpie

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

   ! FIXME #51 densite massique
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   subroutine f_DensiteMassique(iph, p, T, C, f, dfdP, dfdT, dfdC)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdP, dfdT, dfdC(NbComp)

      call f_MolarDensity_with_derivatives(iph, P, T, C, f, dfdP, dfdT, dfdC)

   end subroutine f_DensiteMassique

   ! EnergieInterne
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_EnergieInterne(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_Enthalpie(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_EnergieInterne

   pure subroutine f_specific_enthalpy_gas(p, T, f, dPf, dTf)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f, dPf, dTf

      real(c_double), parameter :: a = 1990.89d+3
      real(c_double), parameter :: b = 190.16d+1

      f = a + T*b
      dPf = 0.d0
      dTf = b

   end subroutine f_specific_enthalpy_gas

   pure subroutine f_specific_enthalpy_liquid(p, T, f, dPf, dTf)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f, dPf, dTf

      real(c_double), parameter :: a = -14.4319d+3
      real(c_double), parameter :: b = 4.70915d+3
      real(c_double), parameter :: cc = -4.87534d0
      real(c_double), parameter :: d = 1.45008d-2
      real(c_double), parameter :: T0 = 273.d0
      real(c_double) :: TdegC

      TdegC = T - T0
      f = a + b*TdegC + cc*(TdegC**2) + d*(TdegC**3)
      dPf = 0.d0
      dTf = b + 2.d0*cc*TdegC + 3.d0*d*(TdegC**2)

   end subroutine f_specific_enthalpy_liquid

   ! Enthalpie
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   ! If Enthalpide depends on the compositon C, change DefFlash.F90
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Enthalpie(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

#ifndef NDEBUG
      if (NbComp /= 1) &
         call CommonMPI_abort("Specific enthalpy formulation is ok for a single component...")
#endif

      if (iph == GAS_PHASE) then
         call f_specific_enthalpy_gas(p, T, f, dPf, dTf)
      else if (iph == LIQUID_PHASE) then
         call f_specific_enthalpy_liquid(p, T, f, dPf, dTf)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_Enthalpie')
#endif
      end if

      dCf = 0.d0

   end subroutine f_Enthalpie

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
