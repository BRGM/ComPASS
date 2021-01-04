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

#ifndef NDEBUG
   use CommonMPI, only: CommonMPI_abort
#endif

   implicit none

   public :: &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy

contains

   ! Fugacity coefficients
   !< rt is the rocktype identifier
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< icp component identifier
   !< P is the reference pressure
   !< T is the temperature
   !< C is the phase molar frcations
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Fugacity(iph, icp, inc, pa, dpadS, f, DPf, DTf, DCf, DSf)
      integer(c_int), intent(in) :: iph, icp
      type(TYPE_IncCVReservoir), intent(in) :: inc
      real(c_double), intent(in) :: pa(NbPhase) ! p^\alpha: phase pressure
      real(c_double), intent(in) :: dpadS(NbPhase)
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

#ifndef NDEBUG
      call CommonMPI_abort("Should never be called with a single phase.")
#endif

   end subroutine f_Fugacity

   ! FIXME #51 densite massique
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMolaire(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: Cs, cw, dcwdp, dcwdT
      real(c_double) :: sc, dscdT, dscdCs
      real(c_double) :: Psat, dPsatdT
      real(c_double) :: rs, cwrs

      real(c_double), parameter :: rho0 = 780.83795d0
      real(c_double), parameter :: a = 1.6269192d0
      real(c_double), parameter :: b = -3.0635410d-3
      real(c_double), parameter :: a1 = 2.4638d-9
      real(c_double), parameter :: a2 = 1.1343d-17
      real(c_double), parameter :: b1 = -1.2171d-11
      real(c_double), parameter :: b2 = 4.8695d-20
      real(c_double), parameter :: c1 = 1.8452d-14
      real(c_double), parameter :: c2 = -5.9978d-23

      ! Pure water part
      rs = 0.d0 ! CHECKME: What is this???? Link to capillary pressure?
      cwrs = 1.d0 + 5.d-2*rs
      cw = cwrs*(a1 + a2*P + T*(b1 + b2*P) + T**2*(c1 + c2*P))
      dcwdp = cwrs*(a2 + b2*T + c2*T**2)
      dcwdT = cwrs*((b1 + b2*P) + T*2.d0*(c1 + c2*P))

      ! Salt correction
      Cs = C(ComPASS_SALT_COMPONENT)
      sc = (rho0 + a*T + b*T**2)*(1.d0 + 6.51d-4*Cs)
      dscdT = (a + 2.d0*b*T)*(1.d0 + 6.51d-4*Cs)
      dscdCs = 6.51d-4*(rho0 + a*T + b*T**2)

      !  pure water
      Psat = (T - 273.d0)**4.d0/1.0d3
      dPsatdT = 4.d0*(T - 273.d0)**3.d0/1.0d3

      f = sc*(1.d0 + cw*(P - Psat))
      dPf = sc*dcwdp*(P - Psat) + sc*cw
      dTf = dscdT*(1.d0 + cw*(P - Psat)) + sc*(dcwdT*(P - Psat) - cw*dPsatdT)
      dCf(ComPASS_SALT_COMPONENT) = dscdCs*(1.d0 + cw*(P - Psat))
      dCf(ComPASS_WATER_COMPONENT) = -dCf(ComPASS_SALT_COMPONENT)

   end subroutine f_DensiteMolaire

   ! FIXME #51 densite massique
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMassique(iph, P, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_DensiteMolaire(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_DensiteMassique

   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Viscosite(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: Tref, ns, dnsdT, b, dbdT, Cs, sc, dscdCs

      Tref = T - 273.d0 - 8.435d0
      b = sqrt(8078.4d0 + Tref**2)
      dbdT = Tref/b
      ! ns : no salt
      ns = 0.021482d0*(Tref + b) - 1.2
      dnsdT = 0.021482d0*(1.d0 + dbdT)
      ! sc: salt correction
      Cs = C(ComPASS_SALT_COMPONENT)
      sc = 1.d0 + Cs*1.34d0 + 6.12d0*Cs**2
      dscdCs = 1.34d0 + 2*6.12d0*Cs

      f = 1.d-3*sc/ns
      dPf = 0.d0
      dTf = -1.d-3*sc*(dnsdT/(ns**2))
      dCf(ComPASS_SALT_COMPONENT) = 1.d-3*(dscdCs)/ns
      dCf(ComPASS_WATER_COMPONENT) = -dCf(ComPASS_SALT_COMPONENT)

   end subroutine f_Viscosite

   ! EnergieInterne
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
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
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
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

   ! Specific Enthalpy (used in FreeFlow)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_SpecificEnthalpy(iph, P, T, f, dPf, dTf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

      real(c_double) :: fv, dPfv, dTfv

      ! FIXME: this should depend on salt concentration
#ifndef NDEBUG
      call CommonMPI_abort("f_SpecificEnthalpy: not implemented correctly.")
#endif

      call f_proxy_enthalpy(P, T, fv, dPfv, dTfv)
      f = fv
      dPf = dPfv
      dTf = dTfv

   end subroutine f_SpecificEnthalpy

end module Thermodynamics
