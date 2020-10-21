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
   use IncCVReservoir, only: TYPE_IncCVReservoir

   implicit none

   public :: &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat, &
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy

contains

   ! Fugacity coefficient
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
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

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0

      if (iph == GAS_PHASE) then
         f = pa(GAS_PHASE)
         dPf = 1.d0
         dSf = dpadS(GAS_PHASE)
      else if (iph == LIQUID_PHASE) then
         call FluidThermodynamics_Psat(inc%Temperature, f, dTf)
#ifndef NDEBUG
      else
         call CommonMPI_abort('Unknow phase in f_Fugacity')
#endif
      end if

   end subroutine f_Fugacity

   pure subroutine f_DensiteMolaire_gas(p, T, C, f, dPf, dTf, dCf)
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: Z, dZ
      real(c_double), parameter :: u = 0.018016d0
      real(c_double), parameter :: R = 8.3145d0

      Z = 1.d0 ! CHECKME: ?
      dZ = 0.d0

      f = p*u/(R*T*Z)
      dPf = u/(R*T*Z)
      dTf = -p*u/(R*T*Z)**2*(R*T*dZ + R*Z)
      dCf = 0.d0

   end subroutine f_DensiteMolaire_gas

   pure subroutine f_DensiteMolaire_liquid(p, T, C, f, dPf, dTf, dCf)
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double), parameter :: rho0 = 780.83795d0
      real(c_double), parameter :: a = 1.6269192d0
      real(c_double), parameter :: b = -3.0635410d-3
      real(c_double), parameter :: a1 = 2.4638d-9
      real(c_double), parameter :: a2 = 1.1343d-17
      real(c_double), parameter :: b1 = -1.2171d-11
      real(c_double), parameter :: b2 = 4.8695d-20
      real(c_double), parameter :: c1 = 1.8452d-14
      real(c_double), parameter :: c2 = -5.9978d-23
      real(c_double) :: Cs, cw, dcwp, dcwt
      real(c_double) :: ds, ss, rs
      real(c_double) :: Psat, dT_Psat, prel

      Cs = 0.d0 ! salinity
      rs = 0.d0

      call FluidThermodynamics_Psat(T, Psat, dT_Psat)

      ss = (rho0 + a*T + b*T**2)*(1.d0 + 6.51d-4*Cs)
      ds = (a + b*T*2.d0)*(1.d0 + 6.51d-4*Cs)
      cw = (1.d0 + 5.d-2*rs) &
           *(a1 + a2*p + T*(b1 + b2*p) + T**2*(c1 + c2*p))
      dcwp = (1.d0 + 5.d-2*rs)*(a2 + T*b2 + T**2*c2)
      dcwt = (1.d0 + 5.d-2*rs)*((b1 + b2*p) + T*2.d0*(c1 + c2*p))
      prel = p - Psat
      f = ss*(1.d0 + cw*prel)
      dPf = ss*dcwp*prel + ss*cw
      dTf = ds*(1.d0 + cw*prel) + ss*dcwt*prel - ss*cw*dT_Psat
      dCf = 0.d0

   end subroutine f_DensiteMolaire_liquid

#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMolaire(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), value, intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      if (iph == GAS_PHASE) then
         call f_DensiteMolaire_gas(p, T, C, f, dPf, dTf, dCf)
      else if (iph == LIQUID_PHASE) then
         call f_DensiteMolaire_liquid(p, T, C, f, dPf, dTf, dCf)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_DensiteMolaire')
#endif
      endif

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
      subroutine f_DensiteMassique(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), value, intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
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
      subroutine f_Viscosite(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: Tref, a, da, b

      if (iph == GAS_PHASE) then
         f = (0.361d0*T - 10.2d0)*1.d-7
         dPf = 0.d0
         dTf = 0.361*1.d-7
         dCf = 0.d0
      else if (iph == LIQUID_PHASE) then
         Tref = T - 273.d0 - 8.435d0
         b = sqrt(8078.4d0 + Tref**2)
         a = 0.021482d0*(Tref + b) - 1.2d0
         da = 0.021482d0*(1.d0 + Tref/b)
         f = 1.d-3/a
         dPf = 0.d0
         dTf = -f*(da/a)
         dCf = 0.d0
#ifndef NDEBUG
      else
         call CommonMPI_abort('Unknow phase in f_Viscosite')
#endif

      end if

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
      subroutine f_EnergieInterne(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), value, intent(in) :: iph
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
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   ! If Enthalpide depends on the compositon C, change DefFlash.F90
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Enthalpie(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")
      integer(c_int), value, intent(in) :: iph
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

   ! Specific Enthalpy (used in FreeFlow)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_SpecificEnthalpy(iph, p, T, f, dPf, dTf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")
      integer(c_int), value, intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

#ifndef NDEBUG
      if (NbComp /= 1) &
         call CommonMPI_abort("Specific enthalpy formulation is ok for a single component...")
#endif

      if (iph == GAS_PHASE) then
         call f_specific_enthalpy_gas(p, T, f(1), dPf(1), dTf(1))
      else if (iph == LIQUID_PHASE) then
         call f_specific_enthalpy_liquid(p, T, f(1), dPf(1), dTf(1))
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_Enthalpie')
#endif
      end if

   end subroutine f_SpecificEnthalpy

   !< T is the Temperature
   pure subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

      real(c_double), value, intent(in) :: T
      real(c_double), intent(out) :: Psat, dT_PSat

      Psat = (T - 273.d0)**4.d0/1.0d3
      dT_PSat = 4.d0*(T - 273.d0)**3.d0/1.0d3

   end subroutine FluidThermodynamics_Psat

   !< P is the Reference Pressure
   pure subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

      real(c_double), value, intent(in) :: P
      real(c_double), intent(out) :: Tsat, dP_Tsat

      Tsat = 100.d0*(P/1.d5)**0.25d0 + 273.d0
      dP_Tsat = 0.25d0*1.d-3*(P/1.d5)**(-0.75d0)

   end subroutine FluidThermodynamics_Tsat

end module Thermodynamics
