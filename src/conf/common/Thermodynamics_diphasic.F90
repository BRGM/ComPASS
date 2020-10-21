!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Thermodynamics

   use iso_c_binding, only: c_double, c_int
   use CommonMPI, only: CommonMPI_abort
   use PhysicalConstants, only: M_H2O, M_air
   use DefModel, only: NbPhase, NbComp, IndThermique, &
                       GAS_PHASE, LIQUID_PHASE, WATER_COMP, AIR_COMP, &
                       get_model_configuration, MCP
   use CommonType, only: ModelConfiguration
   use IncCVReservoir, only: TYPE_IncCVReservoir

   implicit none

   public :: &
      f_Fugacity, &  !< Fucacity
      f_DensiteMolaire, &  !< \xi^alpha(P,T,C,S)
      f_DensiteMassique, &  !< \rho^alpha(P,T,C,S)
      f_Viscosite, &  !< \mu^alpha(P,T,C,S)
      air_henry, &  !< Henry coef for air comp
      air_henry_dT  !< derivative of the Henry coef for air comp

#ifdef _THERMIQUE_
   public :: &
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy, &
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat
#endif

contains

   ! Fugacity coefficient
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< icp component identifier
   !< P is the reference pressure
   !< T is the temperature
   !< C is the phase molar frcations
   !< S is all the saturations
   pure subroutine f_Fugacity(iph, icp, inc, pa, dpadS, f, DPf, DTf, DCf, DSf)
      integer(c_int), intent(in) :: iph, icp
      type(TYPE_IncCVReservoir), intent(in) :: inc
      real(c_double), intent(in) :: pa(NbPhase) ! p^\alpha: phase pressure
      real(c_double), intent(in) :: dpadS(NbPhase)
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      real(c_double) :: T, PSat, dTSat, Pc, dPcdS(NbPhase)
      real(c_double), parameter :: RZetal = 8.314d0*1000.d0/0.018d0

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0

      if (iph == GAS_PHASE) then
         f = pa(GAS_PHASE)
         dPf = 1.d0
         dSf(GAS_PHASE) = dpadS(GAS_PHASE)
      else if (iph == LIQUID_PHASE) then
         T = inc%Temperature
         if (icp == AIR_COMP) then
            call air_henry(T, f)
            call air_henry_dT(dTf)
         else if (icp == WATER_COMP) then
            Pc = pa(GAS_PHASE) - pa(LIQUID_PHASE)
            dPcdS = 0.d0
            dPcdS(GAS_PHASE) = dpadS(GAS_PHASE)
            dPcdS(LIQUID_PHASE) = -dpadS(LIQUID_PHASE)
            call FluidThermodynamics_Psat(T, Psat, dTSat)
            f = Psat*dexp(Pc/(T*RZetal))
            dTf = (dTSat - Psat*Pc/RZetal/(T**2))*dexp(Pc/(T*RZetal))
            dCf = 0.d0
            dSf = dPcdS*f/(T*RZetal)
         endif
      endif
   end subroutine f_Fugacity

   !> \brief Henry coef for air comp
   pure subroutine air_henry(T, H)

      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: H

      real(c_double), parameter :: T1 = 293.d0
      real(c_double), parameter :: T2 = 353.d0
      real(c_double), parameter :: H1 = 6.d+9
      real(c_double), parameter :: H2 = 10.d+9

      H = H1 + (H2 - H1)*(T - T1)/(T2 - T1)

   end subroutine

   !> \brief Derivative of the Henry coef for air comp wrt Temperature
   pure subroutine air_henry_dT(H_dt)

      real(c_double), intent(out) :: H_dt

      real(c_double), parameter :: T1 = 293.d0
      real(c_double), parameter :: T2 = 353.d0
      real(c_double), parameter :: H1 = 6.d+9
      real(c_double), parameter :: H2 = 10.d+9

      H_dt = (H2 - H1)/(T2 - T1)

   end subroutine

   pure subroutine f_DensiteMolaire_gas(p, T, C, f, dPf, dTf, dCf)
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: Pc, DSPc(NbPhase), Pg
      real(c_double), parameter :: Rgp = 8.314d0

      f = p/(Rgp*T)
      dPf = 1/(Rgp*T)
      dTf = -Pg/Rgp/T**2
      dCf = 0.d0

   end subroutine f_DensiteMolaire_gas

   pure subroutine f_DensiteMolaire_liquid(p, T, C, f, dPf, dTf, dCf)
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      f = 1000.d0/M_H2O
      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0

   end subroutine f_DensiteMolaire_liquid

   ! Molar density
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMolaire(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      if (iph == GAS_PHASE) then
         call f_DensiteMolaire_gas( &
            p, T, C, f, dPf, dTf, dCf)
      else if (iph == LIQUID_PHASE) then
         call f_DensiteMolaire_liquid( &
            p, T, C, f, dPf, dTf, dCf)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_DensiteMolaire')
#endif
      endif

   end subroutine f_DensiteMolaire

   ! Massic density
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMassique(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: zeta, M
      integer :: icp

      call f_DensiteMolaire(iph, p, T, C, zeta, dPf, dTf, dCf)

      M = MCP(AIR_COMP, iph)*M_air*C(AIR_COMP) + MCP(WATER_COMP, iph)*M_H2O*C(WATER_COMP)
      f = zeta*M
      dPf = dPf*M
      dTf = dTf*M
      dCf(AIR_COMP) = dCf(AIR_COMP)*M + zeta*MCP(AIR_COMP, iph)*M_air
      dCf(WATER_COMP) = dCf(WATER_COMP)*M + zeta*MCP(WATER_COMP, iph)*M_H2O

   end subroutine f_DensiteMassique

   ! Viscosities
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   pure subroutine f_Viscosite(iph, p, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      if (iph == GAS_PHASE) then
         f = 2.d-5
      else if (iph == LIQUID_PHASE) then
         f = 1.d-3
      endif

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0

   end subroutine f_Viscosite

#ifdef _THERMIQUE_

   ! Internal energy = enthalpie - Pressure (here everything is volumic)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations

#ifdef NDEBUG
   pure &
#endif
      subroutine f_EnergieInterne(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: paovzeta2
      real(c_double) :: zeta, dzetadP, dzetadT, dzetadC(NbComp)
      real(c_double) :: enth, denthdP, denthdT, denthdC(NbComp)

      call f_Enthalpie(iph, p, T, C, enth, denthdP, denthdT, denthdC)
      call f_DensiteMolaire(iph, p, T, C, zeta, dzetadP, dzetadT, dzetadC)
      f = enth - p/zeta
      paovzeta2 = p/(zeta**2)
      dPf = denthdP - 1.d0/zeta + paovzeta2*dzetadP
      dTf = denthdT + paovzeta2*dzetadT
      dCf = denthdC + paovzeta2*dzetadC

   end subroutine f_EnergieInterne

   pure subroutine f_specific_enthalpy_gas(p, T, is_present, f, dPf, dTf)
      real(c_double), intent(in) :: p, T
      integer, intent(in) :: is_present(NbComp)
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

      real(c_double), parameter :: a = 1990.89d+3
      real(c_double), parameter :: b = 190.16d+3
      real(c_double), parameter :: cc = -1.91264d+3
      real(c_double), parameter :: d = 0.2997d+3
      real(c_double) :: Ts, ss, dTss, cp, beta_air, beta_water

      call f_CpGaz(cp)

      Ts = T/100.d0
      ss = a + b*Ts + cc*Ts**2.d0 + d*Ts**3.d0
      dTss = (b + 2.d0*cc*Ts + 3.d0*d*Ts**2.d0)/100.d0
      ! CHECKME: could we assume that MCP(AIR_COMP, iph)==.false. => C(AIR_COMP) = 0
      beta_air = is_present(AIR_COMP)*cp*M_air
      beta_water = is_present(WATER_COMP)*M_H2O
      f(AIR_COMP) = beta_air*T
      f(WATER_COMP) = beta_water*ss
      dPf = 0.d0
      dTf(AIR_COMP) = beta_air
      dTf(WATER_COMP) = beta_water*dTss

   end subroutine f_specific_enthalpy_gas

   pure subroutine f_specific_enthalpy_liquid(p, T, f, dPf, dTf)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

      real(c_double), parameter :: a = -14.4319d+3
      real(c_double), parameter :: b = +4.70915d+3
      real(c_double), parameter :: cc = -4.87534
      real(c_double), parameter :: d = 1.45008d-2
      real(c_double), parameter :: T0 = 273.d0
      real(c_double) :: ss, dTss, TdegC

      TdegC = T - T0
      ss = a + b*TdegC + cc*(TdegC**2) + d*(TdegC**3)
      dTss = b + 2.d0*cc*TdegC + 3*d*(TdegC**2)
      ! FIXME: all components have the same contributions
      f = ss*M_H2O
      dPf = 0.d0
      dTf = dTss*M_H2O

   end subroutine f_specific_enthalpy_liquid

   ! Enthalpie
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_Enthalpie(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: fspec(NbComp), dfspecdP(NbComp), dfspecdT(NbComp)

      if (iph == GAS_PHASE) then
         call f_specific_enthalpy_gas(p, T, MCP(:, GAS_PHASE), fspec, dfspecdP, dfspecdT)
      else if (iph == LIQUID_PHASE) then
         call f_specific_enthalpy_liquid(p, T, fspec, dfspecdP, dfspecdT)
#ifndef NDEBUG
      else
         call CommonMPI_abort("Unknow phase in f_Enthalpie.")
#endif
      endif

      f = fspec(AIR_COMP)*C(AIR_COMP) + fspec(WATER_COMP)*C(WATER_COMP)
      dPf = dfspecdP(AIR_COMP)*C(AIR_COMP) + dfspecdP(WATER_COMP)*C(WATER_COMP)
      dTf = dfspecdT(AIR_COMP)*C(AIR_COMP) + dfspecdT(WATER_COMP)*C(WATER_COMP)
      dCf(AIR_COMP) = fspec(AIR_COMP)
      dCf(WATER_COMP) = fspec(WATER_COMP)

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
      real(c_double), value, intent(in) :: p, T
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

      if (iph == GAS_PHASE) then
         call f_specific_enthalpy_gas(p, T, MCP(:, GAS_PHASE), f, dPf, dTf)
      else if (iph == LIQUID_PHASE) then
         call f_specific_enthalpy_liquid(p, T, f, dPf, dTf)
#ifndef NDEBUG
      else
         call CommonMPI_abort("Unknow phase in f_SpecificEnthalpy.")
#endif
      endif

   end subroutine f_SpecificEnthalpy

   !> \brief Rock volumetric heat capacity
   pure subroutine f_CpGaz(c)

      real(c_double), intent(out) :: c

      c = 1000.d0

   end subroutine

   ! Compute Psat(T)
   !< T is the Temperature
   pure subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

      real(c_double), value, intent(in) :: T
      real(c_double), intent(out) :: Psat, dT_PSat

      ! valid between -50C and 200C
      Psat = 100.d0*dexp(46.784d0 - 6435.d0/T - 3.868d0*dlog(T))
      dT_PSat = (6435.d0/T**2.d0 - 3.868d0/T)*Psat

   end subroutine FluidThermodynamics_Psat

   ! Compute Tsat(P)
   !< P is the Reference Pressure
#ifdef NDEBUG
   pure &
#endif
      subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

      real(c_double), value, intent(in) :: P
      real(c_double), intent(out) :: Tsat, dP_Tsat

      Tsat = 0.d0
      dP_Tsat = 0.d0

#ifndef NDEBUG
      call CommonMPI_abort('entered in Tsat, not implemented')
#endif

   end subroutine FluidThermodynamics_Tsat

#endif

end module Thermodynamics
