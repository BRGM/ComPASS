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
   use PhysicalConstants, only: M_H2O
   use DefModel, only: NbPhase, NbComp, IndThermique, &
                       GAS_PHASE, LIQUID_PHASE, WATER_COMP, CO2_COMP, &
                       get_model_configuration, MCP
   use CommonType, only: ModelConfiguration
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir

   implicit none

   public :: &
      f_Fugacity, &  !< Fugacity
      f_Fugacity_coeff_CO2_gas, & ! phi : correction of the ???
      f_MolarDensity, &  !< \zeta^alpha(P,T,C,S)
      f_MolarDensity_with_derivatives, &  !< \zeta^alpha(P,T,C,S)
      f_VolumetricMassDensity, &  !< \xi^alpha(P,T,C,S)
      f_VolumetricMassDensity_with_derivatives, &  !< \xi^alpha(P,T,C,S)
      f_Viscosity, &  !< \mu^alpha(P,T,C,S)
      f_Viscosity_with_derivatives  !< \mu^alpha(P,T,C,S)
   ! CO2_henry, &  !< Henry coef for CO2 comp
   ! f_PermRel, &  ! k_{r_alpha}(S)

#ifdef _THERMIQUE_
   public :: &
      f_EnergieInterne, &
      f_MolarEnthalpy, &
      f_MolarEnthalpy_with_derivatives, &
      ! f_SpecificEnthalpy, &
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat
#endif

contains

   pure subroutine equilibriumConstantH2O(T, f, dfdT)
      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: f, dfdT

      ! Constants for equilibrium constant calculation
      real(c_double) :: TinC, logk0_H2O
      real(c_double), parameter :: c(4) = [-2.209d0, 3.097d-2, -1.098d-4, 2.048d-7]

      TinC = T - 273.15d0 ! temperature in °C
      logk0_H2O = c(1) + c(2)*TinC + c(3)*TinC**2 + c(4)*TinC**3

      f = 10**logk0_H2O
      ! ddx(10**u) = ln(10)*ddx(u)*10**u
      dfdT = log(10.d0)*f*(c(2) + 2.d0*c(3)*TinC + 3.d0*c(4)*TinC**2)
   end subroutine equilibriumConstantH2O

   pure subroutine equilibriumConstantCO2(T, f, dfdT)
      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: f, dfdT

      ! Constants for equilibrium constant calculation
      real(c_double) :: TinC, logk0_CO2, dlogk0_CO2dT
      ! FIXE ME : add a condition for to switch in KCO2(l) 0 when the temperature
      ! is subcritical and the pressure is above the CO2 saturation pressure
      !c[3] = [ 1.189, 1.304e-2, -5.446e-5 ] ! KCO2(g) for the CO2g
      real(c_double), parameter :: c(3) = [1.169d0, 1.368d-2, -5.380d-5] ! KCO2(l) for the CO2l

      TinC = T - 273.15d0 ! temperature in °C
      logk0_CO2 = c(1) + c(2)*TinC + c(3)*TinC**2
      dlogk0_CO2dT = c(2) + 2.d0*c(3)*TinC

      f = 1.d5*10.d0**logk0_CO2 !converted from bar to Pa
      ! ddx(10**u) = ln(10)*ddx(u)*10**u
      dfdT = log(10.d0)*f*dlogk0_CO2dT
   end subroutine equilibriumConstantCO2

   pure subroutine f_linear_volume(p, T, f, dfdp, dfdT)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f, dfdp, dfdT

      real(c_double) :: MCO2

      call CO2_MasseMolaire(MCO2)

      ! -1.167d-5*(p/1.d6) + ...
      !f = (-1.167d-11*p + 6.539d-6*T - 0.00059d0)*MCO2 ! volume m^3/mol
      !dfdp = -1.167d-11*MCO2
      !dfdT = 6.539d-6*MCO2
      ! test with a mixed term P and T
      !f = (9.682709259446859d-11*p &
      !   + 1.643090334166191d-5*T &
      !   -3.297337244828735d-13*p*T&
      !   - 0.003850396226719975 &
      !)*MCO2
      !dfdp = 9.682709259446859d-11*MCO2 - 3.297337244828735d-13*T*MCO2
      !dfdT = 1.643090334166191d-5*MCO2 - 3.297337244828735d-13*p*MCO2
      ! test with a wider range of T
      f = (8.271046832872992d-11*p &
           + 1.4722010876444494d-5*T &
           - 2.8733291908818404d-13*p*T &
           - 0.00328128373967034 &
           )*MCO2
      dfdp = 8.271046832872992d-11*MCO2 - 2.8733291908818404d-13*T*MCO2
      dfdT = 1.4722010876444494d-5*MCO2 - 2.8733291908818404d-13*p*MCO2
   end subroutine f_linear_volume

   ! Fugacity is f_coeff_CO2_liquid_fug * Cl(CO2)
   !< p is the phase pressure
   !< T is the temperature
   !< does not depend on C (phase molar fractions)
   !< if f_coeff_CO2_liquid_fug depends on C, change the algo in Equilibrium.h
   subroutine f_coeff_CO2_liquid_fug(p, T, f, dfdp, dfdT) &
      bind(C, name="FluidThermodynamics_coeff_CO2_liquid_fug")
      real(c_double), intent(in) :: p, T
      real(c_double), intent(out) :: f, dfdp, dfdT

      real(c_double), parameter :: v_av_CO2 = 3.26d-5 ! [m^3/mol] average partial molar volume of CO2 over the pressure range
      real(c_double), parameter :: R = 8.314d0
      real(c_double) :: k0_CO2, dk0_CO2dT, deltaP, v_av_CO2RT, dv_av_CO2RTdT, c0T, dc0TdT, e_u

      ! IMPORTANT: at subcritical temperatures and pressures above saturation values,
      ! K0_CO2(T) needs to be replaced with another equilibrium constant, K0_CO2(l)
      call equilibriumConstantCO2(T, k0_CO2, dk0_CO2dT) ! equilibrium constant for CO2 at 1 bar
      deltaP = p - 1.d5 ! pressure difference [MPa] between reference pressure p0 = 1bar and pg
      v_av_CO2RT = v_av_CO2/(R*T)
      dv_av_CO2RTdT = -v_av_CO2/(R*T**2)
      c0T = 55.508d0*k0_CO2
      dc0TdT = 55.508d0*dk0_CO2dT
      e_u = exp(deltaP*v_av_CO2RT)
      f = c0T*e_u

      dfdp = c0T*v_av_CO2RT*e_u
      dfdT = e_u*(dc0TdT + c0T*deltaP*dv_av_CO2RTdT)

   end subroutine f_coeff_CO2_liquid_fug

   !> Fugacity coefficient of the CO2 under the gas phase
   ! Fug = Fugacity coefficient * p^g * C^g(CO2)
   pure subroutine f_Fugacity_coeff_CO2_gas(p, T, C, f, dfdp, dfdT, dfdC)

      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdp, dfdT, dfdC(NbComp)
      real(c_double), parameter :: b_CO2 = 2.78d-5 !converted from cm3/mol
      real(c_double), parameter :: R = 8.314d0
      real(c_double) :: a_CO2, da_CO2_dT, V, dVdp, dVdT, p_bar, c0T, c1T, sum_cT, ddTc0, ddTc1, &
                        ln_phi, ddp_ln_phi, ddT_ln_phi, phi

      call f_linear_volume(p, T, V, dVdp, dVdT)
      !p_bar = p/1.d5 ! gas phase pressure in bar
      !FIXME:TEST
      p_bar = p
      ! a_CO2 = (7.54d7 - 4.13d4*T)*1d-12 a_CO2 is originally in bar cm6 K0.5 mol–2
      !a_CO2 = 7.54d-5 - 4.13d-8*T ! mixture parameter of Redlich-Kwong equation [N. Spycher, K. Pruess, and J. Ennis-King 2003]
      !FIXME
      a_CO2 = 7.54 - 4.13d-3*T

      c0T = 2.d0*a_CO2/(R*T**1.5*b_CO2)
      ddTc0 = -2.d0*(4.13d-3*T**1.5 + 1.5d0*a_CO2*T**0.5)/(R*b_CO2*T**3)
      c1T = a_CO2/(R*T**1.5)
      ! ddTc1 = (-4.13d-3*T**1.5 - a_CO2*1.5d0*T**0.5)/(R*T**3)
      ddTc1 = (-4.13d-3/T**1.5 - a_CO2*1.5d0/T**2.5)/R
      sum_cT = c1T/b_CO2 - c0T

      ln_phi = log(V/(V - b_CO2)) &
               + b_CO2/(V - b_CO2) &
               + sum_cT*log((V + b_CO2)/V) &
               - c1T/(V + b_CO2) &
               - log(p_bar*V/(R*T))

      ! ddx(ln(u)) = ddx(u)/u
      ddp_ln_phi = dVdp/V - dVdp/(V - b_CO2) &
                   - b_CO2*dVdp/(V - b_CO2)**2 &
                   + sum_cT*(dVdp/(V + b_CO2) - dVdp/V) &
                   + c1T*dVdp/(V + b_CO2)**2 &
                   - 1.d0/p - dVdp/V

      ddT_ln_phi = dVdT/V - dVdT/(V - b_CO2) &
                   - b_CO2*dVdT/(V - b_CO2)**2 &
                   + (ddTc1/b_CO2 - ddTc0)*log((V + b_CO2)/V) + sum_cT*(dVdT/(V + b_CO2) - dVdT/V) &
                   - (ddTc1*(V + b_CO2) - c1T*dVdT)/(V + b_CO2)**2 &
                   - dVdT/V + 1.d0/T

      f = exp(ln_phi)
      dfdp = ddp_ln_phi*f
      dfdT = ddT_ln_phi*f
      dfdC(:) = 0.d0

   end subroutine f_Fugacity_coeff_CO2_gas

   !< icp component identifier (must be CO2)
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< p is the phase pressure
   !< T is the temperature
   !< C is the phase molar fractions
   subroutine f_Fugacity(icp, iph, p, T, C, f, dfdp, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_fugacity")
      integer(c_int), intent(in) :: icp, iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdp, dfdT, dfdC(NbComp)

      dfdp = 0.d0
      dfdT = 0.d0
      dfdC = 0.d0

#ifndef NDEBUG
      if (icp == WATER_COMP) then
         call CommonMPI_abort('f_Fugacity not implemented for water comp')
      elseif (icp .ne. CO2_COMP) then
         call CommonMPI_abort('f_Fugacity called with a wrong component')
      endif
#endif

      ! icp is CO2
      if (iph == GAS_PHASE) then
         ! Fug = Fugacity coefficient * p^g * C^g(icp)
         call f_Fugacity_coeff_CO2_gas(p, T, C, f, dfdp, dfdT, dfdC)
         ! we compute f = f * C(icp) * p^g in place
         ! so the order of operations below is important
         dfdp = C(icp)*(dfdp*p + f)
         dfdT = dfdT*C(icp)*p
         ! d (f(P,T,C)*C_i*p)/dC_i = p*(f + df/dC_i*C_i)
         ! d (f(P,T,C)*C_i*p)/dC_j = p*(    df/dC_j*C_i), j!=i
         dfdC = dfdC*C(icp) ! df/dC_j*C_i (for all j)
         dfdC(icp) = dfdC(icp) + f ! add f when j=i
         dfdC = dfdC*p ! p*(...) for all j
         f = f*C(icp)*p
      else if (iph == LIQUID_PHASE) then
         ! Fug = ctx * C^l(icp), ctx does not depend on C^l
         call f_coeff_CO2_liquid_fug(p, T, f, dfdp, dfdT)
         ! we compute f = f * C(icp) in place
         ! so the order of operations below is important
         dfdp = dfdp*C(icp)
         dfdT = dfdT*C(icp)
         dfdC = 0.d0
         ! d (f(P,T,C)*C_i)/dC_j
         dfdC(icp) = f ! derivative is f when j=i
         f = f*C(icp)
      endif

   end subroutine f_Fugacity

   pure subroutine f_DensiteMass_gas(p, T, C, f, dPf, dTf, dCf)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: u, usquare

      !u = -1.167d-11*p + 6.539d-6*T - 0.00059d0
      !u = 9.682709259446859d-11*p &
      !   + 1.643090334166191d-5*T &
      !  -3.297337244828735d-13*p*T&
      !   - 0.003850396226719975
      ! test with a wider range of T
      u = 8.271046832872992d-11*p &
          + 1.4722010876444494d-5*T &
          - 2.8733291908818404d-13*p*T &
          - 0.00328128373967034

      usquare = u**2
      f = 1/u
      !dPf = 1.167d-11/usquare
      !dTf = -6.539d-6/usquare
      !dPf = 9.682709259446859d-11 - 3.297337244828735d-13*T
      !dTf = 1.643090334166191d-5 - 3.297337244828735d-13*p
      dPf = -(8.271046832872992d-11 - 2.8733291908818404d-13*T)/usquare
      dTf = -(1.4722010876444494d-5 - 2.8733291908818404d-13*p)/usquare
      dCf(:) = 0.d0

   end subroutine f_DensiteMass_gas

   ! defined in SPE11 benchmark, taken from Garcia 2001
   pure subroutine f_DensiteMass_liquid(p, T, C, f, dPf, dTf, dCf)
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: V_theta, V_H2O, M_CO2, rho_CO2, ddTrho_CO2, rho_H20, ddTrho_H20, &
                        ddprho_H20, one_over_f, ddTone_over_f, ddpone_over_f, ddCO2one_over_f, ddH2Oone_over_f, &
                        one_over_f_square, TinC

      TinC = T - 273.15d0 ! temperature in °C
      call CO2_MasseMolaire(M_CO2)

      V_theta = 1.d-6*(37.51d0 - 9.585d-2*TinC + 8.740d-4*TinC**2 - 5.044d-7*TinC**3) !Volume in 1.d-6*(cm3/mole)
      V_H2O = -4.12258d-13*p + 4.8881d-7*T + 0.00085d0 !Volume (m3/kg)
      rho_CO2 = M_CO2/V_theta
      ddTrho_CO2 = M_CO2*1.d-6*(9.585d-2 - 2.d0*8.740d-4*TinC + 3.d0*5.044d-7*TinC**2)/V_theta**2
      rho_H20 = 1.d0/V_H2O
      ddTrho_H20 = -4.8881d-7/V_H2O**2
      ddprho_H20 = 4.12258d-13/V_H2O**2

      ! C(WATER_COMP) is 1 - C(CO2_COMP)
      one_over_f = C(WATER_COMP)/rho_H20 + C(CO2_COMP)/rho_CO2
      ddTone_over_f = -C(WATER_COMP)*ddTrho_H20/rho_H20**2 - C(CO2_COMP)*ddTrho_CO2/rho_CO2**2
      ddpone_over_f = -C(WATER_COMP)*ddprho_H20/rho_H20**2
      ddCO2one_over_f = 1.d0/rho_CO2
      ddH2Oone_over_f = 1.d0/rho_H20

      f = 1.d0/one_over_f
      one_over_f_square = one_over_f**2
      dTf = -ddTone_over_f/one_over_f_square
      dPf = -ddpone_over_f/one_over_f_square
      dCf(CO2_COMP) = -ddCO2one_over_f/one_over_f_square
      dCf(WATER_COMP) = -ddH2Oone_over_f/one_over_f_square

   end subroutine f_DensiteMass_liquid

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

      real(c_double) :: rho, M, M_CO2
      real(c_double) :: dfdP, dfdT, dfdC(NbComp)

      call f_VolumetricMassDensity_with_derivatives(iph, p, T, C, rho, dfdP, dfdT, dfdC)
      call CO2_MasseMolaire(M_CO2)
      if (iph == GAS_PHASE) then
         M = M_CO2
      else if (iph == LIQUID_PHASE) then
         M = M_CO2*C(CO2_COMP) + M_H2O*C(WATER_COMP)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_MolarDensity')
#endif
      endif
      f = rho/M

   end function f_MolarDensity

   ! Phase molar density and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_MolarDensity_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_molar_density_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f

      real(c_double) :: rho, M_CO2, M

      call f_VolumetricMassDensity_with_derivatives(iph, p, T, C, rho, dfdP, dfdT, dfdC)
      CALL CO2_MasseMolaire(M_CO2)
      if (iph == GAS_PHASE) then
         M = M_CO2
      else if (iph == LIQUID_PHASE) then
         M = M_CO2*C(CO2_COMP) + M_H2O*C(WATER_COMP)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_MolarDensity')
#endif
      endif
      f = rho/M
      dfdP = dfdP/M
      dfdT = dfdT/M
      dfdC(CO2_COMP) = dfdC(CO2_COMP)/M - f*M_CO2/M
      dfdC(WATER_COMP) = dfdC(WATER_COMP)/M - f*M_H2O/M
   end subroutine f_MolarDensity_with_derivatives

   pure subroutine CO2_MasseMolaire(m)

      DOUBLE PRECISION, INTENT(OUT) :: m

      m = 44.01d-3
   end subroutine CO2_MasseMolaire

   ! Volumetric mass density (no derivatives)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      function f_VolumetricMassDensity(iph, p, T, C) result(f) &
      bind(C, name="FluidThermodynamics_volumetric_mass_density")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double) :: f
      ! not necessary, just to use f_DensiteMolaire_gas and
      ! f_DensiteMolaire_liquid (no distinction without der)
      real(c_double) :: dfdP, dfdT, dfdC(NbComp)

      if (iph == GAS_PHASE) then
         call f_DensiteMass_gas( &
            p, T, C, f, dfdP, dfdT, dfdC)
      else if (iph == LIQUID_PHASE) then
         call f_DensiteMass_liquid( &
            p, T, C, f, dfdP, dfdT, dfdC)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_MolarDensity')
#endif
      endif
   end function f_VolumetricMassDensity

   ! Volumetric mass density (with derivatives)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_VolumetricMassDensity_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_volumetric_mass_density_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f

      if (iph == GAS_PHASE) then
         call f_DensiteMass_gas( &
            p, T, C, f, dfdP, dfdT, dfdC)
      else if (iph == LIQUID_PHASE) then
         call f_DensiteMass_liquid( &
            p, T, C, f, dfdP, dfdT, dfdC)
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_MolarDensity_with_derivatives')
#endif
      endif

   end subroutine f_VolumetricMassDensity_with_derivatives

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

      if (iph == GAS_PHASE) then
         !f = 15.d-6
         !f = 1.d-4 ! Value for CO2
         ! Linear model for CO2
         f = 8.150657587292315d-14*p - 9.483270717014074d-07*T &
             + 4.250176084439237d-15*p*T + 0.0003481359413684678d0
      else if (iph == LIQUID_PHASE) then
         f = 1.d-3
      endif

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

      if (iph == GAS_PHASE) then
         !f = 15.d-6
         !f = 1.d-4 ! Value for CO2
         ! Linear model for CO2
         f = 8.150657587292315d-14*p - 9.483270717014074d-07*T &
             + 4.250176084439237d-15*p*T + 0.0003481359413684678d0
         dfdP = 8.150657587292315d-14 + 4.250176084439237d-15*T
         dfdT = -9.483270717014074d-07 + 4.250176084439237d-15*p
         dfdC(:) = 0.d0
      else if (iph == LIQUID_PHASE) then
         f = 1.d-3
         dfdP = 0.d0
         dfdT = 0.d0
         dfdC = 0.d0
      endif

   end subroutine f_Viscosity_with_derivatives

#ifdef _THERMIQUE_

   ! Internal energy = enthalpie - Pressure (here everything is volumic)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   pure subroutine f_linear_cp(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      if (iph == GAS_PHASE) then
         f = -1.2896d-6*p + 0.1907d0*T + 67.0627d0
         dPf = -1.2896d-6
         dTf = 0.1907d0
         dCf(:) = 0.d0
      else if (iph == LIQUID_PHASE) then
         f = -0.0345d-6*P + 0.0105d0*T + 71.8364d0
         dPf = -0.0345d-6
         dTf = 0.0105d0
         dCf(:) = 0.d0
      end if
   end subroutine f_linear_cp

   ! Phase molar enthalpy (J/mol)
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

      !real(c_double) :: Cp, dCpdP, dCpdT, dCpdC(NbComp)

      !call f_linear_cp(iph, p, T, C, Cp, dCpdP, dCpdT, dCpdC)
      !f = Cp*T
      ! test of linear function for enthalpy (No enthalpy for mixture!!)
      if (iph == GAS_PHASE) then
         ! convert Phase molar enthalpy to (J/mol)
         f = 0.38243394899581634d-3*p + 129.1528158488924d0*T &
             - 1.2765670455401057d-6*T*p - 27991.94543799301d0
      else if (iph == LIQUID_PHASE) then
         f = 0.015156424118806432d-3*p + 0.07425560400438522d3*T &
             - 20.215678254623178d3
      end if
   end function f_MolarEnthalpy

   ! Phase molar enthalpy and its derivatives
   !< iph is an identifier for each phase
   !< P is the phase Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
#ifdef NDEBUG
   pure &
#endif
      subroutine f_MolarEnthalpy_with_derivatives(iph, p, T, C, f, dfdP, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_molar_enthalpy_with_derivatives")
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: dfdP, dfdT, dfdC(NbComp)
      real(c_double), intent(out), target :: f

      ! linear approximation from NIST (No enthalpy for mixture!!)
      if (iph == GAS_PHASE) then
         ! convert Phase molar enthalpy to (J/mol)
         f = 0.38243394899581634d-3*p + 129.1528158488924d0*T &
             - 1.2765670455401057d-6*T*p - 27991.94543799301d0
         dfdP = 0.38243394899581634d-3 - 1.2765670455401057d-6*T
         dfdT = 129.1528158488924d0 - 1.2765670455401057d-6*p
         dfdC(:) = 0.
      else if (iph == LIQUID_PHASE) then
         f = 0.015156424118806432d-3*p + 0.07425560400438522d3*T &
             - 20.215678254623178d3
         dfdP = 0.015156424118806432d-3
         dfdT = 0.07425560400438522d3
         dfdC(:) = 0.
      end if
   end subroutine f_MolarEnthalpy_with_derivatives

   subroutine f_EnergieInterne(iph, p, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      real(c_double) :: paovzeta2
      real(c_double) :: H, dPH, dTH, dCH(NbComp)
      real(c_double) :: zeta, dzetadP, dzetadT, dzetadC(NbComp)

      call f_MolarEnthalpy_with_derivatives(iph, p, T, C, H, dPH, dTH, dCH)
      call f_MolarDensity_with_derivatives(iph, p, T, C, zeta, dzetadP, dzetadT, dzetadC)
      f = H - p/zeta

      paovzeta2 = p/(zeta**2)
      dPf = dPH - 1.d0/zeta + paovzeta2*dzetadP
      dTf = dTH + paovzeta2*dzetadT
      dCf = dCH + paovzeta2*dzetadC

   end subroutine f_EnergieInterne

   !> \brief Rock volumetric heat capacity
   pure subroutine f_CpGaz(c)

      real(c_double), intent(out) :: c

      c = 1000.d0

   end subroutine

   ! Compute Psat(T)
   !< T is the Temperature
   pure subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: Psat, dT_PSat

      ! valid between -50C and 200C
      ! Psat = 100.d0*dexp(46.784d0 - 6435.d0/T - 3.868d0*dlog(T))
      ! dT_PSat = (6435.d0/T**2.d0 - 3.868d0/T)*Psat
      Psat = 1.013d5*dexp(13.7d0 - 5120.d0/T)
      dT_PSat = 5120.d0*Psat/T**2

   end subroutine FluidThermodynamics_Psat

   ! Compute Tsat(P)
   !< P is the Reference Pressure
#ifdef NDEBUG
   pure &
#endif
      subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

      real(c_double), intent(in) :: P
      real(c_double), intent(out) :: Tsat, dP_Tsat

      Tsat = 0.d0
      dP_Tsat = 0.d0

#ifndef NDEBUG
      call CommonMPI_abort('entered in Tsat, not implemented')
#endif

   end subroutine FluidThermodynamics_Tsat

#endif

end module Thermodynamics
