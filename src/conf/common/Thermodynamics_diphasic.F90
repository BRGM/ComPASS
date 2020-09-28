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
                       get_model_configuration
   use CommonType, only: ModelConfiguration
   use RelativePermeabilities, only: f_PermRel
   use CapillaryPressure, only: f_PressionCapillaire

   implicit none

   public :: &
      f_Fugacity, &  !< Fucacity
      f_DensiteMolaire, &  !< \xi^alpha(P,T,C,S)
      f_DensiteMassique, &  !< \rho^alpha(P,T,C,S)
      f_Viscosite, &  !< \mu^alpha(P,T,C,S)
      air_henry, &  !< Henry coef for air comp
      air_henry_dT, &  !< derivative of the Henry coef for air comp
      liquid_pressure          !< compute liquid pressure with hydrostatic pressure

#ifdef _THERMIQUE_
   public :: &
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy, &
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat
#endif

contains

   pure function liquid_pressure(z_ref, p_ref, rho, g, z)
      real(c_double), intent(in) :: z_ref
      real(c_double), intent(in) :: p_ref
      real(c_double), intent(in) :: rho
      real(c_double), intent(in) :: g
      real(c_double), intent(in) :: z
      real(c_double) :: liquid_pressure

      liquid_pressure = p_ref - rho*g*(z - z_ref)
   end function liquid_pressure

   ! *** Physics *** !

   ! Fugacity coefficient
   !< rt is the rocktype identifier
   !< iph is the phase identifier : GAS_PHASE or LIQUID_PHASE
   !< icp component identifier
   !< P is the reference pressure
   !< T is the temperature
   !< C is the phase molar frcations
   !< S is all the saturations
   pure subroutine f_Fugacity(rt, iph, icp, P, T, C, S, f, DPf, DTf, DCf, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph, icp
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      real(c_double) :: PSat, dTSat, Pc, DSPc(NbPhase)
      real(c_double) :: RZetal

      RZetal = 8.314d0*1000.d0/0.018d0

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0

      if (iph == GAS_PHASE) then
         call f_PressionCapillaire(rt, iph, S, Pc, DSPc)  ! Pg=Pref + Pc
         f = P + Pc

         dPf = 1.d0
         dSf = DSPc
      else if (iph == LIQUID_PHASE) then
         if (icp == AIR_COMP) then
            call air_henry(T, f)
            call air_henry_dT(dTf)
         else if (icp == WATER_COMP) then
            call f_PressionCapillaire(rt, iph, S, Pc, DSPc)
            ! FIXME: Pl = Pref + f_PressionCapillaire, so Pc = -Pc
            Pc = -Pc
            DSPc = -DSPc
            call FluidThermodynamics_Psat(T, Psat, dTSat)

            f = Psat*dexp(Pc/(T*RZetal))

            dTf = (dTSat - Psat*Pc/RZetal/(T**2))*dexp(Pc/(T*RZetal))
            dSf = DSPc*f/(T*RZetal)
         endif
      endif
   end subroutine f_Fugacity

   !> \brief Henry coef for air comp
   pure subroutine air_henry(T, H)

      real(c_double), intent(in) :: T
      real(c_double), intent(out) :: H

      real(c_double) :: T1, T2
      real(c_double) :: H1, H2

      T1 = 293.d0
      T2 = 353.d0

      H1 = 6.d+9
      H2 = 10.d+9

      H = H1 + (H2 - H1)*(T - T1)/(T2 - T1)

   end subroutine

   !> \brief Derivative of the Henry coef for air comp wrt Temperature
   pure subroutine air_henry_dT(H_dt)

      real(c_double), intent(out) :: H_dt

      real(c_double) :: T1, T2
      real(c_double) :: H1, H2

      T1 = 293.d0
      T2 = 353.d0

      H1 = 6.d+9
      H2 = 10.d+9

      H_dt = (H2 - H1)/(T2 - T1)

   end subroutine

   ! Molar density
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
#ifdef NDEBUG
   pure &
#endif
      subroutine f_DensiteMolaire(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_density")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: Rgp, Pc, DSPc(NbPhase), Pg
      integer(c_int) :: rt(IndThermique + 1)

      Rgp = 8.314d0

      if (iph == GAS_PHASE) then
         rt = 0 ! FIXME: rt is not used because Pref=Pg so Pc=0.
         call f_PressionCapillaire(rt, iph, S, Pc, DSPc)
#ifndef NDEBUG
         if (Pc .ne. 0.d0) &
            call CommonMPI_abort('possible error in f_DensiteMolaire (change rt)')
#endif
         Pg = P + Pc
         f = Pg/(Rgp*T)

         dPf = 1/(Rgp*T)
         dTf = -Pg/Rgp/T**2
         dCf = 0.d0
         dSf = DSPc(iph)/(Rgp*T)
      else if (iph == LIQUID_PHASE) then
         f = 1000.d0/M_H2O

         dPf = 0.d0
         dTf = 0.d0
         dCf = 0.d0
         dSf = 0.d0
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
      subroutine f_DensiteMassique(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: zeta, comp_m(NbComp), m
      integer :: icp
      type(ModelConfiguration) :: configuration
      configuration = get_model_configuration()

      call f_DensiteMolaire(iph, P, T, C, S, zeta, dPf, dTf, dCf, dSf)   ! P is Reference Pressure

      comp_m(AIR_COMP) = M_air
      comp_m(WATER_COMP) = M_H2O

      m = 0.d0
      ! loop of component in phase iph
      do icp = 1, NbComp
         ! configuration%MCP(icp,iph)=1 if component in phase iph
         m = m + configuration%MCP(icp, iph)*C(icp)*comp_m(icp)
      enddo

      f = zeta*m

      dPf = dPf*m
      dTf = dTf*m
      ! loop of component in phase iph
      do icp = 1, NbComp
         dCf(icp) = dCf(icp)*m + zeta*configuration%MCP(icp, iph)*comp_m(icp)
      enddo
      dSf = dSf*m

   end subroutine f_DensiteMassique

   ! Viscosities
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   pure subroutine f_Viscosite(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      if (iph == GAS_PHASE) then
         f = 2.d-5
      else if (iph == LIQUID_PHASE) then
         f = 1.d-3
      endif

      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0
      dSf = 0.d0

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
      subroutine f_EnergieInterne(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: Piph, Pc, DSPc(NbPhase)
      real(c_double) :: zeta, dzetadP, dzetadT, dzetadC(NbComp), dzetadS(NbPhase)
      real(c_double) :: enth, denthdP, denthdT, denthdC(NbComp), denthdS(NbPhase)
      integer(c_int) :: rt(IndThermique + 1)

      rt = 0
      call f_PressionCapillaire(rt, iph, S, Pc, DSPc)
      Piph = P + Pc

      CALL f_Enthalpie(iph, P, T, C, S, enth, denthdP, denthdT, denthdC, denthdS) ! called with reference pressure
      CALL f_DensiteMolaire(iph, P, T, C, S, zeta, dzetadP, dzetadT, dzetadC, dzetadS) ! called with reference pressure
      f = enth - Piph/zeta ! P is phase pressure
      dPf = denthdP - 1.d0/zeta + Piph/zeta**2*dzetadP
      dTf = denthdT + Piph/zeta**2*dzetadT
      dCf = denthdC + Piph/zeta**2*dzetadC
      dSf = denthdS - DSPc/zeta + Piph/zeta**2*dzetadS

   end subroutine f_EnergieInterne

   ! Enthalpie
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   pure subroutine f_Enthalpie(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, b, cc, d, Ts, T0, ss, dTss, cp
      type(ModelConfiguration) :: configuration

      if (iph == GAS_PHASE) then
         a = 1990.89d+3
         b = 190.16d+3
         cc = -1.91264d+3
         d = 0.2997d+3

         Ts = T/100.d0

         ss = a + b*Ts + cc*Ts**2.d0 + d*Ts**3.d0
         dTss = (b + 2.d0*cc*Ts + 3.d0*d*Ts**2.d0)/100.d0

         call f_CpGaz(cp)

         configuration = get_model_configuration()

         f = configuration%MCP(AIR_COMP, iph)*C(AIR_COMP)*cp*M_air*T + &
             configuration%MCP(WATER_COMP, iph)*C(WATER_COMP)*ss*M_H2O

         dPf = 0.d0
         dTf = configuration%MCP(AIR_COMP, iph)*C(AIR_COMP)*cp*M_air + &
               configuration%MCP(WATER_COMP, iph)*C(WATER_COMP)*M_H2O*dTss
         dCf(AIR_COMP) = configuration%MCP(AIR_COMP, iph)*cp*M_air*T
         dCf(WATER_COMP) = configuration%MCP(WATER_COMP, iph)*ss*M_H2O
         dSf = 0.d0

      else if (iph == LIQUID_PHASE) then

         a = -14.4319d+3
         b = +4.70915d+3
         cc = -4.87534
         d = 1.45008d-2
         T0 = 273.d0

         ss = a + b*(T - T0) + cc*(T - T0)**2.d0 + d*(T - T0)**3.d0
         dTss = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2.d0

         f = ss*M_H2O

         dPf = 0.d0
         dTf = dTss*M_H2O
         dCf = 0.d0
         dSf = 0.d0
      endif

   end subroutine f_Enthalpie

   ! Specific Enthalpy (used in FreeFlow)
   !< iph is an identifier for each phase: GAS_PHASE or LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   pure subroutine f_SpecificEnthalpy(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp), &
                                     dCf(NbComp, NbComp), dSf(NbComp, NbPhase)

      real(c_double) :: a, b, cc, d, Ts, T0, ss, dTss, cp
      type(ModelConfiguration) :: configuration

      if (iph == GAS_PHASE) then
         a = 1990.89d+3
         b = 190.16d+3
         cc = -1.91264d+3
         d = 0.2997d+3

         Ts = T/100.d0

         ss = a + b*Ts + cc*Ts**2.d0 + d*Ts**3.d0
         dTss = (b + 2.d0*cc*Ts + 3.d0*d*Ts**2.d0)/100.d0

         call f_CpGaz(cp)

         configuration = get_model_configuration()

         f(AIR_COMP) = configuration%MCP(AIR_COMP, iph)*cp*M_air*T
         f(WATER_COMP) = configuration%MCP(WATER_COMP, iph)*ss*M_H2O

         dPf = 0.d0
         dTf(AIR_COMP) = configuration%MCP(AIR_COMP, iph)*cp*M_air
         dTf(WATER_COMP) = configuration%MCP(WATER_COMP, iph)*M_H2O*dTss
         dCf = 0.d0
         dSf = 0.d0

      else if (iph == LIQUID_PHASE) then

         a = -14.4319d+3
         b = +4.70915d+3
         cc = -4.87534
         d = 1.45008d-2
         T0 = 273.d0

         ss = a + b*(T - T0) + cc*(T - T0)**2.d0 + d*(T - T0)**3.d0
         dTss = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2.d0

         f(:) = ss*M_H2O

         dPf = 0.d0
         dTf(:) = dTss*M_H2O
         dCf = 0.d0
         dSf = 0.d0
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
