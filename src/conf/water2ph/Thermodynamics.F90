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

   use, intrinsic :: iso_c_binding

   use DefModel

   implicit none

   public :: &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      f_PermRel, & ! k_{r_alpha}(S)
      f_PressionCapillaire, & ! P_{c,alpha}(S)
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat

#ifdef _THERMIQUE_
   public :: &
      f_EnergieInterne, &
      f_Enthalpie
#endif

contains

   ! Fugacity
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_Fugacity(rt, iph, icp, P, T, C, S, f, DPf, DTf, DCf, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph, icp
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      real(c_double) :: PSat, dTSat, Pc, DSPc(NbPhase)

      if (iph == 1) then

         f = P
         dPf = 1.d0
         dTf = 0.d0

      else if (iph == 2) then

         call FluidThermodynamics_Psat(T, Psat, dTSat)
         f = Psat
         dPf = 0.d0
         dTf = dTSat
      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Fugacity

   subroutine check_array_interop(ina, outa) &
      bind(C, name="check_array_interop")

      real(c_double), intent(in)  :: ina(NbPhase)
      real(c_double), intent(out) :: outa(NbPhase)

      print *, "Array input", ina
      outa = ina + 0.5

   end subroutine check_array_interop

   ! Densite molaire
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_DensiteMolaire(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_density")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: Cs, rho0, a, b, a1, a2, b1, b2, c1, c2, cw, dcwp, dcwt
      real(c_double) :: u, R, T0, d, Z, ds, ss, rs, dZ
      real(c_double) :: Psat, dT_Psat

      if (iph == 1) then

         u = 0.018016d0
         R = 8.3145d0
         T0 = 273.d0

         a = 1.102d0
         b = -4.461d-4

         Z = 1.d0
         dZ = 0.d0

         f = P*u/(R*T*Z)
         dPf = u/(R*T*Z)
         dTf = -P*u/(R*T*Z)**2*(R*T*dZ + R*Z)

      else if (iph == 2) then

         rs = 0.d0

         call FluidThermodynamics_Psat(T, Psat, dT_Psat)

         Cs = 35.d0
         rho0 = 780.83795d0
         a = 1.6269192d0
         b = -3.0635410d-3

         a1 = 2.4638d-9
         a2 = 1.1343d-17
         b1 = -1.2171d-11
         b2 = 4.8695d-20
         c1 = 1.8452d-14
         c2 = -5.9978d-23

         ss = (rho0 + a*T + b*T**2)*(1.d0 + 6.51d-4*Cs)
         ds = (a + b*T*2.d0)*(1.d0 + 6.51d-4*Cs)

         cw = (1.d0 + 5.d-2*rs) &
              *(a1 + a2*P + T*(b1 + b2*P) + T**2*(c1 + c2*P))

         dcwp = (1.d0 + 5.d-2*rs)*(a2 + T*b2 + T**2*c2)
         dcwt = (1.d0 + 5.d-2*rs)*((b1 + b2*P) + T*2.d0*(c1 + c2*P))

         f = ss*(1.d0 + cw*(P - Psat))
         dPf = ss*dcwp*(P - Psat) + ss*cw
         dTf = ds*(1.d0 + cw*(P - Psat)) + ss*dcwt*(P - Psat) - ss*cw*dT_Psat
      else
         print *, "densite molaire error: f_DensiteMolaire"
      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_DensiteMolaire

   ! Densite Massique
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_DensiteMassique(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      call f_DensiteMolaire(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

   end subroutine f_DensiteMassique

   ! Viscosities
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_Viscosite(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, ss, ds, Cs, T1

      ! outputs
      if (iph == 1) then

         f = (0.361d0*T - 10.2d0)*1.d-7
         dPf = 0.d0
         dTf = 0.361*1.d-7

      else if (iph == 2) then

         T1 = 300.d0
         Cs = 0.04d0
         a = 1.d0 + Cs*1.34d0 + 6.12d0*Cs**2
         ss = 0.021482*(T - 273.d0 - 8.435d0 &
                        + sqrt(8078.4d0 + (T - 273.d0 - 8.435d0)**2))
         ss = ss - 1.2
         ds = 0.021482d0*(1.d0 + (T - 273.d0 - 8.435d0) &
                          /sqrt(8078.4d0 + (T - 273.d0 - 8.435d0)**2))

         f = 1.d-3*a/ss
         dPf = 0.d0
         dTf = -a*1.d-3*ds/ss**2
      else
         print *, "viscosite error: f_Viscosite"
      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Viscosite

   ! Permeabilites = S**2
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_PermRel(rt, iph, S, f, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DSf(NbPhase)

      if (iph == 1) then
         f = S(1)**2
         dSf(1) = 2.d0*S(1)
         dSf(2) = 0.d0
      else if (iph == 2) then
         f = S(2)**2
         dSf(1) = 0.d0
         dSf(2) = 2.d0*S(2)
      else
         print *, "Perm Rel error"
      end if

   end subroutine f_PermRel

   ! Pressions Capillaires des Phases et leurs derivees
   subroutine f_PressionCapillaire(rt, iph, S, f, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DSf(NbPhase)

      f = 0.d0
      dSf(:) = 0.d0

   end subroutine f_PressionCapillaire

#ifdef _THERMIQUE_

   ! EnergieInterne
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   subroutine f_EnergieInterne(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      call f_Enthalpie(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

   end subroutine f_EnergieInterne

   ! Enthalpie
   ! iph is an identificator for each phase:
   ! PHASE_GAS = 1; PHASE_WATER = 2
   ! If Enthalpide depends on the compositon C, change DefFlash.F90
   subroutine f_Enthalpie(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

      ! input
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, b, cc, d, T0

      if (iph == PHASE_GAS) then
         a = 1990.89d+3
         b = 190.16d+3

         f = a + T*b/100.d0
         dPf = 0.d0
         dTf = b/100.d0

      else if (iph == PHASE_WATER) then
         a = -14.4319d+3
         b = 4.70915d+3
         cc = -4.87534d0
         d = 1.45008d-2
         T0 = 273.d0

         f = a + b*(T - T0) + cc*(T - T0)**2 + d*(T - T0)**3
         dPf = 0.d0
         dTf = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2
      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Enthalpie

#endif

   subroutine FluidThermodynamics_Psat(T, Psat, dT_PSat) &
      bind(C, name="FluidThermodynamics_Psat")

      real(c_double), value, intent(in) :: T
      real(c_double), intent(out) :: Psat, dT_PSat

      Psat = (T - 273.d0)**4.d0/1.0d3
      dT_PSat = 4.d0*(T - 273.d0)**3.d0/1.0d3

   end subroutine FluidThermodynamics_Psat

   subroutine FluidThermodynamics_Tsat(P, Tsat, dP_Tsat) &
      bind(C, name="FluidThermodynamics_Tsat")

      real(c_double), value, intent(in) :: P
      real(c_double), intent(out) :: Tsat, dP_Tsat

      Tsat = 100.d0*(P/1.d5)**0.25d0 + 273.d0
      dP_Tsat = 0.25d0*1.d-3*(P/1.d5)**(-0.75d0)

   end subroutine FluidThermodynamics_Tsat

end module Thermodynamics
