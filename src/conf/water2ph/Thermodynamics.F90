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

   implicit none

   public :: &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      f_PermRel, & ! k_{r_alpha}(S)
      f_PressionCapillaire, & ! P_{c,alpha}(S)
      FluidThermodynamics_Psat, &
      FluidThermodynamics_Tsat, &
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy

contains

   subroutine f_Fugacity(rt, iph, icp, P, T, C, S, f, DPf, DTf, DCf, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph, icp
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      real(c_double) :: PSat, dTSat, Pc, DSPc(NbPhase)

      if (iph == GAS_PHASE) then
         f = P
         dPf = 1.d0
         dTf = 0.d0

      else if (iph == LIQUID_PHASE) then
         call FluidThermodynamics_Psat(T, Psat, dTSat)
         f = Psat
         dPf = 0.d0
         dTf = dTSat

#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_Fugacity')
#endif

      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Fugacity

   ! FIXME: To be removed
   subroutine check_array_interop(ina, outa) &
      bind(C, name="check_array_interop")

      real(c_double), intent(in)  :: ina(NbPhase)
      real(c_double), intent(out) :: outa(NbPhase)

      print *, "Array input", ina
      outa = ina + 0.5

   end subroutine check_array_interop

   ! FIXME #51 densite massique
   subroutine f_DensiteMolaire(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: Cs, rho0, a, b, a1, a2, b1, b2, c1, c2, cw, dcwp, dcwt
      real(c_double) :: u, R, T0, d, Z, ds, ss, rs, dZ
      real(c_double) :: Psat, dT_Psat

      if (iph == GAS_PHASE) then
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

      else if (iph == LIQUID_PHASE) then
         rs = 0.d0

         call FluidThermodynamics_Psat(T, Psat, dT_Psat)

         Cs = 0.d0 ! salinity
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

#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_DensiteMolaire')
#endif

      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_DensiteMolaire

   ! FIXME #51 densite massique
   subroutine f_DensiteMassique(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      call f_DensiteMolaire(iph, P, T, C, S, f, dPf, dTf, dCf, dSf)

   end subroutine f_DensiteMassique

   subroutine f_Viscosite(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, ss, ds, Cs, T1

      if (iph == GAS_PHASE) then
         f = (0.361d0*T - 10.2d0)*1.d-7
         dPf = 0.d0
         dTf = 0.361*1.d-7

      else if (iph == LIQUID_PHASE) then
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
         
#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_Viscosite')
#endif

      end if

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Viscosite

   ! Permeabilites = S**2
   ! iph is an identificator for each phase:
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
   subroutine f_PermRel(rt, iph, S, f, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DSf(NbPhase)

      if (iph == GAS_PHASE) then
         f = S(1)**2
         dSf(1) = 2.d0*S(1)
         dSf(2) = 0.d0

      else if (iph == LIQUID_PHASE) then
         f = S(2)**2
         dSf(1) = 0.d0
         dSf(2) = 2.d0*S(2)

#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_DensiteMolaire')
#endif

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

   ! EnergieInterne
   ! iph is an identificator for each phase:
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
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
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
   ! If Enthalpide depends on the compositon C, change DefFlash.F90
   subroutine f_Enthalpie(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, b, cc, d, T0

      if (iph == GAS_PHASE) then
         a = 1990.89d+3
         b = 190.16d+3

         f = a + T*b/100.d0
         dTf = b/100.d0

      else if (iph == LIQUID_PHASE) then
         a = -14.4319d+3
         b = 4.70915d+3
         cc = -4.87534d0
         d = 1.45008d-2
         T0 = 273.d0

         f = a + b*(T - T0) + cc*(T - T0)**2 + d*(T - T0)**3
         dTf = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2

#ifndef NDEBUG
      else
         call CommonMPI_abort('unknow phase in f_Enthalpie')
#endif

      end if

      dPf = 0.d0
      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Enthalpie

  ! Specific Enthalpy (used in FreeFlow)
  ! iph is an identificator for each phase:
  ! GAS_PHASE or LIQUID_PHASE
  subroutine f_SpecificEnthalpy(iph,P,T,C,S,f,dPf,dTf,dCf,dSf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")

    ! input
    integer(c_int), value, intent(in) :: iph
    real(c_double), value, intent(in) :: P, T
    real(c_double), intent(in) :: C(NbComp), S(NbPhase)

    ! output
    real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp), &
                                  dCf(NbComp, NbComp), dSf(NbComp, NbPhase)

      real(c_double) :: a, b, cc, d, T0

    if(iph == GAS_PHASE)then
      a = 1990.89d+3
      b = 190.16d+3

      f(:) = a + T*b/100.d0
      dTf(:) = b/100.d0

    else if(iph == LIQUID_PHASE)then
      a = -14.4319d+3
      b = 4.70915d+3
      cc = -4.87534d0
      d = 1.45008d-2
      T0 = 273.d0

      f(:) = a + b*(T - T0) + cc*(T - T0)**2 + d*(T - T0)**3
      dTf(:) = b + 2.d0*cc*(T - T0) + 3.d0*d*(T - T0)**2
    endif

    dPf = 0.d0
    dCf = 0.d0
    dSf = 0.d0

  end subroutine f_SpecificEnthalpy


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
