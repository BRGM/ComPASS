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
   use mpi, only: MPI_Abort
   use CommonMPI, only: ComPASS_COMM_WORLD
   use DefModel, only: NbPhase, NbComp, IndThermique
   
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
      f_Enthalpie

contains

   ! Fugacity
   ! iph is an identificator for each phase:
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
   subroutine f_Fugacity(rt, iph, icp, P, T, C, S, f, DPf, DTf, DCf, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph, icp
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      integer :: errcode, Ierr
      
		write(*,*) "Should never be called with a single component."
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      
   end subroutine f_Fugacity

   ! Densite molaire
   ! iph is an identificator for each phase:
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
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

      rs = 0.d0

      call FluidThermodynamics_Psat(T, Psat, dT_Psat)

      rho0 = 1.d3

      f = rho0
      dPf = 0.d0
      dTf = 0.d0

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_DensiteMolaire

   ! Densite Massique
   ! iph is an identificator for each phase:
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
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
   ! GAS_PHASE = 1; LIQUID_PHASE = 2
   subroutine f_Viscosite(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      real(c_double) :: a, ss, ds, Cs, T1

      f = 3.d-4
      dPf = 0.d0
      dTf = 0.d0

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

         f = 1.d0
         dSf = 0.d0

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

      a = -14.4319d+3
      b = 4.70915d+3
      T0 = 273.d0

      f = a + b*(T - T0)
      dPf = 0.d0
      dTf = b

      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Enthalpie

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
