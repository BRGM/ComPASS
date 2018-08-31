!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: a single fluid phase with *linear* behavior

module Thermodynamics

   use, intrinsic :: iso_c_binding

   use DefModel
   use CommonMPI

   implicit none

   !isobar thermal expansivity (K-1) and isothermal compressibility (Pa-1)
   type, bind(c) :: fluid_properties_type
      real(c_double) :: specific_mass
      real(c_double) :: compressibility
      real(c_double) :: thermal_expansivity
      real(c_double) :: volumetric_heat_capacity
      real(c_double) :: dynamic_viscosity
   end type

   ! OPTIMIZE: is there a loss of performance uing the target keyword?
   type(fluid_properties_type), target :: fluid_properties

   public :: &
      get_fluid_properties, &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      f_PermRel, & ! k_{r_alpha}(S)
      f_PressionCapillaire, & ! P_{c,alpha}(S)
      f_EnergieInterne, &
      f_Enthalpie

contains

   function get_fluid_properties() result(properties_p) &
      bind(C, name="get_fluid_properties")

      type(c_ptr) :: properties_p

      properties_p = c_loc(fluid_properties)

   end function get_fluid_properties

   subroutine f_Fugacity(rt, iph, icp, P, T, C, S, f, DPf, DTf, DCf, DSf)

      ! input
      integer(c_int), intent(in) :: rt(IndThermique + 1)
      integer(c_int), intent(in) :: iph, icp
      real(c_double), intent(in) :: P, T, C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, DPf, DTf, DCf(NbComp), DSf(NbPhase)

      integer :: errcode, Ierr

      write (*, *) "Should never be called with a single component."
      call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)

   end subroutine f_Fugacity

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

      f = fluid_properties%specific_mass
      dPf = fluid_properties%specific_mass*fluid_properties%compressibility
      dTf = fluid_properties%specific_mass*fluid_properties%thermal_expansivity
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

   subroutine f_Viscosite(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: &
         f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      f = fluid_properties%dynamic_viscosity
      dPf = 0.d0
      dTf = 0.d0

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
   subroutine f_Enthalpie(iph, P, T, C, S, f, dPf, dTf, dCf, dSf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")

      ! input
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp), S(NbPhase)

      ! output
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp), dSf(NbPhase)

      f = fluid_properties%volumetric_heat_capacity * T
      dPf = 0.d0
      dTf = fluid_properties%volumetric_heat_capacity
      dCf(:) = 0.d0
      dSf(:) = 0.d0

   end subroutine f_Enthalpie

end module Thermodynamics
