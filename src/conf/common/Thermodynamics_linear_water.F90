!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Model: a single fluid phase with *linear* behavior

module Thermodynamics

   use, intrinsic :: iso_c_binding, only: c_double, c_int, c_loc, c_ptr
   use mpi, only: MPI_Abort
   use CommonMPI, only: ComPASS_COMM_WORLD
   use DefModel, only: NbPhase, NbComp, IndThermique
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir

#ifndef NDEBUG
   use CommonMPI, only: CommonMPI_abort
#endif

   implicit none

   !isobar thermal expansivity (K-1) and isothermal compressibility (Pa-1)
   type, bind(c) :: fluid_properties_type
      real(c_double) :: specific_mass
      real(c_double) :: compressibility
      real(c_double) :: thermal_expansivity
      real(c_double) :: volumetric_heat_capacity
      real(c_double) :: dynamic_viscosity
      real(c_double) :: reference_pressure
      real(c_double) :: reference_temperature
   end type

   ! OPTIMIZE: is there a loss of performance uing the target keyword?
   type(fluid_properties_type), target :: fluid_properties

   public :: &
      get_fluid_properties, &
      f_Fugacity, & ! Fucacity
      f_DensiteMolaire, & ! \xi^alpha(P,T,C,S)
      f_DensiteMassique, & ! \rho^alpha(P,T,C,S)
      f_Viscosite, & ! \mu^alpha(P,T,C,S)
      f_EnergieInterne, &
      f_Enthalpie, &
      f_SpecificEnthalpy

contains

   function get_fluid_properties() result(properties_p) &
      bind(C, name="get_fluid_properties")

      type(c_ptr) :: properties_p

      properties_p = c_loc(fluid_properties)

   end function get_fluid_properties

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
      subroutine f_Fugacity(icp, iph, p, T, C, f, dfdp, dfdT, dfdC) &
      bind(C, name="FluidThermodynamics_fugacity")
      integer(c_int), intent(in) :: icp, iph
      real(c_double), intent(in) :: p, T, C(NbComp)
      real(c_double), intent(out) :: f, dfdp, dfdT, dfdC(NbComp)

#ifndef NDEBUG
      call CommonMPI_abort("Fugacity should never be called with a single component.")
#endif

   end subroutine f_Fugacity

   ! Densite molaire
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_DensiteMolaire(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_density")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      f = fluid_properties%specific_mass &
          *exp(fluid_properties%compressibility*(P - fluid_properties%reference_pressure)) &
          *exp(fluid_properties%thermal_expansivity*(T - fluid_properties%reference_temperature))
      dPf = f*fluid_properties%compressibility
      dTf = f*fluid_properties%thermal_expansivity
      dCf = 0.d0

   end subroutine f_DensiteMolaire

   ! Densite Massique
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_DensiteMassique(iph, P, T, C, f, dPf, dTf, dCf)
      integer(c_int), intent(in) :: iph
      real(c_double), intent(in) :: P, T, C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_DensiteMolaire(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_DensiteMassique

   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_Viscosite(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_dynamic_viscosity")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      f = fluid_properties%dynamic_viscosity
      dPf = 0.d0
      dTf = 0.d0
      dCf = 0.d0

   end subroutine f_Viscosite

   ! EnergieInterne
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_EnergieInterne(iph, P, T, C, f, dPf, dTf, dCf)
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_Enthalpie(iph, P, T, C, f, dPf, dTf, dCf)

   end subroutine f_EnergieInterne

   subroutine f_enthalpy_proxy(T, f, dTf)
      real(c_double), value, intent(in) :: T
      real(c_double), intent(out) :: f, dTf

      f = fluid_properties%volumetric_heat_capacity*T
      dTf = fluid_properties%volumetric_heat_capacity

   end subroutine f_enthalpy_proxy

   ! Enthalpie
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_Enthalpie(iph, P, T, C, f, dPf, dTf, dCf) &
      bind(C, name="FluidThermodynamics_molar_enthalpy")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(in) :: C(NbComp)
      real(c_double), intent(out) :: f, dPf, dTf, dCf(NbComp)

      call f_enthalpy_proxy(T, f, dTf)
      dPf = 0.d0
      dCf = 0.d0

   end subroutine f_Enthalpie

   ! Specific Enthalpy (used in FreeFlow)
   !< iph is an identifier for each phase, here only one phase: LIQUID_PHASE
   !< P is the Reference Pressure
   !< T is the Temperature
   !< C is the phase molar fractions
   !< S is all the saturations
   subroutine f_SpecificEnthalpy(iph, P, T, f, dPf, dTf) &
      bind(C, name="FluidThermodynamics_molar_specific_enthalpy")
      integer(c_int), value, intent(in) :: iph
      real(c_double), value, intent(in) :: P, T
      real(c_double), intent(out) :: f(NbComp), dPf(NbComp), dTf(NbComp)

      real(c_double) :: fv, dTfv

#ifndef NDEBUG
      if (NbComp /= 1) &
         call CommonMPI_abort("f_SpecificEnthalpy: single component assumed!")
#endif

      call f_enthalpy_proxy(T, fv, dTfv)
      f = fv
      dPf = 0.d0
      dTf = dTfv

   end subroutine f_SpecificEnthalpy

end module Thermodynamics
