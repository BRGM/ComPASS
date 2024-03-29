!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Thermodynamics_register

   use, intrinsic :: iso_c_binding, only: &
      c_ptr, c_funptr, c_f_procpointer, c_f_pointer
   use DefModel, only: NbPhase
   use Thermodynamics_interface, only: &
      viscosity_with_derivatives, viscosity, &
      molar_density_with_derivatives, molar_density, &
      volumetric_mass_density_with_derivatives, volumetric_mass_density, &
      molar_enthalpy_with_derivatives, molar_enthalpy

   implicit none

contains

   subroutine register_c_viscosities_with_derivatives(a) &
      bind(C, name="register_c_viscosities_with_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), viscosity_with_derivatives(k)%f)
      enddo

   end subroutine register_c_viscosities_with_derivatives

   subroutine register_c_viscosities_without_derivatives(a) &
      bind(C, name="register_c_viscosities_without_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), viscosity(k)%f)
      enddo

   end subroutine register_c_viscosities_without_derivatives

   subroutine register_c_molar_densities_with_derivatives(a) &
      bind(C, name="register_c_molar_densities_with_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), molar_density_with_derivatives(k)%f)
      enddo

   end subroutine register_c_molar_densities_with_derivatives

   subroutine register_c_molar_densities_without_derivatives(a) &
      bind(C, name="register_c_molar_densities_without_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), molar_density(k)%f)
      enddo

   end subroutine register_c_molar_densities_without_derivatives

   subroutine register_c_volumetric_mass_densities_with_derivatives(a) &
      bind(C, name="register_c_volumetric_mass_densities_with_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), volumetric_mass_density_with_derivatives(k)%f)
      enddo

   end subroutine register_c_volumetric_mass_densities_with_derivatives

   subroutine register_c_volumetric_mass_densities_without_derivatives(a) &
      bind(C, name="register_c_volumetric_mass_densities_without_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), volumetric_mass_density(k)%f)
      enddo

   end subroutine register_c_volumetric_mass_densities_without_derivatives

   subroutine register_c_molar_enthalpies_with_derivatives(a) &
      bind(C, name="register_c_molar_enthalpies_with_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), molar_enthalpy_with_derivatives(k)%f)
      enddo

   end subroutine register_c_molar_enthalpies_with_derivatives

   subroutine register_c_molar_enthalpies_without_derivatives(a) &
      bind(C, name="register_c_molar_enthalpies_without_derivatives")
      type(c_ptr), value, intent(in) :: a

      integer :: k
      type(c_funptr), pointer :: all_pointers(:)

      call c_f_pointer(a, all_pointers, [NbPhase])
      do k = 1, NbPhase
         call c_f_procpointer(all_pointers(k), molar_enthalpy(k)%f)
      enddo

   end subroutine register_c_molar_enthalpies_without_derivatives

end module Thermodynamics_register
