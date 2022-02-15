!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module FreeFlowTypes

   use iso_c_binding, only: c_double
   use DefModel, only: NbPhase, NbComp

   !> \brief Storage for Freeflow values
   !! DOF are Node.
   !! if this Type is modified, mandatory to modify file wrappers/FreeFlow_wrappers.cpp
   !! and src/core/include/StateObjects.h
   !! to have the same structures in C++ and Python.
   type, bind(C) :: TYPE_FFfarfield
      real(c_double) :: & !< values of
         Pressure, & !< gas pressure
         Temperature(NbPhase), & !< gas far-field and liquid (rain) Temperature
         Comp(NbComp, NbPhase), & !< Molar composition of the element
         Imposed_flux(NbPhase), & !< imposed flux at the boundary (for the rain for example)
         Hm(NbPhase), & !< gas (liquid = 0) Convective-diffusive const of the FreeFlow boundary layer
         HT !< Thermal Convective const of the FreeFlow boundary layer
   end TYPE TYPE_FFfarfield

!> \brief  to allow = between two TYPE_FFfarfield
   interface assignment(=)
      module procedure assign_type_FFstate
   end interface assignment(=)

contains

   !> \brief Define operator = between two TYPE_FFfarfield:  ff2=ff1
   subroutine assign_type_FFstate(ff2, ff1)

      type(TYPE_FFfarfield), intent(in) :: ff1
      type(TYPE_FFfarfield), intent(out) :: ff2

      ff2%Pressure = ff1%Pressure
      ff2%Temperature = ff1%Temperature
      ff2%Comp = ff1%Comp
      ff2%Imposed_flux = ff1%Imposed_flux
   end subroutine assign_type_FFstate

end module FreeFlowTypes
