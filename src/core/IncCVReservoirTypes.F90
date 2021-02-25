!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncCVReservoirTypes

   use iso_c_binding, only: c_int, c_double
   use DefModel, only: NbPhase, NbComp, NbCompThermique

   !> \brief Unknown for Degree Of Freedom (including thermal).
   !! DOF can be Cell, Fracture Face or Node.
   !! if this Type is modified, mandatory to modify file wrappers/IncCV_wrappers.cpp
   !! to have the same structures in C++ and Python.
   type, bind(C) :: TYPE_IncCVReservoir
      integer(c_int) :: ic !< context: index of the set of present phase(s)
      real(c_double) :: & !< values of Inc
         Pression, & !< Reference Pressure of the element
         phase_pressure(NbPhase), &
         Temperature, & !< Temperature of the element
         Comp(NbComp, NbPhase), & !< Molar composition of the element
         Saturation(NbPhase), & !< Saturation of the element
         AccVol(NbCompThermique) !< Accumulation term integrated over volume
#ifdef _WIP_FREEFLOW_STRUCTURES_
      ! values of Inc for the soil-atmosphere boundary coundition
      real(c_double) :: FreeFlow_flowrate(NbPhase) !< molar flowrate in the freeflow (atmosphere) at the interface
#endif
   end TYPE TYPE_IncCVReservoir

!> \brief  to allow = between two TYPE_IncCVReservoir
   interface assignment(=)
      module procedure assign_type_inccv
   end interface assignment(=)

contains

   !> \brief Define operator = between two TYPE_IncCV:  inc2=inc1
   subroutine assign_type_inccv(inc2, inc1)

      type(TYPE_IncCVReservoir), intent(in) :: inc1
      type(TYPE_IncCVReservoir), intent(out) :: inc2

      inc2%ic = inc1%ic

      inc2%Pression = inc1%Pression
      inc2%phase_pressure = inc1%phase_pressure
#ifdef _THERMIQUE_
      inc2%Temperature = inc1%Temperature
#endif
      inc2%Comp = inc1%Comp
      inc2%Saturation = inc1%Saturation
      inc2%AccVol = inc1%AccVol
#ifdef _WIP_FREEFLOW_STRUCTURES_
      inc2%FreeFlow_flowrate = inc1%FreeFlow_flowrate
#endif
   end subroutine assign_type_inccv

end module IncCVReservoirTypes
