!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncCV

    use iso_c_binding, only: c_double
    use Newton, only: &
      Newton_pointers_to_values, Newton_increments_pointers, Newton_increments
    use IncCVReservoir, only: &
      IncCVReservoir_LoadIncPreviousTimeStep, IncCVReservoir_SaveIncPreviousTimeStep, &
      IncCVReservoir_free, IncCVReservoir_allocate, IncCVReservoir_NewtonIncrement
    use IncCVWells, only: &
      IncCVWells_allocate, IncCVWells_free, IncCVWells_NewtonIncrement, &
      IncCVWells_SaveIncPreviousTimeStep, IncCVWells_LoadIncPreviousTimeStep

    implicit none


    public :: &
        IncCV_allocate, &
        IncCV_NewtonIncrement, &
        IncCV_LoadIncPreviousTimeStep, &
        IncCV_SaveIncPreviousTimeStep, &
        IncCV_free

    contains

    subroutine IncCV_allocate

    call IncCVReservoir_allocate
    call IncCVWells_allocate
    
    end subroutine IncCV_allocate

    subroutine IncCV_free

    call IncCVReservoir_free
    call IncCVWells_free

  end subroutine IncCV_free

  subroutine IncCV_NewtonIncrement_C(increments_pointers, relaxation) &
     bind(C, name="IncCV_NewtonIncrement")

    type(Newton_increments_pointers), intent(in), value :: increments_pointers
    real(c_double), intent(in), value :: relaxation
    type(Newton_increments) :: increments
    
    call Newton_pointers_to_values(increments_pointers, increments)
    call IncCV_NewtonIncrement( &
       increments%nodes, increments%fractures, increments%cells, &
       increments%injectors, increments%producers, relaxation &
    )

  end subroutine IncCV_NewtonIncrement_C


    !> \brief Newton increment of Nodes, Fracture Faces, Cells and Wells unknows.
    subroutine IncCV_NewtonIncrement( &
        NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
        NewtonIncreWellInj, NewtonIncreWellProd, relax &
    )

    real(c_double), dimension(:, :), intent(in) :: &
        NewtonIncreNode, &
        NewtonIncreFrac, &
        NewtonIncreCell
    real(c_double), dimension(:), intent(in) :: &
        NewtonIncreWellInj, &
        NewtonIncreWellProd

    real(c_double), intent(in), value :: relax

    call IncCVReservoir_NewtonIncrement(NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax)
    call IncCVWells_NewtonIncrement(NewtonIncreWellInj, NewtonIncreWellProd, relax)

    end subroutine IncCV_NewtonIncrement

    !> \brief Save current status if it is necessary to start again current time iteration.
    !!
    !! Copy IncObj to IncObjPreviousTimeStep
    subroutine IncCV_SaveIncPreviousTimeStep() &
        bind(C, name="IncCV_SaveIncPreviousTimeStep")

    call IncCVReservoir_SaveIncPreviousTimeStep
    call IncCVWells_SaveIncPreviousTimeStep

    end subroutine IncCV_SaveIncPreviousTimeStep

    !> \brief Load previous status to start again current time iteration.
    !!
    !! Copy IncObjPreviousTimeStep to IncObj
    subroutine IncCV_LoadIncPreviousTimeStep() &
        bind(C, name="IncCV_LoadIncPreviousTimeStep")

    call IncCVReservoir_LoadIncPreviousTimeStep
    call IncCVWells_LoadIncPreviousTimeStep

    end subroutine IncCV_LoadIncPreviousTimeStep

    end module IncCV
