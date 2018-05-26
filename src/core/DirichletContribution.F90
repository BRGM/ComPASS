!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DirichletContribution

    use MeshSchema
    use DefModel
    use Thermodynamics

    use CommonMPI
    use Physics
    use SchemeParameters
    use IncCVReservoir
    use IncCVWells

    use iso_c_binding

    implicit none

    TYPE(TYPE_IncCVReservoir), allocatable, dimension(:), target, public :: &
        IncNodeDirBC !< Dirichlet boundary unknowns for current time step (size NbNodeLocal)

    public :: &
        DirichletContribution_allocate, &
        DirichletContribution_free, &
        DirichletContribution_update

     ! The following subroutines are defined in:
    ! DefInitBCvalues.F90
    public :: &
        IncCV_SetInitialValue, &
        IncCV_SetDirBCValue

    contains

       ! IncCV_SetInitialvalue and
    ! IncCV_SetDirBCvalue are defined in:
#include "DefInitBCvalues.F90"

    subroutine DirichletContribution_allocate

        integer :: Nb, Nnz

        allocate (IncNodeDirBC(NbNodeLocal_Ncpus(commRank + 1)))

    end subroutine DirichletContribution_allocate

    subroutine DirichletContribution_free

        deallocate (IncNodeDirBC)

    end subroutine DirichletContribution_free

    subroutine DirichletContribution_update

        integer :: k

        do k = 1, NbNodeLocal_Ncpus(commRank + 1)

            if (IdNodeLocal(k)%P == "d") then
            IncNode(k)%ic = IncNodeDirBC(k)%ic
            IncNode(k)%Pression = IncNodeDirBC(k)%Pression
            IncNode(k)%Saturation(:) = IncNodeDirBC(k)%Saturation(:)
            IncNode(k)%Comp(:, :) = IncNodeDirBC(k)%Comp(:, :)
            end if

#ifdef _THERMIQUE_

            if (IdNodeLocal(k)%T == "d") then
            IncNode(k)%ic = IncNodeDirBC(k)%ic
            IncNode(k)%Temperature = IncNodeDirBC(k)%Temperature
            end if

#endif

        end do

    end subroutine DirichletContribution_update

end module DirichletContribution
