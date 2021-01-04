!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DirichletContribution

   use CommonMPI, only: commRank
   use IncCVReservoirTypes, only: TYPE_IncCVReservoir
   use IncCVReservoir, only: IncNode
   use MeshSchema, only: IdNodeLocal, NbNodeLocal_Ncpus

   implicit none

   TYPE(TYPE_IncCVReservoir), allocatable, dimension(:), target, public :: &
      IncNodeDirBC !< Dirichlet boundary unknowns for current time step (size NbNodeLocal)

   public :: &
      DirichletContribution_allocate, &
      DirichletContribution_free, &
      DirichletContribution_update

contains

   !> \brief Allocate vector with Dirichlet contributions
   subroutine DirichletContribution_allocate

      allocate (IncNodeDirBC(NbNodeLocal_Ncpus(commRank + 1)))

   end subroutine DirichletContribution_allocate

   !> \brief Deallocate vector with Dirichlet contributions
   subroutine DirichletContribution_free

      deallocate (IncNodeDirBC)

   end subroutine DirichletContribution_free

   subroutine DirichletContribution_update() &
      bind(C, name="DirichletContribution_update")

      integer :: k

      do k = 1, NbNodeLocal_Ncpus(commRank + 1)

         ! Can not use "=" of the two structures
         ! because Dirichlet can be on Pressure only or Temperature only
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
