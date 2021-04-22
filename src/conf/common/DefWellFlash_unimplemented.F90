!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefWellFlash

   use CommonMPI, only: CommonMPI_abort
   use IncCVWells, only: NodebyWellProdLocal

   implicit none

   public :: DefWellFlash_TimeFlashWellProd

contains

   subroutine DefWellFlash_TimeFlashWellProd() &
      bind(C, name="DefWellFlash_TimeFlashWellProd")

      integer :: nWell

      ! does nothing if there are no production wells
      do nWell = 1, NodebyWellProdLocal%Nb
         call CommonMPI_abort("Well flash is not implemented for this physics.")
      end do

   end subroutine DefWellFlash_TimeFlashWellProd

end module DefWellFlash
