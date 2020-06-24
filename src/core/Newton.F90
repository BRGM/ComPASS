!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Newton

   use, intrinsic :: iso_c_binding, only: c_double, c_ptr, c_f_pointer
   use CommonMPI, only: commRank

   use DefModel, only: NbIncTotalMax
   use MeshSchema, only: &
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus, NbCellLocal_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus

   implicit none

   type, bind(C) :: Newton_increments_pointers
      type(c_ptr) :: nodes
      type(c_ptr) :: fractures
      type(c_ptr) :: cells
      type(c_ptr) :: injectors
      type(c_ptr) :: producers
   end type Newton_increments_pointers

   type Newton_increments
      real(c_double), pointer :: nodes(:, :)
      real(c_double), pointer :: fractures(:, :)
      real(c_double), pointer :: cells(:, :)
      real(c_double), pointer :: injectors(:)
      real(c_double), pointer :: producers(:)
   end type Newton_increments

   public :: &
      Newton_pointers_to_values

contains

   subroutine Newton_pointers_to_values(increment_pointers, increment_values)

      type(Newton_increments_pointers), intent(in), value :: increment_pointers
      type(Newton_increments), intent(out) :: increment_values

      call c_f_pointer( &
         increment_pointers%nodes, increment_values%nodes, &
         shape=[NbIncTotalMax, NbNodeLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%fractures, increment_values%fractures, &
         shape=[NbIncTotalMax, NbFracLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%cells, increment_values%cells, &
         shape=[NbIncTotalMax, NbCellLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%injectors, increment_values%injectors, &
         shape=[NbWellInjLocal_Ncpus(commRank + 1)] &
         )
      call c_f_pointer( &
         increment_pointers%producers, increment_values%producers, &
         shape=[NbWellProdLocal_Ncpus(commRank + 1)] &
         )

   end subroutine Newton_pointers_to_values

end module Newton
