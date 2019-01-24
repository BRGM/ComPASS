!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module SolvePetscWrapper

   use, intrinsic :: iso_c_binding

   use SolvePetsc
   use StringWrapper

   implicit none

   public :: &
      SolvePetscWrapper_dump_system_C

contains

   subroutine SolvePetscWrapper_dump_system_C(basename) &
      bind(C, name="SolvePetsc_dump_system")

      type(cpp_string_wrapper), intent(in) :: basename

      call SolvePetsc_dump_system(fortran_string(basename))

   end subroutine SolvePetscWrapper_dump_system_C

end module SolvePetscWrapper
