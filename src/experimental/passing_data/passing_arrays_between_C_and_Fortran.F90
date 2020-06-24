!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Note: well is only consided in cpramg with pressure
!       i.e. SolvePetsc_cpramgPCApply_P_multiplicative
!       well is not added in SolvePetsc_At

! This model is just for testing

module Foo

   use, intrinsic :: iso_c_binding

   implicit none

   public :: &
      pass_and_dump_dim_array, &
      increment_first_column

contains

   subroutine pass_and_dump_dim_array(cp, n, scp) &
      bind(C, name="pass_and_dump_dim_array")

      type(c_ptr), intent(in), value :: cp, scp
      integer(c_int), intent(in), value :: n
      real(c_double), pointer :: p(:, :)
      integer(c_int), pointer :: s(:)
      integer :: i

      call c_f_pointer(scp, s, [n])
      call c_f_pointer(cp, p, s)
      do i = 1, size(p, 2)
         write (*, *) p(:, i)
      end do

   end subroutine pass_and_dump_dim_array

   subroutine increment_first_column(cp, n, scp) &
      bind(C, name="increment_first_column")

      type(c_ptr), intent(in), value :: cp, scp
      integer(c_int), intent(in), value :: n
      real(c_double), pointer :: p(:, :)
      integer(c_int), pointer :: s(:)
      integer :: i

      call c_f_pointer(scp, s, [n])
      call c_f_pointer(cp, p, s)
      p(1, :) = p(1, :) + 1

   end subroutine increment_first_column

end module Foo
