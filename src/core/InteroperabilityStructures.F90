!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module InteroperabilityStructures

   use, intrinsic :: iso_c_binding

   implicit none

   type, bind(C) :: cpp_array_wrapper
      type(c_ptr)       :: p
      integer(c_size_t) :: n
   end type cpp_array_wrapper

end module InteroperabilityStructures

