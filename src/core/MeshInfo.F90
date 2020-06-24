!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module MeshInfo

   use iso_c_binding, only: c_size_t

   type, bind(C) :: PartElement
      integer(c_size_t) :: nb_owns
      integer(c_size_t) :: nb
   end type PartElement

   type, bind(C) :: PartInfo
      type(PartElement) :: nodes
      type(PartElement) :: fractures
      type(PartElement) :: injectors
      type(PartElement) :: producers
   end type PartInfo

contains

   function number_of_ghosts(element) result(n)
      type(PartElement), intent(in) :: element
      integer(c_size_t) :: n

      n = element%nb - element%nb_owns

   end function number_of_ghosts

end module MeshInfo
