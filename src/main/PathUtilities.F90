!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module PathUtilities

      implicit none

      public :: &
        make_directory

    contains

      subroutine make_directory(path)

        character(len=*), intent(in) :: path
        character(len(path) + 10)      :: cmd

        write (cmd, '(A)') "mkdir -p "//trim(path)
        call system(cmd)

      end subroutine make_directory

    end module PathUtilities
