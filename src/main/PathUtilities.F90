
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
