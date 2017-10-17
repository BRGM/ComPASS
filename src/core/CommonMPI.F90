!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module CommonMPI

  use mpi

  implicit none

  integer, public :: ComPASS_COMM_WORLD !< Communicator (all CPUs)

  integer, public :: &
       commRank, & !< Rank of the calling process
       commSize, & !< Total number of processors (given by MPI)
       Ncpus       !< Total number of processors

  public :: &
      CommonMPI_init

contains
  !> \brief Initialize MPI constants:
  !! ComPASS_COMM_WORLD, commSize = Ncpus, commRank
  subroutine CommonMPI_init(comm)

    integer, intent(in) :: comm ! PETSC_COMM_WORLD
    integer :: Ierr

    ! ComPASS communicator
    call MPI_Comm_dup(comm, ComPASS_COMM_WORLD, Ierr)
    if(Ierr /= MPI_SUCCESS) then
      print*, "MPI_Comm_dup error"
    end if

    ! Init CommonMPI: commRank, commSize,
    call MPI_Comm_size(ComPASS_COMM_WORLD, commSize, Ierr)
    if(Ierr /= MPI_SUCCESS) then
      print*, "MPI_Comm_size error"
    end if

    call MPI_Comm_rank(ComPASS_COMM_WORLD, commRank, Ierr)
    if(Ierr /= MPI_SUCCESS) then
      print*, "MPI_Comm_rank error"
    end if

    Ncpus = commSize

  end subroutine CommonMPI_init

end module CommonMPI
