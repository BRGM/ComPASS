!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module PartitionMesh

  use CommonType
  use GlobalMesh

  use iso_c_binding, only : C_CHAR, C_NULL_CHAR

  implicit none

  ! Output
  integer, allocatable, dimension(:) :: ProcbyCell

  interface

    subroutine Metis_C ( nMaille, ptMaillebyMaille, numMaillebyMaille, npart, objval, procbyMaille ) bind ( C, name = "Metis_C" )
      use iso_c_binding, only : C_INT
      integer ( c_int ) :: ptMaillebyMaille(*), numMaillebyMaille(*)
      integer ( c_int ) :: procbyMaille(*)
      integer ( c_int ), value :: nMaille, npart, objval
    end subroutine Metis_C

  end interface

  public :: PartitionMesh_Metis, PartitionMesh_Dealloc

contains

  ! ------------------------------------------------------------------ !

  subroutine PartitionMesh_Metis(Ncpus)

    integer, intent(in) :: Ncpus
    integer :: objval, i, pt

    objval = 0

    allocate(ProcbyCell(NbCell)) ! deallocation ds deallocate global mesh
    ProcbyCell(:) = 0

    if (Ncpus==1) then
      ProcbyCell(:) = 0
    else
      call Metis_C ( NbCell, CellbyCell%Pt, CellbyCell%Num-1, Ncpus, objval, ProcbyCell )
    endif

  end subroutine PartitionMesh_Metis


  subroutine PartitionMesh_Dealloc

    deallocate(ProcbyCell)

  end subroutine PartitionMesh_Dealloc

end module PartitionMesh
