!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module PartitionMesh

  use iso_c_binding, only : c_char, c_null_char, c_int
  use CommonType, only: CSR
  use CommonMPI, only: CommonMPI_abort

  implicit none

  interface

    subroutine Metis_C ( nMaille, ptMaillebyMaille, numMaillebyMaille, npart, objval, procbyMaille ) bind ( C, name = "Metis_C" )
      use iso_c_binding, only : c_int
      integer ( c_int ) :: ptMaillebyMaille(*), numMaillebyMaille(*)
      integer ( c_int ) :: procbyMaille(*)
      integer ( c_int ), value :: nMaille, npart, objval
    end subroutine Metis_C

  end interface

  public :: PartitionMesh_Metis

contains

  subroutine PartitionMesh_Metis(Ncpus, CellbyCell, ProcbyCell)

    integer, intent(in) :: Ncpus
    type(CSR), intent(in) :: CellbyCell
    integer, allocatable, dimension(:), intent(inout) :: ProcbyCell
    integer :: NbCell, objval

    NbCell = CellbyCell%Nb
    objval = 0
    
    if(.not.allocated(ProcbyCell)) call CommonMPI_abort('ProcbyCell should be allocated')
    if(size(ProcbyCell)/=NbCell) call CommonMPI_abort('inconsistent numer of cells')

    ProcbyCell(:) = 0

    if (Ncpus==1) then
      ProcbyCell(:) = 0
    else
      call Metis_C ( NbCell, CellbyCell%Pt, CellbyCell%Num-1, Ncpus, objval, ProcbyCell )
    endif

  end subroutine PartitionMesh_Metis

end module PartitionMesh
