!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefMSWell
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use mpi, only: &
      MPI_CHARACTER, &
      MPI_INT, &
      MPI_DOUBLE, &
      MPI_Type_Create_Struct, &
      MPI_ADDRESS_KIND, &
      MPI_Abort

   use mpi_f08, only: &
      MPI_Get_address

   use, intrinsic :: iso_c_binding, only: &
      c_ptr, c_size_t, c_null_ptr, c_loc, c_int, c_double, c_char, c_bool

   use CommonType, only: CSR
   use CommonMPI, only: ComPASS_COMM_WORLD, CommonMPI_abort
   use DefModel
   use Physics
   use DefWell

#else

   use, intrinsic :: iso_c_binding, only: &
      c_ptr, c_size_t, c_null_ptr, c_loc, c_int, c_double, c_char, c_bool

   use mpi, only: &
      MPI_CHARACTER, &
      MPI_INT, &
      MPI_DOUBLE, &
      MPI_Type_Create_Struct, &
      MPI_ADDRESS_KIND, &
      MPI_Abort

   ! could also be `use petscmpi` in latest versions of PETSc
   use petscsys, only: &
      MPI_Get_address

   use CommonType, only: CSR
   use CommonMPI, only: ComPASS_COMM_WORLD, CommonMPI_abort
   use DefModel, only: NbComp
   use Physics, only: thickness
   use DefWell, only: &
      TYPE_CSRDataNodeWell, DefWell_WellIndex
#endif

   implicit none

   type, bind(c) :: MSWellData_type
      integer(c_int) :: Id
      real(c_double) :: &
         Radius, & ! both well types
         PressionMax, & ! injector only
         PressionMin, & ! producer only
         ImposedFlowrate, & ! both well types (>=0 for producer <0 for injector)
         CompTotal(NbComp), & ! injector only
         InjectionTemperature ! injector only
      ! WARNING: we put character at the end of the structure
      ! because of "memory padding" when creating mpi well data structure
      ! cf. DefWell_mpi_register_well_data_description
      character(c_char) :: &
         IndWell, & ! both well types 'p' for pressure mode ; 'f' for flowrate mode; 'c' for closed
         Well_type  !< 'p' for producer, 'i' for injector.
   end type MSWellData_type

   ! FIXME: switch to pointer
   type(MSWellData_type), allocatable, public, target, dimension(:) :: &
      DataMSWell

   type(TYPE_CSRDataNodeWell), public :: &
      NodeDatabyMSWell !< CSR store data about Parent and Well index of nodes of each MSWell
   public:: &
      DefMSWell_Make_ComputeWellIndex, &
      DefMSWell_mpi_register_mswell_data_description

contains

   subroutine DefMSWell_Make_ComputeWellIndex( &
      NbNode, XNode, CellbyNode, NodebyCell, &
      FracbyNode, NodebyFace, PermCell, PermFrac)

      integer, intent(in) :: NbNode
      double precision, allocatable, dimension(:, :), intent(in) :: XNode
      type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
      !type(FractureInfoCOC), intent(in) :: FracbyNode
      type(CSR), intent(in) :: FracbyNode

      double precision, allocatable, dimension(:, :, :), intent(in) :: PermCell
      double precision, allocatable, dimension(:), intent(in) :: PermFrac

      integer :: NbMSWell
      double precision, allocatable, dimension(:) :: WellRadius
      ! FIXME: set consistent values to error codes
      integer :: errcode, Ierr

      if (.NOT. allocated(DataMSWell)) then
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         write (*, *) "ERROR DataMSWell is not allocated."
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
      NbMSWell = size(DataMSWell)

      allocate (WellRadius(NbMSWell))

      WellRadius(:) = 0
      WellRadius(1:NbMSWell) = DataMSWell(:)%Radius
      call DefWell_WellIndex(NodeDatabyMSWell, NbMSWell, WellRadius, &
                             NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, &
                             NodebyFace, PermCell, PermFrac)

      deallocate (WellRadius)
   end subroutine DefMSWell_Make_ComputeWellIndex

   subroutine DefMSWell_mpi_register_mswell_data_description(mpi_id)

      integer, intent(out) :: mpi_id

      integer, parameter :: count = 9
      integer :: blocklengths(count)
      integer(kind=MPI_ADDRESS_KIND) :: begin, offset, displacements(count)
      integer :: types(count)
      type(MSWellData_type) :: dummy
      integer :: Ierr

      call MPI_Get_address(dummy, begin, Ierr)
      call MPI_Get_address(dummy%Id, offset, Ierr)
      displacements(1) = offset - begin
      call MPI_Get_address(dummy%Radius, offset, Ierr)
      displacements(2) = offset - begin
      call MPI_Get_address(dummy%PressionMax, offset, Ierr)
      displacements(3) = offset - begin
      call MPI_Get_address(dummy%PressionMin, offset, Ierr)
      displacements(4) = offset - begin
      call MPI_Get_address(dummy%ImposedFlowrate, offset, Ierr)
      displacements(5) = offset - begin
      call MPI_Get_address(dummy%CompTotal, offset, Ierr)
      displacements(6) = offset - begin
      call MPI_Get_address(dummy%InjectionTemperature, offset, Ierr)
      displacements(7) = offset - begin
      call MPI_Get_address(dummy%IndWell, offset, Ierr)
      displacements(8) = offset - begin
      call MPI_Get_address(dummy%Well_type, offset, Ierr)
      displacements(9) = offset - begin

      types(:) = MPI_DOUBLE
      types(1) = MPI_INT
      types(count - 1) = MPI_CHARACTER
      types(count) = MPI_CHARACTER

      blocklengths(:) = 1
      blocklengths(6) = NbComp

      call MPI_Type_Create_Struct(count, blocklengths, displacements, types, mpi_id, Ierr)
      if (Ierr /= 0) call CommonMPI_abort('Couldt not create well data MPI structure.')
      call MPI_Type_commit(mpi_id, Ierr)
      if (Ierr /= 0) call CommonMPI_abort('Couldt not commit well data MPI structure.')

   end subroutine DefMSWell_mpi_register_mswell_data_description

end module DefMSWell
