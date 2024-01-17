!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

#include <ComPASS_PETSc_definitions.h>

module SyncPetscMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_f_pointer, c_size_t, c_long
   use CommonMPI
   use CommonType
   use DefModel
   use MeshSchema
   use SyncPetsc
#else
   use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_f_pointer, c_size_t, c_long
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort
   use CommonType, only: FamilyDOFIdCOC
   use DefModel, only: NbCompThermique !, psprim, NbContexte
   use MeshSchema, only: &
      NbMSWellNodeLocal_Ncpus, NbMSWellNodeOwn_Ncpus, NumMSWellNodebyProc
   use SyncPetsc, only: &
      MatrixSize, SyncPetsc_colnum_dof_family, check_dof_family
#endif
   implicit none

   public:: &
      SyncPetscMSWells_colnum

contains

   subroutine SyncPetscMSWells_global_matrix_size(matrix_size) &
      bind(C, name="SyncPetscMSWells_global_matrix_size")
      type(MatrixSize), intent(out) :: matrix_size
      integer(c_size_t) :: nr, nc
      integer :: i

   end subroutine SyncPetscMSWells_global_matrix_size

   subroutine SyncPetscMSWells_local_matrix_size(matrix_size) &
      bind(C, name="SyncPetscMSWells_local_matrix_size")
      type(MatrixSize), intent(out) :: matrix_size

      matrix_size%nbrows = (NbMSWellNodeLocal_Ncpus(commRank + 1))*NbCompThermique

      matrix_size%nbcols = (NbMSWellNodeOwn_Ncpus(commRank + 1))*NbCompThermique

   end subroutine SyncPetscMSWells_local_matrix_size

   ! FIXME: this could be refactored (several times the same thing...)
   subroutine SyncPetscMSWells_colnum(ColNum)
      integer, dimension(:), intent(out) :: ColNum

      integer :: i, local_col_offset, global_col_offset(NCpus)
      integer :: NbMSWellNodeLocal, NbMSWellNodeOwn

      NbMSWellNodeLocal = NbMSWellNodeLocal_Ncpus(commRank + 1)
      NbMSWellNodeOwn = NbMSWellNodeOwn_Ncpus(commRank + 1)

      local_col_offset = 0
      global_col_offset = 0
      if (Ncpus > 1) then
         do i = 1, Ncpus - 1
            global_col_offset(i + 1) = global_col_offset(i) &
                                       + (NbMSWellNodeOwn_Ncpus(i))*NbCompThermique
         end do
      end if

#ifndef NDEBUG
      if (size(NumMSWellNodebyProc%ids) /= NbMSWellNodeLocal) &
         call CommonMPI_abort("Inconsistent number of mswell nodes dofs in SyncPetsMSWells")
      if (.not. check_dof_family(NumMSWellNodebyProc, NbMSWellNodeOwn)) &
         call CommonMPI_abort("Inconsistent mswell nodes family in SyncPetsMSWells")
#endif
      call SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, NumMSWellNodebyProc, NbCompThermique, ColNum)
      local_col_offset = local_col_offset + NbMSWellNodeLocal*NbCompThermique
      global_col_offset = global_col_offset + NbCompThermique*NbMSWellNodeOwn_Ncpus

   end subroutine SyncPetscMSWells_colnum

   subroutine SyncPetscMSWells_colnum_dof_family(local_col_offset, global_col_offset, dof_family, dof_size, col_map)
      integer, intent(in) :: local_col_offset ! offset in terms of local columns
      integer, dimension(Ncpus), intent(in) :: global_col_offset
      type(FamilyDOFIdCOC), intent(in) :: dof_family
      integer, intent(in) :: dof_size
      integer, dimension(:), intent(inout) :: col_map

   end subroutine SyncPetscMSWells_colnum_dof_family

   subroutine SyncPetscMSWells_colnum_from_C(ColNum, nColNum) &
      bind(C, name="SyncPetscMSWells_colnum")
      type(c_ptr), value, intent(in) :: ColNum
      integer(c_size_t), value, intent(in) :: nColNum
      type(MatrixSize) :: matrix_size
      integer(c_int), pointer :: pColNum(:)

      call SyncPetscMSWells_local_matrix_size(matrix_size)
      if (matrix_size%nbrows /= nColNum) then
         call CommonMPI_abort('Wrong local matrix sizes in SyncPetsc_colnum')
      end if
      call c_f_pointer(ColNum, pColNum, (/nColNum/))
      call SyncPetscMSWells_colnum(pColNum)

   end subroutine SyncPetscMSWells_colnum_from_C

end module SyncPetscMSWells

subroutine syncpetscmswells_getsolmswell(x_s, increments_pointers)

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petsc.h>
#endif

   use petsc

   use CommonMPI, only: commRank
   use Newton, only: Newton_increments, Newton_pointers_to_values, Newton_increments_pointers
   use DefModel, only: NbCompThermique !, psprim, NbContexte
   use MeshSchema, only: &
      NbWellInjLocal_Ncpus, NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus, &
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus, &
      NbMSWellNodeLocal_Ncpus, NbMSWellNodeOwn_Ncpus

   implicit none

   !The Newton_increments vectors contain the solution of the whole complete system, i.e.,
   !nodes, fracs, well_prod, well_inj and mswell_nodes. But x_s contains only the mswell_nodes sol

   Vec, intent(inout) :: x_s
   type(Newton_increments_pointers), intent(in), value :: increments_pointers
   type(Newton_increments) :: increments
   integer :: i, j, start
   integer :: NbMSWellNodeLocal
   PetscErrorCode :: Ierr
   double precision, pointer :: ptr(:)

   NbMSWellNodeLocal = NbMSWellNodeLocal_Ncpus(commRank + 1)

   call Newton_pointers_to_values(increments_pointers, increments)

   ! get values from x_s
   call VecGetArrayF90(x_s, ptr, Ierr)
   CMP_PETSC_CHECK(Ierr)

   ! increment mswell_nodes
   do i = 1, NbMSWellNodeLocal
      start = (i - 1)*NbCompThermique
      do j = 1, NbCompThermique
         increments%mswell_nodes(j, i) = ptr(start + j)
      end do
   end do

   call VecRestoreArrayF90(x_s, ptr, Ierr)
   CMP_PETSC_CHECK(Ierr)

end subroutine syncpetscmswells_getsolmswell
