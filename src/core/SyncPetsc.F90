!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module SyncPetsc

   use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_f_pointer, c_size_t, c_long
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort
   use CommonType, only: FamilyDOFIdCOC
   use DefModel, only: NbCompThermique !, psprim, NbContexte
   use MeshSchema, only: &
      NumNodebyProc, NumFracbyProc, NumWellInjbyProc, NumWellProdbyProc, &
      NbWellInjOwn_Ncpus, NbWellInjLocal_Ncpus, NbWellInjLocal_Ncpus, &
      NbWellProdLocal_Ncpus, NbWellProdOwn_Ncpus, &
      NbNodeOwn_Ncpus, NbFracOwn_Ncpus, NbNodeLocal_Ncpus, &
      NbFracLocal_Ncpus, NbCellLocal_Ncpus

   implicit none

   type, bind(C) :: MatrixSize
      integer(c_size_t) :: nbrows, nbcols
   end type MatrixSize

   public:: &
      SyncPetsc_colnum

contains

   subroutine SyncPetsc_global_matrix_size(matrix_size) &
      bind(C, name="SyncPetsc_global_matrix_size")
      type(MatrixSize), intent(out) :: matrix_size
      integer(c_size_t) :: nr, nc
      integer :: i

      nr = 0
      nc = 0
      do i = 1, Ncpus
         nr = nr + (NbNodeLocal_Ncpus(i) + NbFracLocal_Ncpus(i))*NbCompThermique &
              + NbWellInjLocal_Ncpus(i) + NbWellProdLocal_Ncpus(i)
         nc = nc + (NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i))*NbCompThermique &
              + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
      end do
      matrix_size%nbrows = nr
      matrix_size%nbcols = nc

   end subroutine SyncPetsc_global_matrix_size

   subroutine SyncPetsc_local_matrix_size(matrix_size) &
      bind(C, name="SyncPetsc_local_matrix_size")
      type(MatrixSize), intent(out) :: matrix_size

      matrix_size%nbrows = (NbNodeLocal_Ncpus(commRank + 1) &
                            + NbFracLocal_Ncpus(commRank + 1))*NbCompThermique &
                           + NbWellInjLocal_Ncpus(commRank + 1) + NbWellProdLocal_Ncpus(commRank + 1)
      matrix_size%nbcols = (NbNodeOwn_Ncpus(commRank + 1) &
                            + NbFracOwn_Ncpus(commRank + 1))*NbCompThermique &
                           + NbWellInjOwn_Ncpus(commRank + 1) + NbWellProdOwn_Ncpus(commRank + 1)

   end subroutine SyncPetsc_local_matrix_size

   ! ! FIXME: this could be refactored (several times the same thing...)
   ! subroutine SyncPetsc_rownum(RowNum)

   !    integer, dimension(:), intent(out) :: RowNum
   !    integer, allocatable, dimension(:) :: NbSumRow
   !    integer :: i, ipc, start, s
   !    integer :: NbNodeLocal, NbFracLocal, NbWellInjLocal, NbWellProdLocal

   !    NbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)
   !    NbFracLocal = NbFracLocal_Ncpus(commRank + 1)
   !    NbWellInjLocal = NbWellInjLocal_Ncpus(commRank + 1)
   !    NbWellProdLocal = NbWellProdLocal_Ncpus(commRank + 1)

   !    ! number of node and frac in the procs before commRank
   !    allocate(NbSumRow(Ncpus))

   !    NbSumRow(:) = 0
   !    do i=1, Ncpus-1
   !       NbSumRow(i+1) = NbSumRow(i) &
   !          + (NbNodeLocal_Ncpus(i) + NbFracLocal_Ncpus(i)) * NbCompThermique &
   !          + NbWellInjLocal_Ncpus(i) + NbWellProdLocal_Ncpus(i)
   !    end do

   !    ! RowNum, node
   !    do i=1, NbNodeLocal * NbCompThermique
   !       RowNum(i) = i + NbSumRow(commRank+1)
   !    end do

   !    ! RowNum, frac
   !    start = NbNodeLocal * NbCompThermique
   !    do i=1, NbFracLocal * NbCompThermique
   !       RowNum(i+start) = i + start + NbSumRow(commRank+1) ! (node, frac, well inj, well prod)
   !    end do

   !    ! RowNum, well inj
   !    start = (NbNodeLocal + NbFracLocal) * NbCompThermique
   !    do i=1, NbWellInjLocal
   !       RowNum(i+start) = i + start + NbSumRow(commRank+1) ! (node, frac, well inj, well prod)
   !    end do

   !    ! RowNum, well prod
   !    start = (NbNodeLocal + NbFracLocal) * NbCompThermique + NbWellInjLocal
   !    do i=1, NbWellProdLocal
   !       RowNum(i+start) = i + start + NbSumRow(commRank+1) ! (node, frac, well inj, well prod)
   !    end do

   !    deallocate(NbSumRow)

   ! end subroutine SyncPetsc_rownum

   function check_dof_family(dof_family, nb_owns) result(ok)
      type(FamilyDOFIdCOC), intent(in) :: dof_family
      integer, intent(in) :: nb_owns
      logical :: ok

      integer i
      ok = .false.
      if (size(dof_family%ids) < nb_owns) return
      do i = 1, nb_owns
         if (dof_family%ids(i)%proc /= commRank) then
            write (*, *) "inconsistent node proc"
            return
         end if
      end do
      ok = .true.
   end function check_dof_family

   ! FIXME: this could be refactored (several times the same thing...)
   subroutine SyncPetsc_colnum(ColNum)

      integer, dimension(:), intent(out) :: ColNum

      integer :: i, local_col_offset, global_col_offset(NCpus)
      integer :: NbNodeLocal, NbNodeOwn, NbFracLocal, NbFracOwn
      integer :: NbWellInjLocal, NbWellInjOwn, NbWellProdLocal, NbWellProdOwn

      NbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)
      NbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)
      NbFracLocal = NbFracLocal_Ncpus(commRank + 1)
      NbFracOwn = NbFracOwn_Ncpus(commRank + 1)
      NbWellInjLocal = NbWellInjLocal_Ncpus(commRank + 1)
      NbWellInjOwn = NbWellInjOwn_Ncpus(commRank + 1)
      NbWellProdLocal = NbWellProdLocal_Ncpus(commRank + 1)
      NbWellProdOwn = NbWellProdOwn_Ncpus(commRank + 1)

      local_col_offset = 0
      global_col_offset = 0
      if (Ncpus > 1) then
         do i = 1, Ncpus - 1
            global_col_offset(i + 1) = global_col_offset(i) &
                                       + (NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i))*NbCompThermique &
                                       + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
         end do
      end if

#ifndef NDEBUG
      if (size(NumNodebyProc%ids) /= NbNodeLocal) &
         call CommonMPI_abort("Inconsistent number of node dofs")
      if (.not. check_dof_family(NumNodebyProc, NbNodeOwn)) &
         call CommonMPI_abort("Inconsistent node family")
#endif
      call SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, NumNodebyProc, NbCompThermique, ColNum)
      local_col_offset = local_col_offset + NbNodeLocal*NbCompThermique
      global_col_offset = global_col_offset + NbCompThermique*NbNodeOwn_Ncpus

      if (NbFracLocal > 0) then
#ifndef NDEBUG
         if (size(NumFracbyProc%ids) /= NbFracLocal) &
            call CommonMPI_abort("Inconsistent number of Frac dofs")
         if (.not. check_dof_family(NumFracbyProc, NbFracOwn)) &
            call CommonMPI_abort("Inconsistent Frac family")
#endif
         call SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, NumFracbyProc, NbCompThermique, ColNum)
         local_col_offset = local_col_offset + NbFracLocal*NbCompThermique
      end if
      global_col_offset = global_col_offset + NbCompThermique*NbFracOwn_Ncpus

      if (NbWellInjLocal > 0) then
#ifndef NDEBUG
         if (size(NumWellInjbyProc%ids) /= NbWellInjLocal) &
            call CommonMPI_abort("Inconsistent number of WellInj dofs")
         if (.not. check_dof_family(NumWellInjbyProc, NbWellInjOwn)) &
            call CommonMPI_abort("Inconsistent WellInj family")
#endif
         call SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, NumWellInjbyProc, 1, ColNum)
         local_col_offset = local_col_offset + 1*NbWellInjLocal
      endif
      global_col_offset = global_col_offset + 1*NbWellInjOwn_Ncpus

      if (NbWellProdLocal > 0) then
#ifndef NDEBUG
         if (size(NumWellProdbyProc%ids) /= NbWellProdLocal) &
            call CommonMPI_abort("Inconsistent number of WellProd dofs")
         if (.not. check_dof_family(NumWellProdbyProc, NbWellProdOwn)) &
            call CommonMPI_abort("Inconsistent WellProd family")
#endif
         call SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, NumWellProdbyProc, 1, ColNum)
         ! local_col_offset = local_col_offset + 1 * NbWellProdLocal
      end if
      ! global_col_offset = global_col_offset + 1 * NbWellProdOwn_Ncpus

   end subroutine SyncPetsc_colnum

   subroutine SyncPetsc_colnum_dof_family(local_col_offset, global_col_offset, dof_family, dof_size, col_map)
      integer, intent(in) :: local_col_offset ! offset in terms of local columns
      integer, dimension(Ncpus), intent(in) :: global_col_offset
      type(FamilyDOFIdCOC), intent(in) :: dof_family
      integer, intent(in) :: dof_size
      integer, dimension(:), intent(inout) :: col_map

      integer :: i, s, n, proc, local_id

      n = size(dof_family%ids)
      do i = 1, n
         proc = dof_family%ids(i)%proc
         local_id = dof_family%ids(i)%local_id
         do s = 1, dof_size
            col_map(local_col_offset + (i - 1)*dof_size + s) = global_col_offset(proc + 1) + (local_id - 1)*dof_size + s
         end do
      end do

   end subroutine SyncPetsc_colnum_dof_family

   subroutine SyncPetsc_colnum_from_C(ColNum, nColNum) &
      bind(C, name="SyncPetsc_colnum")
      type(c_ptr), value, intent(in) :: ColNum
      integer(c_size_t), value, intent(in) :: nColNum
      type(MatrixSize) :: matrix_size
      integer(c_int), pointer :: pColNum(:)

      call SyncPetsc_local_matrix_size(matrix_size)
      if (matrix_size%nbrows /= nColNum) then
         call CommonMPI_abort('Wrong local matrix sizes in SyncPetsc_colnum')
      end if
      call c_f_pointer(ColNum, pColNum, (/nColNum/))
      call SyncPetsc_colnum(pColNum)

   end subroutine SyncPetsc_colnum_from_C

end module SyncPetsc

subroutine syncpetsc_getsolnodefracwell(x_s, increments_pointers)

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
      NbNodeLocal_Ncpus, NbFracLocal_Ncpus

   implicit none

   Vec, intent(inout) :: x_s
   type(Newton_increments_pointers), intent(in), value :: increments_pointers
   type(Newton_increments) :: increments
   integer :: i, j, start
   integer :: NbNodeLocal, NbFracLocal
   PetscErrorCode :: Ierr
   double precision, pointer :: ptr(:)

   NbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)
   NbFracLocal = NbFracLocal_Ncpus(commRank + 1)

   call Newton_pointers_to_values(increments_pointers, increments)

   ! get values from x_s
   call VecGetArrayF90(x_s, ptr, Ierr)
   CHKERRQ(Ierr)

   ! increment node
   do i = 1, NbNodeLocal
      start = (i - 1)*NbCompThermique
      do j = 1, NbCompThermique
         increments%nodes(j, i) = ptr(start + j)
      end do
   end do

   ! increment frac
   do i = 1, NbFracLocal
      start = (i + NbNodeLocal - 1)*NbCompThermique
      do j = 1, NbCompThermique
         increments%fractures(j, i) = ptr(start + j)
      end do
   end do

   ! increment well inj
   start = (NbNodeLocal + NbFracLocal)*NbCompThermique
   do i = 1, NbWellInjLocal_Ncpus(commRank + 1)
      increments%injectors(i) = ptr(start + i)
   end do

   ! increment well prod
   start = start + NbWellInjLocal_Ncpus(commRank + 1)
   do i = 1, NbWellProdLocal_Ncpus(commRank + 1)
      increments%producers(i) = ptr(start + i)
   end do

   call VecRestoreArrayF90(x_s, ptr, Ierr)
   CHKERRQ(Ierr)

end subroutine syncpetsc_getsolnodefracwell
