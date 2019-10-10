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

   ! FIXME: this could be refactored (several times the same thing...)
   subroutine SyncPetsc_colnum(ColNum)

         integer, dimension(:), intent(out) :: ColNum
         integer, allocatable, dimension(:) :: NbSumCol
         integer :: i, ipc, start, s
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
   
      ! number of node own and frac own in the procs before commRank
         allocate(NbSumCol(Ncpus))
         NbSumCol(:) = 0
         do i=1, Ncpus-1
            NbSumCol(i+1) = NbSumCol(i) &
               + (NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i)) * NbCompThermique &
               + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
         end do
   
         ! ColNum, node
      do i=1, NbNodeLocal
         ipc = NumNodebyProc%Val(i) ! this node is in proc ipc
#ifndef NDEBUG
         if(i<=NbNodeOwn.and.ipc/=commRank) then
            call CommonMPI_abort('inconsistent node proc')
         end if
#endif
         ! NumNodebyProc%Num(i) is the num of this node in the proc that it's own
         do s=1, NbCompThermique
            ColNum((i-1)*NbCompThermique+s) = (NumNodebyProc%Num(i)-1) * NbCompThermique + s + NbSumCol(ipc+1)
            ! if(commRank==0) then
            !    print*, (i-1)*NbCompThermique+s, ColNum((i-1)*NbCompThermique+s), NumNodebyProc%Num(i), ipc, NbSumCol(ipc+1)
            ! end if
         end do
      end do

      ! ColNum, frac
      start = NbNodeLocal * NbCompThermique
      do i=1, NbFracLocal
         ipc = NumFracbyProc%Val(i) ! this frac is in proc ipc
#ifndef NDEBUG
         if(i<=NbFracOwn.and.ipc/=commRank) then
            call CommonMPI_abort('inconsistent frac proc')
         end if
#endif
         ! NumFracbyProc%Num(i) is the num of this frac in the proc that it's own
         do s=1, NbCompThermique
            ColNum(start+(i-1)*NbCompThermique+s) = (NumFracbyProc%Num(i)-1) * NbCompThermique + s &
               + NbSumCol(ipc+1) + NbNodeOwn_Ncpus(ipc+1) * NbCompThermique
         end do
      end do

      ! ColNum, well inj
      start = (NbNodeLocal + NbFracLocal) * NbCompThermique
      do i=1, NbWellInjLocal
         ipc = NumWellInjbyProc%Val(i) ! this well inj is in proc ipc
#ifndef NDEBUG
         if(i<=NbWellInjOwn.and.ipc/=commRank) then
            call CommonMPI_abort('inconsistent injector proc')
         end if
#endif
         ! NumWellInjbyProc%Num(i) is the num of this well inj in the proc that it's own
         ColNum(start+i) = NumWellInjbyProc%Num(i) + NbSumCol(ipc+1) &
            + (NbNodeOwn_Ncpus(ipc+1) + NbFracOwn_Ncpus(ipc+1)) * NbCompThermique
         ! if(commRank==1) then
         !    print*, start+i, ColNum(start+i), ipc, NumWellInjbyProc%Num(i)
         ! end if
      end do

      ! ColNum, well prod
      start = (NbNodeLocal + NbFracLocal) * NbCompThermique + NbWellInjLocal
      do i=1, NbWellProdLocal
         ipc = NumWellProdbyProc%Val(i) ! this well inj is in proc ipc
#ifndef NDEBUG
         if(i<=NbWellProdOwn.and.ipc/=commRank) then
            call CommonMPI_abort('inconsistent producer proc')
         end if
#endif
         ! NumWellInjbyProc%Num(i) is the num of this well inj in the proc that it's own
         ColNum(start+i) = NumWellProdbyProc%Num(i) + NbSumCol(ipc+1) &
            + (NbNodeOwn_Ncpus(ipc+1) + NbFracOwn_Ncpus(ipc+1)) * NbCompThermique + NbWellInjOwn_Ncpus(ipc+1)
      end do

      deallocate(NbSumCol)

   end subroutine SyncPetsc_colnum
   
   ! ! FIXME: this could be refactored (several times the same thing...)
   ! subroutine SyncPetsc_rowcolnum(RowNum, ColNum)
      
   !    integer, dimension(:), intent(out) :: RowNum, ColNum

   !    call SyncPetsc_rownum(RowNum)
   !    call SyncPetsc_colnum(ColNum)

   ! end subroutine SyncPetsc_rowcolnum

   ! subroutine SyncPetsc_rowcolnum_from_C(RowNum, nRowNum, ColNum, nColNum) &
   !    bind(C, name="SyncPetsc_rowcolnum")
   !    type(c_ptr), value, intent(in) :: RowNum, ColNum
   !    integer(c_size_t), value, intent(in) :: nRowNum, nColNum
   !    type(MatrixSize) :: matrix_size
   !    integer(c_int), pointer :: pRowNum(:), pColNum(:)
      
   !    call SyncPetsc_local_matrix_size(matrix_size)
   !    if(matrix_size%nbrows/=nRowNum.or.nRowNum/=nColNum) then
   !       call CommonMPI_abort('Wrong local matrix sizes in SyncPetsc_rowcolnum')
   !    end if
   !    call c_f_pointer(RowNum, pRowNum, (/ nRowNum /))
   !    call c_f_pointer(ColNum, pColNum, (/ nColNum /))
   !    call SyncPetsc_rowcolnum(pRowNum, pColNum)

   ! end subroutine SyncPetsc_rowcolnum_from_C

   subroutine SyncPetsc_colnum_from_C(ColNum, nColNum) &
      bind(C, name="SyncPetsc_colnum")
      type(c_ptr), value, intent(in) :: ColNum
      integer(c_size_t), value, intent(in) :: nColNum
      type(MatrixSize) :: matrix_size
      integer(c_int), pointer :: pColNum(:)
      
      call SyncPetsc_local_matrix_size(matrix_size)
      if(matrix_size%nbrows/=nColNum) then
         call CommonMPI_abort('Wrong local matrix sizes in SyncPetsc_colnum')
      end if
      call c_f_pointer(ColNum, pColNum, (/ nColNum /))
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
