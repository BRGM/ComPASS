! This is based on PETSc example #7:
! https://www.mcs.anl.gov/petsc/petsc-current/src/vec/vec/examples/tutorials/ex7.c
! https://www.mcs.anl.gov/petsc/petsc-current/src/vec/vec/examples/tutorials/ex7f.F
! Passing PETSc objects between C and Fortran does not rely on iso C binding
! but on opaque pointers and compiler dependant name mangling hence some restrictions
! on the case of the subroutine names and the fact that subroutines cannot
! be encapsulated in a Fortran module.

! WARNING: subroutine name must be lower case
subroutine dump_from_fortran(mat)

#include <petsc/finclude/petsc.h>

   use petsc

   implicit none

   Mat, intent(in) :: mat
   PetscErrorCode :: ierr

   call MatView(mat, PETSC_VIEWER_STDOUT_WORLD, ierr)

end subroutine dump_from_Fortran
