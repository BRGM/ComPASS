! This file is to be included in fortran files that need to check petsc error
! As of PETSc 3.17 PetscCall supersedes CHKERRQ
! (cf. https://petsc.org/release/docs/changes/317)
! Yet both are just the same macro definition so for the time being
! we just define our own macro defintion

#include <petsc/finclude/petscsys.h>
#ifndef CMP_PETSC_CHECK
#if defined(PETSC_HAVE_FORTRAN_FREE_LINE_LENGTH_NONE)
#define CMP_PETSC_CHECK(ierr)      \
   if (ierr /= 0) then;            \
   call PetscErrorF(ierr, 0, "X"); \
   return;                         \
   endif
#else
#define CMP_PETSC_CHECK(ierr) \
   if (ierr /= 0) then;       \
   call PetscErrorF(ierr);    \
   return;                    \
   endif
#endif
#endif
