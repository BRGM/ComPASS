!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

! Note: well is only consided in cpramg with pressure
!       i.e. SolvePetsc_cpramgPCApply_P_multiplicative
!       well is not added in SolvePetsc_At

module SolvePetsc

  use iso_c_binding, only: c_bool, c_int, c_double, c_ptr, c_f_pointer

  use CommonMPI, only: commRank, ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort
  use CommonType, only: FamilyDOFIdCOC
  use DefModel, only: psprim, NbCompThermique, NbContexte
  use IncCVReservoir, only: IncNode, IncFrac
  use MeshSchema, only: &
     NumNodebyProc, NumFracbyProc, NumWellInjbyProc, NumWellProdbyProc, &
     NbWellInjOwn_Ncpus, NbWellInjLocal_Ncpus, NbWellInjLocal_Ncpus, &
     NbWellProdLocal_Ncpus, NbWellProdOwn_Ncpus, &
     NbNodeOwn_Ncpus, NbFracOwn_Ncpus, NbNodeLocal_Ncpus, &
#ifdef _WIP_FREEFLOW_STRUCTURES_
     IdFFNodeLocal, &
#endif
     NbFracLocal_Ncpus, NbCellLocal_Ncpus
  use Jacobian, only: JacA, Sm

  ! tmp
  use IncPrimSecd, only: NbIncTotalPrim_ctx 


#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h> 
#else
#include <petsc/finclude/petsc.h>
#endif

  use petsc

  implicit none

  ! Solver: A_mpi x_mpi = Sm_mpi
  Mat, private :: A_mpi
  Vec, private :: Sm_mpi
!   Vec, private :: x_mpi
  Vec, private :: y_mpi

  KSP, private :: ksp_mpi ! ksp
  PC,  private :: pc_mpi  ! pc

  ! Preconditioner CPR-AMG
  Mat, private :: Ap       ! Pression part of A_mpi
  PC,  private :: pcamg_p  ! pc amg for pression
  PC,  private :: pcilu0   ! pc ilu0 for all
  Vec, private :: v_pt, P1v_pt, P1v, AP1v ! tmp vectors

#ifdef _THERMIQUE_
  Mat, private :: At       ! Temperature part of A_mpi
  PC,  private :: pcamg_t  ! pc amg for temperature
  Vec, private :: AP1v_t, P2AP1v, AP2AP1v, P2AP1v_t

  logical, allocatable, dimension(:) :: &
       IsTprimNodeFracOwn, & ! T is prim/secd, element in (node own, frac own)
       IsTprimNodeFracLocal  ! T is prim/secd, element in (node local, frac local)
#endif

  integer, private :: kspitmax          ! max nb of iterations
  double precision, private :: PetscKspTol   ! tolerance

  ! vector contains history of convergence
  double precision, allocatable, dimension(:), private :: &
       kspHistory

  ! size of A_mpi
  integer, private :: &
       NBlockrowL, & ! number of local (block) rows
       NBlockcolL, & ! number of local (block) cols
       NBlockrowG, & ! number of global (block) rows
       NBlockcolG    ! number of global (block) cols

  integer, private :: &
       NrowL, & ! number of local (point) rows
       NcolL, & ! number of local (point) cols
       NrowG, & ! number of global (point) rows
       NcolG    ! number of global (point) cols


  ! Blockrowstart(i): sum of block rows (node own + frac own + well own) before proc i
  ! Blockcolstart(i): sum of block cols (node own + frac own + well own) before proc i = Blockrowstart
  integer, allocatable, dimension(:), private :: &
       Blockrowstart, &
       Blockcolstart

  ! rowstart(i): sum of rows before proc i: (node own+frac own) * NbCompThermique + well own
  ! colstart(i): sum of cols before proc i: (node own+frac own) * NbCompThermique + well own = rowstart
  integer, allocatable, dimension(:), private :: &
       rowstart, &
       colstart

  integer, allocatable, dimension(:), private :: &
       RowLToRowGBlock, &
       ColLToColGBlock, &
       RowLToRowG, &
       ColLToColG

  ! tmp values
  integer, private ::   &
       NbNodeOwn,       &
       NbNodeLocal,     &
       NbFracOwn,       &
       NbFracLocal,     &
       NbCellLocal,     &
       NbWellInjOwn,    &
       NbWellProdOwn,   &
       NbWellInjLocal,  &
       NbWellProdLocal

  public :: &
       SolvePetsc_Init,               &
       SolvePetsc_SetUp,              &
       SolvePetsc_KspSolve,           &
       SolvePetsc_check_solution,     & 
       SolvePetsc_free,               &
       SolvePetsc_cpramgFree,         &
       SolvePetsc_Ksp_configuration,  &
       SolvePetsc_KspSolveIterations

  private :: &
       SolvePetsc_Init_cpramg_specific,  &
       SolvePetsc_CreateKsp,  &
       SolvePetsc_CreateAmpi, &
       SolvePetsc_CreateSm,   &
       SolvePetsc_SetAmpi,    &
       SolvePetsc_SetSm,      &
                                !
       SolvePetsc_cpramgCreateApAt,  &
       SolvePetsc_cpramgCreateKsp,   &
       SolvePetsc_cpramgPCSetUp,     &
                                !
       SolvePetsc_cpramgPCApply_P_multiplicative,  &
       SolvePetsc_cpramgPCApply_P_additive,  &
                                !
       SolvePetsc_SetAp,       &
       SolvePetsc_LtoG,        &
       SolvePetsc_RowColStart

#ifdef _THERMIQUE_

  private :: &
       SolvePetsc_SetAt, &
       SolvePetsc_cpramgPCApply_PT_multiplicative, &
       SolvePetsc_cpramgPCApply_PT_additive, &
       SolvePetsc_cpramgPCApply_T_multiplicative
#endif

contains

!< Create structure of mat and solver
  subroutine SolvePetsc_Init(kspitmax_in, ksptol_in, &
                             activate_cpramg, activate_direct_solver)

    integer, intent(in) :: kspitmax_in
    double precision, intent(in) :: ksptol_in
    logical(c_bool), intent(in) :: activate_direct_solver, activate_cpramg

    integer :: i

    kspitmax = kspitmax_in
    PetscKspTol = ksptol_in

    ! set tmp values
    NbNodeOwn = NbNodeOwn_Ncpus(commRank+1)
    NbNodeLocal = NbNodeLocal_Ncpus(commRank+1)
    NbFracOwn = NbFracOwn_Ncpus(commRank+1)
    NbFracLocal = NbFracLocal_Ncpus(commRank+1)
    NbCellLocal = NbCellLocal_Ncpus(commRank+1)

    NbWellInjOwn = NbWellInjOwn_Ncpus(commRank+1)
    NbWellInjLocal = NbWellInjLocal_Ncpus(commRank+1)
    NbWellProdOwn = NbWellProdOwn_Ncpus(commRank+1)
    NbWellProdLocal = NbWellProdLocal_Ncpus(commRank+1)

    ! local row/col size:  node own and frac own and well own
    NBlockrowL = NbNodeOwn + NbFracOwn + NbWellInjOwn + NbWellProdOwn ! =JacA%Nb
    NBlockcolL = NBlockrowL

    ! global row/col size: sum of all procs
    NBlockrowG = 0
    do i=1, Ncpus
       NBlockrowG = NBlockrowG + NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i) &
            + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
    end do
    NBlockcolG = NBlockrowG

    call SolvePetsc_RowColStart
    call SolvePetsc_CreateAmpi
    call SolvePetsc_CreateSm
    
    if(activate_cpramg.and..not.activate_direct_solver) then
        call SolvePetsc_Init_cpramg_specific(kspitmax_in, ksptol_in)
    else
        if(activate_direct_solver) then
            call SolvePetsc_CreateKsp_direct_solver
        else 
            call SolvePetsc_CreateKsp
        endif
    endif

    ! compute RowLToRowG and ColLToColG
    call SolvePetsc_LtoG
    call SolvePetsc_LtoGBlock

    ! delete row/col start
    deallocate(rowstart)
    deallocate(colstart)
    deallocate(Blockrowstart)
    deallocate(Blockcolstart)

  end subroutine SolvePetsc_Init

  subroutine SolvePetsc_Init_cpramg_specific(kspitmax_in, ksptol_in)

    integer, intent(in) :: kspitmax_in
    double precision, intent(in) :: ksptol_in

    integer :: i
    PetscErrorCode :: Ierr

    call SolvePetsc_cpramgCreateApAt
    call SolvePetsc_cpramgCreateKsp

    ! create tmp vectors for CPR-AMG
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NBlockrowL, NBlockrowG, v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecSet(v_pt, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(v_pt, Ierr)
    CHKERRQ(Ierr)

    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NBlockrowL, NBlockrowG, P1v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecSet(P1v_pt, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(P1v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(P1v_pt, Ierr)
    CHKERRQ(Ierr)

    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NrowL, NrowG, P1v, Ierr)
    CHKERRQ(Ierr)
    call VecSet(P1v, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(P1v, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(P1v, Ierr)
    CHKERRQ(Ierr)

    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NrowL, NrowG, AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecSet(AP1v, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(AP1v, Ierr)
    CHKERRQ(Ierr)

#ifdef _THERMIQUE_

    ! AP1v_t
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NBlockrowL, NBlockrowG, AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecSet(AP1v_t, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(AP1v_t, Ierr)
    CHKERRQ(Ierr)

    ! P2AP1v_t
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NBlockrowL, NBlockrowG, P2AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecSet(P2AP1v_t, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(P2AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(P2AP1v_t, Ierr)
    CHKERRQ(Ierr)

    ! P2AP1v
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NrowL, NrowG, P2AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecSet(P2AP1v, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(P2AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(P2AP1v, Ierr)
    CHKERRQ(Ierr)

    ! P2AP1v
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NrowL, NrowG, AP2AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecSet(AP2AP1v, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(AP2AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(AP2AP1v, Ierr)
    CHKERRQ(Ierr)

    allocate(IsTprimNodeFracOwn(NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)))
    allocate(IsTprimNodeFracLocal(NbNodeLocal_Ncpus(commRank+1)+NbFracLocal_Ncpus(commRank+1)))

#endif

  end subroutine SolvePetsc_Init_cpramg_specific


  ! compute Blockrowstart and Blockcolstart
  ! compute rowstart and colstart
  subroutine SolvePetsc_RowColStart

    integer :: i

    ! Blockrowstart(i) and Blockcolstart(i):
    !     sum of node own and frac own and well own before proc i,
    !     the block element(node/frac) is considered as 1 row/col
    ! the rows between Blockrowstart(i)+1 and Blockrowstart(i+1)
    !     are in proc i-1, where i=1,2,...,Ncpus
    allocate(Blockrowstart(Ncpus))
    allocate(Blockcolstart(Ncpus))

    Blockrowstart(1) = 0
    do i=1, Ncpus-1 !
       Blockrowstart(i+1) = Blockrowstart(i) &
            + NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i)  &
            + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
    end do
    Blockcolstart(:) = Blockrowstart(:)

    ! rowstart(i) and colstart(i):
    !    sum of node own and frac own and well own before proc i
    !    the block element (node/frac) is considered as NbCompthermique row/col
    ! the rows between Blockrowstart(i)+1 and Blockrowstart(i+1)
    !    are in proc i-1, where i=1,2,...,Ncpus
    allocate(rowstart(Ncpus))
    allocate(colstart(Ncpus))

    rowstart(1) = 0
    do i=1, Ncpus-1 !
       rowstart(i+1) = rowstart(i) &
            + (NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i)) * NbCompThermique  &
            + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
    end do
    colstart(:) = rowstart(:)

  end subroutine SolvePetsc_RowColStart


  ! Setup: set values mat and second member
  subroutine SolvePetsc_SetUp() &
    bind(C, name="SolvePetsc_SetUp")

    PetscErrorCode :: Ierr

    ! set values: Ampi and Sm
    call SolvePetsc_SetAmpi
    call SolvePetsc_SetSm

    call PCSetUp(pc_mpi, Ierr)
    CHKERRQ(Ierr)
    call KSPSetUp(ksp_mpi, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_SetUp


  ! do setup required for shell pc cpr-amg
  subroutine SolvePetsc_cpramgPCSetUp(pcin, Ierr)

    PC, intent(inout) :: pcin
    PetscErrorCode, intent(out) :: Ierr

    ! set values: Ap, At
    call SolvePetsc_SetAp

! #ifdef _THERMIQUE_
!     call SolvePetsc_SetAt
! #endif

    call PCSetUp(pcamg_p, Ierr)
    CHKERRQ(Ierr)

    call PCSetUp(pcilu0, Ierr)
    CHKERRQ(Ierr)

! #ifdef _THERMIQUE_

!     call PCSetUp(pcamg_t, Ierr)
!     CHKERRQ(Ierr)
! #endif

  end subroutine SolvePetsc_cpramgPCSetUp


  ! Create Mat A_mpi in AIJ format (the default parallel PETSc format)
  !
  ! A_mpi for 3 mpi procs example:
  !     node own/frac own(mpi=0), (mpi=1),  (mpi=2)
  !         |   a11,                a12,     a13    | mpi 0, node own/frac own
  ! A_mpi = |   a21,                a22,     a22    | mpi 1
  !         |   a31,                a32,     a33    | mpi 2
  !
  ! In mpi 0, insert its JacA to (a11, a12, a13)
  ! In mpi 1,  ...            to (a21, a22, a23)
  !    ...
  ! Rq important:
  !    the rows of JacA is composed with node own and frac own
  !    the cols of JacA is composed with node local and frac local
  subroutine SolvePetsc_CreateAmpi

    integer, allocatable, dimension(:) :: &
         d_nnz, & ! number of nonzeros per row in diagonal portion
         o_nnz    ! number of nonzeros per row in the off-diagonal portion

    ! some tmp values, explained when setting values
    integer :: Nl_Fo, Nl_Fl, Nl_Fl_WIo, Nl_Fl_WIl, Nl_Fl_WIl_WPo, Nl_Fl_WIl_WPl


    ! tmp values
    integer :: i, li, j, s, lj, start, Ierr, colend

    ! local row/col size:  node own and frac own
    NrowL = (NbNodeOwn + NbFracOwn) * NbCompThermique &
         + NbWellInjOwn + NbWellProdOwn
    NcolL = NrowL

    ! global row/col size: sum of all procs
    NrowG = 0
    do i=1, Ncpus
       NrowG = NrowG &
            + (NbNodeOwn_Ncpus(i) + NbFracOwn_Ncpus(i)) * NbCompThermique &
            + NbWellInjOwn_Ncpus(i) + NbWellProdOwn_Ncpus(i)
    end do
    NcolG = NrowG

    ! number of nonzeros per row in diag or off-diag portion
    allocate(d_nnz(NrowL))
    allocate(o_nnz(NrowL))
    d_nnz(:) = 0
    o_nnz(:) = 0

    ! the diagonal portion is between colstart and colend
    colend = colstart(commRank+1) + NcolL

    ! N: Node, F: Frac, WI: well inj, WP: well prod
    ! l: local, o: own
    Nl_Fo = NbNodeLocal + NbFracOwn
    Nl_Fl = NbNodeLocal + NbFracLocal
    Nl_Fl_WIo = NbNodeLocal + NbFracLocal + NbWellInjOwn
    Nl_Fl_WIl = NbNodeLocal + NbFracLocal + NbWellInjLocal
    Nl_Fl_WIl_WPo = NbNodeLocal + NbFracLocal + NbWellInjLocal + NbWellProdOwn
    Nl_Fl_WIl_WPl = NbNodeLocal + NbFracLocal + NbWellInjLocal + NbWellProdLocal

    ! rows associated with node own and frac own
    do i=1, NbNodeOwn + NbFracOwn

       do j=JacA%Pt(i)+1, JacA%Pt(i+1)

          lj = JacA%Num(j)

          if (lj<=NbNodeOwn) then ! node own
             do s=1, NbCompThermique
                d_nnz((i-1)*(NbCompThermique)+s) = d_nnz((i-1)*(NbCompThermique)+s) + NbCompThermique
             end do

          else if (NbNodeOwn<lj .and. lj<=NbNodeLocal) then ! node ghost
             do s=1, NbCompThermique
                o_nnz((i-1)*(NbCompThermique)+s) = o_nnz((i-1)*(NbCompThermique)+s) + NbCompThermique
             end do

          else if (NbNodeLocal<lj .and. lj<=Nl_Fo) then
             do s=1, NbCompThermique ! frac own
                d_nnz((i-1)*(NbCompThermique)+s) = d_nnz((i-1)*(NbCompThermique)+s) + NbCompThermique
             end do

          else if(Nl_Fo<lj .and. lj<= Nl_Fl) then
             do s=1, NbCompThermique ! frac ghost
                o_nnz((i-1)*(NbCompThermique)+s) = o_nnz((i-1)*(NbCompThermique)+s) + NbCompThermique
             end do

          else if(Nl_Fl<lj .and. lj<= Nl_Fl_WIo) then
             do s=1, NbCompThermique ! well inj own
                d_nnz((i-1)*(NbCompThermique)+s) = d_nnz((i-1)*(NbCompThermique)+s) + 1
             end do

          else if(Nl_Fl_WIo<lj .and. lj<= Nl_Fl_WIl) then
             do s=1, NbCompThermique ! well inj ghost
                o_nnz((i-1)*(NbCompThermique)+s) = o_nnz((i-1)*(NbCompThermique)+s) + 1
             end do

          else if(Nl_Fl_WIl<lj .and. lj<= Nl_Fl_WIl_WPo) then
             do s=1, NbCompThermique ! well prod own
                d_nnz((i-1)*(NbCompThermique)+s) = d_nnz((i-1)*(NbCompThermique)+s) + 1
             end do
          else
             do s=1, NbCompThermique ! well prod ghost
                o_nnz((i-1)*(NbCompThermique)+s) = o_nnz((i-1)*(NbCompThermique)+s) + 1
             end do
          end if
       end do
    end do

    ! rows associated with well inj own and well prod own
    start = (NbNodeOwn+NbFracOwn) * NbCompThermique

    do i=1, NbWellInjOwn+NbWellProdOwn

       li = i + NbNodeOwn + NbFracOwn

       do j=JacA%Pt(li)+1, JacA%Pt(li+1)

          lj = JacA%Num(j)

          if (lj<=NbNodeOwn) then ! node own
             d_nnz(start+i) = d_nnz(start+i) + NbCompThermique

          else if (NbNodeOwn<lj .and. lj<=NbNodeLocal) then ! node ghost
             o_nnz(start+i) = o_nnz(start+i) + NbCompThermique

          else if (NbNodeLocal<lj .and. lj<=Nl_Fo) then ! frac own
             d_nnz(start+i) = d_nnz(start+i) + NbCompThermique

          else if(Nl_Fo<lj .and. lj<= Nl_Fl) then ! frac ghost
             o_nnz(start+i) = o_nnz(start+i) + NbCompThermique

          else if(Nl_Fl<lj .and. lj<= Nl_Fl_WIo) then ! well inj own
             d_nnz(start+i) = d_nnz(start+i) + 1

          else if(Nl_Fl_WIo<lj .and. lj<= Nl_Fl_WIl) then ! well inj ghost
             o_nnz(start+i) = o_nnz(start+i) + 1

          else if(Nl_Fl_WIl<lj .and. lj<= Nl_Fl_WIl_WPo) then ! well prod own
             d_nnz(start+i) = d_nnz(start+i) + 1

          else if(Nl_Fl_WIl_WPo<lj .and. lj<=Nl_Fl_WIl_WPl) then ! well prod ghost
             o_nnz(start+i) = o_nnz(start+i) + 1
          else
             print*, "error in create ampi"
          end if
       end do
    end do

    !write(*,*) 'proc', commRank, 'has', &
    !            NbNodeOwn, 'nodes own', &
    !            NbFracOwn, 'fractures own', &
    !            NbCompThermique, 'components + thermal'
    !write(*,*) 'create sparse matrix on proc', commRank, &
    !           'with structure', &
    !            NrowL, '(out of', NrowG, ') x', &
    !            NcolL, '(out of', NcolG, ')'
    call MatCreateAIJ(ComPASS_COMM_WORLD, &
         NrowL, NcolL, &
         NrowG, NcolG, &
         0, d_nnz, 0, o_nnz, A_mpi, Ierr)
    CHKERRQ(Ierr)

    deallocate(d_nnz)
    deallocate(o_nnz)

  end subroutine SolvePetsc_CreateAmpi


  ! Set values A_mpi
  subroutine SolvePetsc_SetAmpi

    integer :: i, j, s
    integer :: row, col
    PetscErrorCode :: Ierr

    integer :: &
         m, idxm(NbCompThermique), & ! number of rows and their global indices
         n, idxn(NbCompThermique)    ! number of cols and their global indices

    m = NbCompThermique
    n = NbCompThermique

    ! cols of JacA are (node, frac, wellinj, wellprod)
    ! insert JacA(i,j), i is node own or frac own or wellinj own or wellprod own;
    !                   j is node or frac or wellinj or wellprod
    ! insert value JacA%Val(:,:,*), the index of "*" in A_mpi is (row, col)

    ! rows associated with node own and frac own
    do i=1, NbNodeOwn + NbFracOwn

       do j=JacA%Pt(i)+1, JacA%Pt(i+1)

          row = RowLToRowG(i) - 1             ! 0-based in petsc
          col = ColLToColG( JacA%Num(j) ) - 1 ! 0-based in petsc

          ! if(commRank==0 .and. i==10) then
          !    print*, i, j, row+1, col+1, JacA%Num(j)
          ! end if

          ! col is node or frac, insert JacA%Val(:,:,j)
          if(JacA%Num(j) <= (NbNodeLocal+NbFracLocal)) then

             do s=1, NbCompThermique
                idxm(s) = row + s - 1
                idxn(s) = col + s - 1
             end do

             call MatSetValues(A_mpi, m, idxm, n, idxn, JacA%Val(:,:,j), INSERT_VALUES, Ierr)
             CHKERRQ(Ierr)

          else ! col is wellinj or wellprod, insert JacA%Val(1,:,j)

             do s=1, NbCompThermique
                idxm(s) = row + s - 1
             end do
             idxn(1) = col

             ! index order of matrix JacBigA(:,:,j) is (col,row)
            !write(*,*) 'Setting values on', commRank
             call MatSetValues(A_mpi, m, idxm, 1, idxn, JacA%Val(1,:,j), INSERT_VALUES, Ierr)
            !write(*,*) 'Setting values done on', commRank
             CHKERRQ(Ierr)
          end if
       end do
    end do

    ! rows associated with  wellinj own and wellprod own
    do i=(NbNodeOwn+NbFracOwn)+1, (NbNodeOwn+NbFracOwn+NbWellInjOwn+NbWellProdOwn)

       do j=JacA%Pt(i)+1, JacA%Pt(i+1)
          !write(*,*) 'Well jacobian at local index', j, 'size', size(JacA%Val, 1), 'x', size(JacA%Val, 2)
          !write(*,*) JacA%Val(:,:,j)
          row = RowLToRowG(i) - 1             ! 0-based in petsc
          col = ColLToColG( JacA%Num(j) ) - 1 ! 0-based in petsc

          ! col is node or frac, insert JacA%Val(:,1,j)
          if(JacA%Num(j) <= (NbNodeLocal+NbFracLocal)) then

             do s=1, NbCompThermique
                idxn(s) = col + s - 1
             end do
             idxm(1) = row

             call MatSetValues(A_mpi, 1, idxm, n, idxn, JacA%Val(:,1,j), INSERT_VALUES, Ierr)
             CHKERRQ(Ierr)

          else ! col is wellinj or wellprod, insert JacA%Val(1,1,j)

             idxm(1) = row
             idxn(1) = col

             call MatSetValues(A_mpi, 1, idxm, 1, idxn, JacA%Val(1,1,j), INSERT_VALUES, Ierr)
             CHKERRQ(Ierr)
          end if
       end do
    end do

    call MatAssemblyBegin(A_mpi,MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)
    call MatAssemblyEnd(A_mpi,MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)

    !write(*,*) '>>>>>>>>>> A is set <<<<<<<<<<'
    !call MatView(A_mpi,PETSC_VIEWER_STDOUT_WORLD,Ierr)

  end subroutine SolvePetsc_SetAmpi


  ! Set values Ap: part pression of A_mpi
  subroutine SolvePetsc_SetAp

    integer :: i, j
    integer :: row, col
    PetscErrorCode :: Ierr

    ! cols of JacA are (node, frac)
    ! insert value JacA%Val(1,1,*)
    do i=1, JacA%Nb

       row = RowLToRowGBlock(i) - 1             ! 0-based in petsc

       do j= JacA%Pt(i)+1, JacA%Pt(i+1)

          col = ColLToColGBlock( JacA%Num(j) ) - 1 ! 0-based in petsc

          call MatSetValue(Ap, row, col, JacA%Val(1,1,j), INSERT_VALUES, Ierr)
          CHKERRQ(Ierr)
       end do
    end do

    call MatAssemblyBegin(Ap, MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)
    call MatAssemblyEnd(Ap, MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)

    ! call MatView(Ap,PETSC_VIEWER_STDOUT_WORLD,Ierr)

  end subroutine SolvePetsc_SetAp


#ifdef _THERMIQUE_

  ! Set values At, temperature part of A_mpi
  subroutine SolvePetsc_SetAt

    ! TODO: add well

    integer :: i, j
    integer :: row, col
    PetscErrorCode :: Ierr

    logical :: IsTprimbyContext(NbContexte)

    do i=1, NbContexte
       do j=1, NbIncTotalPrim_ctx(i)
          if(psprim(j,i)==2) then
             IsTprimbyContext(i) = .true. ! For context i, T is prim
          else
             IsTprimbyContext(i) = .false. ! For context i, T is secd
          end if
       enddo
    end do

    ! T is prim or secd for unknowns (node and frac)
    ! IsTprimNodeFracOwn: unknowns (node own, frac own)
    ! IsTprimNodeFracLocal: unknowns (node local, frac local)
    do i=1, NbNodeOwn_Ncpus(commRank+1)
       IsTprimNodeFracOwn(i) = IsTprimbyContext(IncNode(i)%ic)
    end do
    do i=1, NbFracOwn_Ncpus(commRank+1)
       IsTprimNodeFracOwn(NbNodeOwn_Ncpus(commRank+1)+i) &
            = IsTprimbyContext(IncFrac(i)%ic)
    end do

    do i=1, NbNodeLocal_Ncpus(commRank+1)
       IsTprimNodeFracLocal(i) = IsTprimbyContext(IncNode(i)%ic)
    end do
    do i=1, NbFracLocal_Ncpus(commRank+1)
       IsTprimNodeFracLocal(NbNodeLocal_Ncpus(commRank+1)+i) &
            = IsTprimbyContext(IncFrac(i)%ic)
    end do

    ! cols of JacA are (node, frac)
    ! insert value JacA%Val(2,2,*)
    ! if T is send for i, At(i,:)=0, At(:,i)=0, At(i,i)=1
    do i=1, JacA%Nb

       row = RowLToRowGBlock(i) - 1 ! 0-based in petsc

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then ! for i, T is prim

          do j= JacA%Pt(i)+1, JacA%Pt(i+1)

             col = ColLToColGBlock( JacA%Num(j) ) - 1 ! 0-based in petsc

             if(IsTprimNodeFracLocal( JacA%Num(j) ) .eqv. .true.) then ! for j, T is prim
#ifdef _WIP_FREEFLOW_STRUCTURES_
                if( JacA%Num(j)<=NbNodeLocal .and. IdFFNodeLocal(JacA%Num(j)) ) then  ! FIXME: FreeFlow node, T is first inc (not always)
                   call MatSetValue(At, row, col, JacA%Val(1,1,j), INSERT_VALUES, Ierr)
                   CHKERRQ(Ierr)
                   call CommonMPI_abort('in solvePetsc entered in new loop (_WIP_FREEFLOW_STRUCTURES_)')
                else ! reservoir node, T is second inc
                   call MatSetValue(At, row, col, JacA%Val(2,2,j), INSERT_VALUES, Ierr)
                   CHKERRQ(Ierr)
                endif
#else
                ! reservoir node, T is second inc
                call MatSetValue(At, row, col, JacA%Val(2,2,j), INSERT_VALUES, Ierr)
                CHKERRQ(Ierr)
#endif
             else
                call MatSetValue(At, row, col, 0.d0, INSERT_VALUES, Ierr)
                CHKERRQ(Ierr)
             end if
          end do

       else ! for i, T is secd

          do j= JacA%Pt(i)+1, JacA%Pt(i+1)

             col = ColLToColGBlock( JacA%Num(j) ) - 1 ! 0-based in petsc

             if(row==col) then
                call MatSetValue(At, row, col, 1.d0, INSERT_VALUES, Ierr)
                CHKERRQ(Ierr)
             else
                call MatSetValue(At, row, col, 0.d0, INSERT_VALUES, Ierr)
                CHKERRQ(Ierr)
             end if
          end do

       end if
    end do

    call MatAssemblyBegin(At, MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)
    call MatAssemblyEnd(At, MAT_FINAL_ASSEMBLY,Ierr)
    CHKERRQ(Ierr)

    ! call MatView(At,PETSC_VIEWER_STDOUT_WORLD,Ierr)

  end subroutine SolvePetsc_SetAt

#endif

  subroutine SolvePetsc_Ksp_configuration(&
    tolerance, maximum_nb_iterations, restart_iteration) &
  bind(C, name="SolvePetsc_Ksp_configuration")

    real(c_double), intent(in), value :: tolerance
    integer(c_int), intent(in), value :: maximum_nb_iterations
    integer(c_int), intent(in), value :: restart_iteration
    PetscErrorCode :: Ierr

    ! CHECKME: are there border effects in assigning these global variables?
    PetscKspTol = tolerance
    kspitmax = maximum_nb_iterations
    ! CHECKME: Cf. PETSc doc
    ! We are supposed to use a left PC ?
    ! So the good choice would be KSP_NORM_PRECONDITIONED 
    ! call KSPSetNormType(ksp_mpi, KSP_NORM_PRECONDITIONED, Ierr)
    call KSPSetNormType(ksp_mpi, KSP_NORM_UNPRECONDITIONED, Ierr)
    CHKERRQ(Ierr)
    ! CHECKME: why no atol???
    call KSPSetTolerances(ksp_mpi, PetscKspTol, PETSC_DEFAULT_REAL, 1.d10, kspitmax, Ierr)
    CHKERRQ(Ierr)
    if(.not. allocated(kspHistory)) then
       allocate(kspHistory(kspitmax+1))
    else
        if(size(kspHistory)<kspitmax+1) then
           deallocate(kspHistory)
           allocate(kspHistory(kspitmax+1))
        end if
    end if
    call KSPSetResidualHistory(ksp_mpi, kspHistory, kspitmax, PETSC_TRUE, Ierr)
    CHKERRQ(Ierr)
    ! CHECKME: no restart!!!
    call KSPGMRESSetRestart(ksp_mpi, restart_iteration, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_Ksp_configuration

  subroutine SolvePetsc_CreateKsp_direct_solver
  
    PetscErrorCode :: Ierr

    call KSPCreate(ComPASS_COMM_WORLD, ksp_mpi, Ierr); CHKERRQ(Ierr)
    call KSPSetOperators(ksp_mpi, A_mpi, A_mpi, Ierr); CHKERRQ(Ierr)
    call KSPSetType(ksp_mpi, KSPPREONLY, Ierr); CHKERRQ(Ierr)
    call KSPGetPC(ksp_mpi, pc_mpi, Ierr); CHKERRQ(Ierr)
    call PCSetType(pc_mpi, PCLU, Ierr); CHKERRQ(Ierr)

    call KSPSetFromOptions(ksp_mpi, Ierr); CHKERRQ(Ierr)

  end subroutine SolvePetsc_CreateKsp_direct_solver

  subroutine SolvePetsc_CreateKsp
  
    PetscErrorCode :: Ierr

    call KSPCreate(ComPASS_COMM_WORLD, ksp_mpi, Ierr); CHKERRQ(Ierr)
    call KSPSetOperators(ksp_mpi, A_mpi, A_mpi, Ierr); CHKERRQ(Ierr)
    
    ! CHECKME: Cf. PETSc doc
    ! Normally, it is best to use the KSPSetFromOptions() command and
    ! then set the KSP type from the options database
    call KSPSetType(ksp_mpi, KSPGMRES, Ierr); CHKERRQ(Ierr)
    call SolvePetsc_Ksp_configuration(PetscKspTol, kspitmax, kspitmax)

    ! call KSPMonitorSet(ksp_mpi, KSPMonitorDefault, &
    !      PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, Ierr); CHKERRQ(Ierr)

    call KSPGetPC(ksp_mpi, pc_mpi, Ierr); CHKERRQ(Ierr)

    ! call PCSetType(pc_mpi,PCNONE,Ierr); CHKERRQ(Ierr)

    ! call PCSetType(pc_mpi, PCLU, Ierr); CHKERRQ(Ierr)

    call PCSetType(pc_mpi, PCHYPRE, Ierr); CHKERRQ(Ierr)
    ! ! Euclid in hypre: ILU(k), k is level
    ! call PCHYPRESetType(pc_mpi, "euclid", Ierr); CHKERRQ(Ierr)
    ! call PetscOptionsSetValue("-pc_hypre_euclid_levels", "0", Ierr); CHKERRQ(Ierr)
    ! call PetscOptionsSetValue("-pc_hypre_euclid_bj","1",Ierr); CHKERRQ(Ierr)
    call PCHYPRESetType(pc_mpi,'boomeramg',Ierr); CHKERRQ(Ierr)

    call KSPSetFromOptions(ksp_mpi, Ierr); CHKERRQ(Ierr)

  end subroutine SolvePetsc_CreateKsp


  ! Set up ksp for CPR-AMG preconditioner
  subroutine SolvePetsc_cpramgCreateKsp

    PetscErrorCode :: Ierr

    ! 1. KSP
    ! 2. PC Pression
    ! 3. PC ilu0
    ! 4. PC Temperature (if thermal)

    ! 1. ksp solver create

    ! FIXME: the following lines are duplicated in SolvePetsc_CreateKsp
 
    ! ksp solver create
    call KSPCreate(ComPASS_COMM_WORLD, ksp_mpi, Ierr); CHKERRQ(Ierr)
    call KSPSetOperators(ksp_mpi, A_mpi, A_mpi, Ierr); CHKERRQ(Ierr)
    
    ! CHECKME: Cf. PETSc doc
    ! Normally, it is best to use the KSPSetFromOptions() command and
    ! then set the KSP type from the options database
    call KSPSetType(ksp_mpi, KSPGMRES, Ierr); CHKERRQ(Ierr)
    call SolvePetsc_Ksp_configuration(PetscKspTol, kspitmax, kspitmax)
    
    ! ! monitor
    ! call KSPMonitorSet(ksp_mpi, KSPMonitorDefault, &
    !      PETSC_NULL_OBJECT, PETSC_NULL_FUNCTION, Ierr)
    ! CHKERRQ(Ierr)

    ! set preconditioner CPR-AMG
    call KSPGetPC(ksp_mpi, pc_mpi, Ierr)
    CHKERRQ(Ierr)
    call PCSetType(pc_mpi, PCSHELL, Ierr)
    CHKERRQ(Ierr)

    ! apply shell pc
#ifdef _THERMIQUE_


    ! call PCShellSetApply(pc_mIpi, SolvePetsc_cpramgPCApply_P_additive, Ierr)
    call PCShellSetApply(pc_mpi, SolvePetsc_cpramgPCApply_P_multiplicative, Ierr)

    ! call PCShellSetApply(pc_mpi, SolvePetsc_cpramgPCApply_PT_additive, Ierr)

    ! call PCShellSetApply(pc_mpi, SolvePetsc_cpramgPCApply_PT_multiplicative, Ierr)
    ! call PCShellSetApply(pc_mpi, SolvePetsc_cpramgPCApply_T_multiplicative, Ierr)
    CHKERRQ(Ierr)
#else
    call PCShellSetApply(pc_mpi, SolvePetsc_cpramgPCApply_P_multiplicative, Ierr)
    CHKERRQ(Ierr)
#endif

    ! setup shell pc and solver
    call PCShellSetSetUp(pc_mpi, SolvePetsc_cpramgPCSetUp,Ierr)
    CHKERRQ(Ierr)

    ! Set KSP options from the options database
    call KSPSetFromOptions(ksp_mpi, Ierr)
    CHKERRQ(Ierr)


    ! 2. PC Pression
    call PCCreate(ComPASS_COMM_WORLD, pcamg_p, Ierr)
    CHKERRQ(Ierr)
    call PCSetOperators(pcamg_p, Ap, Ap, Ierr)
    CHKERRQ(Ierr)
    call PCSetType(pcamg_p, PCHYPRE, Ierr)
    CHKERRQ(Ierr)
    call PCHYPRESetType(pcamg_p, "boomeramg", Ierr)
    CHKERRQ(Ierr)

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
    call PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.5",Ierr)
#else
    call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_hypre_boomeramg_strong_threshold","0.5",Ierr)
#endif
    CHKERRQ(Ierr)

    call PCSetFromOptions(pcamg_p, Ierr)
    CHKERRQ(Ierr)


    ! 3. PC ilu0
    call PCCreate(ComPASS_COMM_WORLD, pcilu0, Ierr)
    CHKERRQ(Ierr)
    call PCSetOperators(pcilu0, A_mpi, A_mpi, Ierr)
    CHKERRQ(Ierr)

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
    call PCSetType(pcilu0, PCHYPRE, Ierr)
    CHKERRQ(Ierr)
    call PCHYPRESetType(pcilu0, "euclid", Ierr)
    CHKERRQ(Ierr)
    call PetscOptionsSetValue("-pc_hypre_euclid_levels", "0", Ierr)
    CHKERRQ(Ierr)
    call PetscOptionsSetValue("-pc_hypre_euclid_bj","1",Ierr)
    CHKERRQ(Ierr)    
#else
    call PCSetType(pcilu0, PCBJACOBI, Ierr) 
    CHKERRQ(Ierr)
#endif

    call PCSetFromOptions(pcilu0, Ierr)
    CHKERRQ(Ierr)

! #ifdef _THERMIQUE_

    ! ! 4. PC Temperature
    ! call PCCreate(ComPASS_COMM_WORLD, pcamg_t, Ierr)
    ! CHKERRQ(Ierr)
    ! call PCSetOperators(pcamg_t, At, At, Ierr)
    ! CHKERRQ(Ierr)

    ! call PCSetType(pcamg_t, PCHYPRE, Ierr)
    ! CHKERRQ(Ierr)
    ! call PCHYPRESetType(pcamg_t, "boomeramg", Ierr)
    ! CHKERRQ(Ierr)
    ! call PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.5",Ierr)
    ! CHKERRQ(Ierr)
    ! call PCSetFromOptions(pcamg_t, Ierr)
    ! CHKERRQ(Ierr)

    ! call PCSetType(pcamg_t, PCLU, Ierr)
    ! CHKERRQ(Ierr)
    ! call PCFactorSetMatSolverPackage(pcamg_t, MATSOLVERSUPERLU_DIST, Ierr)
    ! CHKERRQ(Ierr)
    ! call PCFactorSetMatSolverPackage(pcamg_t, MATSOLVERMUMPS, Ierr)
    ! CHKERRQ(Ierr)

! #endif

  end subroutine SolvePetsc_cpramgCreateKsp


  ! *** AMG for pression first, ILU0 second *** !
  !   P1^{-1}: amg for pression
  !   P2^{-1}: ilu0 for all

  ! Preconditioner: v -> P^{-1}v, no thermal
  ! P^{-1}v = P1^{-1}v + P2^{-1}(v-A_mpi*P_1^{-1}*v)
  ! Ps: v must not be modified
  subroutine SolvePetsc_cpramgPCApply_P_multiplicative(pcin, v, Pv, Ierr)

    PC, intent(inout) :: pcin
    Vec, intent(inout) :: v, Pv
    PetscErrorCode, intent(inout) :: Ierr

    double precision, pointer :: ptr1(:)
    double precision, pointer :: ptr2(:)
    integer :: NbNodeFrac, NbWell, i, iv

    NbNodeFrac = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)
    NbWell = NbWellInjOwn_Ncpus(commRank+1) + NbWellProdOwn_Ncpus(commRank+1)

    ! Pression part of v:
    ! v_pt = R1*v, R1: restriction matrix to pression
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, NbNodeFrac
       iv = (i-1)*NbCompThermique + 1
       ptr2(i) = ptr1(iv)
    end do

    do i=1, NbWell
       iv = NbNodeFrac*NbCompThermique + i
       ptr2(i+NbNodeFrac) = ptr1(iv)
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_1^{-1} v_pt, AMG
    call PCApply(pcamg_p, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    ! P1v has been initialized as 0
    ! The part hors pression is always 0
    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, NbNodeFrac
       iv = (i-1)*NbCompThermique + 1
       ptr2(iv) = ptr1(i)
    end do

    do i=1, NbWell
       iv = NbNodeFrac*NbCompThermique + i
       ptr2(iv) = ptr1(i+NbNodeFrac)
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = A*P1v
    call MatMult(A_mpi, P1v, AP1v, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = v - AP1v
    call VecAYPX(AP1v, -1.d0, v, Ierr)
    CHKERRQ(Ierr)

    ! Pv = P_2^{-1} AP1v, ILU(0)
    call PCApply(pcilu0, AP1v, Pv, Ierr)
    CHKERRQ(Ierr)

    ! Pv = Pv + P1v
    call VecAXPY(Pv, 1.d0, P1v, Ierr)
    CHKERRQ(Ierr)


    ! ! *** ILU0 first, AMG second *** !
    ! !   P1^{-1}: ilu0 for all
    ! !   P2^{-1}: amg for pression

    ! ! Pv = P_2^{-1} v, ILU(0)
    ! call PCApply(pcilu0, v, Pv, Ierr)
    ! CHKERRQ(Ierr)

    ! call MatMult(A_mpi, Pv, AP1v, Ierr)
    ! CHKERRQ(Ierr)

    ! call VecAYPX(AP1v, -1.d0, v, Ierr)
    ! CHKERRQ(Ierr)

    ! call VecGetArrayReadF90(AP1v, ptr1, Ierr)
    ! CHKERRQ(Ierr)
    ! call VecGetArrayF90(v_pt, ptr2, Ierr)
    ! CHKERRQ(Ierr)

    ! do i=1, Nb
    !    iv = (i-1)*NbCompThermique + 1
    !    ptr2(i) = ptr1(iv)
    ! end do

    ! call VecRestoreArrayReadF90(AP1v, ptr1, Ierr)
    ! CHKERRQ(Ierr)
    ! call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    ! CHKERRQ(Ierr)

    ! ! P1v_pt = P_1^{-1} v_pt, AMG
    ! call PCApply(pcamg_p, v_pt, P1v_pt, Ierr)
    ! CHKERRQ(Ierr)


    ! call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    ! CHKERRQ(Ierr)
    ! call VecGetArrayF90(Pv, ptr2, Ierr)
    ! CHKERRQ(Ierr)

    ! do i=1, Nb
    !       iv = (i-1)*NbCompThermique + 1
    !       ptr2(iv) = ptr2(iv) + ptr1(i)
    ! end do

    ! call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    ! CHKERRQ(Ierr)
    ! call VecRestoreArrayF90(Pv, ptr2, Ierr)
    ! CHKERRQ(Ierr)

  end subroutine SolvePetsc_cpramgPCApply_P_multiplicative


  ! *** AMG for pression and ILU0 second (additive) *** !
  !   P1^{-1}: amg for pression
  !   P2^{-1}: ilu0 for all
  ! Preconditioner: v -> P^{-1}v, thermal
  ! Ps: v must not be modified
  subroutine SolvePetsc_cpramgPCApply_P_additive(pcin, v, pv, Ierr)

    PC, intent(inout) :: pcin
    Vec, intent(inout) :: v, pv
    PetscErrorCode, intent(inout) :: Ierr

    double precision, pointer :: ptr1(:)
    double precision, pointer :: ptr2(:)
    integer :: Nb, i, iv

    Nb = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)

    ! Pv = P_2^{-1} v, ILU(0)
    call PCApply(pcilu0, v, Pv, Ierr)
    CHKERRQ(Ierr)

    ! Pression part of v:
    ! v_pt = R1*v, R1: restriction matrix to pression
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(i) = ptr1(iv)
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_1^{-1} v_pt, AMG
    call PCApply(pcamg_p, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    ! P1v has been initialized as 0
    ! The part hors pression and temperature is always 0
    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(Pv, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(iv) = ptr2(iv) + ptr1(i)
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_cpramgPCApply_P_additive


#ifdef _THERMIQUE_


  ! *** AMG for pression first, AMG for temperature second, ILU0 third *** !
  !   P1^{-1}: amg for pression
  !   P2^{-1}: amg for temperature
  !   P3^{-1}: ilu0 for all
  ! Preconditioner: v -> P^{-1}v, thermal
  ! Ps: v must not be modified
  subroutine SolvePetsc_cpramgPCApply_PT_multiplicative(pcin, v, pv, Ierr)

    PC, intent(inout) :: pcin
    Vec, intent(inout) :: v, pv
    PetscErrorCode, intent(inout) :: Ierr

    double precision, pointer :: ptr1(:)
    double precision, pointer :: ptr2(:)
    integer :: Nb, i, iv

    Nb = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)

    ! Pression part of v
    ! v_pt = R1*v, R1: restriction matrix to pression
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(i) = ptr1(iv)
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_1^{-1} v_pt, AMG
    call PCApply(pcamg_p, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    ! P1v has been initialized as 0
    ! The part hors pression is always 0
    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(iv) = ptr1(i)
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = A*P1v
    call MatMult(A_mpi, P1v, AP1v, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = v - AP1v
    call VecAYPX(AP1v, -1.d0, v, Ierr)
    CHKERRQ(Ierr)

    ! Temperature part of AP1v
    ! AP1v_t = R2*AP1v, R2: restriction matrix to temperature
    call VecGetArrayReadF90(AP1v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(AP1v_t, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(i) = ptr1(iv)
       else
          ptr2(i) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(AP1v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(AP1v_t, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P2AP1v_t = P_2^{-1} AP1v_t, AMG
    call PCApply(pcamg_t, AP1v_t, P2AP1v_t, Ierr)
    CHKERRQ(Ierr)

    ! P2AP1v has been initialized as 0
    ! The part hors temperature is always 0
    call VecGetArrayReadF90(P2AP1v_t, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P2AP1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(iv) = ptr1(i)
       else
          ptr2(iv) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(P2AP1v_t, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P2AP1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! AP2AP1v = A * P2AP1v
    call MatMult(A_mpi, P2AP1v, AP2AP1v, Ierr)
    CHKERRQ(Ierr)

    ! AP2AP1v = AP1v - AP2AP1v
    call VecAYPX(AP2AP1v, -1.d0, AP1v, Ierr)
    CHKERRQ(Ierr)

    ! Pv = P_3^{-1} AP2AP1v, ILU(0)
    call PCApply(pcilu0, AP2AP1v, Pv, Ierr)
    CHKERRQ(Ierr)

    ! Pv = Pv + P2AP1v
    call VecAXPY(Pv, 1.d0, P2AP1v, Ierr)
    CHKERRQ(Ierr)

    ! Pv = Pv + P1v
    call VecAXPY(Pv, 1.d0, P1v, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_cpramgPCApply_PT_multiplicative


  ! *** AMG for pression and AMG for temperature first (additive), ILU0 second (multiplicative) *** !
  !   P1^{-1}: amg for pression
  !   P2^{-1}: amg for temperature
  !   P3^{-1}: ilu0 for all
  ! Preconditioner: v -> P^{-1}v, thermal
  ! Ps: v must not be modified
  subroutine SolvePetsc_cpramgPCApply_PT_additive(pcin, v, pv, Ierr)

    PC, intent(inout) :: pcin
    Vec, intent(inout) :: v, pv
    PetscErrorCode, intent(inout) :: Ierr

    double precision, pointer :: ptr1(:)
    double precision, pointer :: ptr2(:)
    integer :: Nb, i, iv

    Nb = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)

    ! Pression part of v:
    ! v_pt = R1*v, R1: restriction matrix to pression
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(i) = ptr1(iv)
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_1^{-1} v_pt, AMG
    call PCApply(pcamg_p, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    ! P1v has been initialized as 0
    ! The part hors pression and temperature is always 0
    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 1
       ptr2(iv) = ptr1(i)
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)


    ! Temperature part of v:
    ! v_pt = R2*v, R2: restriction matrix to temperature
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(i) = ptr1(iv)
       else
          ptr2(i) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_2^{-1} v_pt, AMG
    call PCApply(pcamg_t, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(iv) = ptr1(i)
       else
          ptr2(iv) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = A*P1v
    call MatMult(A_mpi, P1v, AP1v, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = v - AP1v
    call VecAYPX(AP1v, -1.d0, v, Ierr)
    CHKERRQ(Ierr)

    ! Pv = P_3^{-1} AP1v, ILU(0)
    call PCApply(pcilu0, AP1v, Pv, Ierr)
    CHKERRQ(Ierr)

    ! Pv = Pv + P1v
    call VecAXPY(Pv, 1.d0, P1v, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_cpramgPCApply_PT_additive


  ! *** AMG for temperature, ILU0 second (multiplicative) *** !
  !   P2^{-1}: amg for temperature
  !   P3^{-1}: ilu0 for all
  ! Preconditioner: v -> P^{-1}v, thermal
  ! Ps: v must not be modified
  subroutine SolvePetsc_cpramgPCApply_T_multiplicative(pcin, v, pv, Ierr)

    PC, intent(inout) :: pcin
    Vec, intent(inout) :: v, pv
    PetscErrorCode, intent(inout) :: Ierr

    double precision, pointer :: ptr1(:)
    double precision, pointer :: ptr2(:)
    integer :: Nb, i, iv

    Nb = NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)

    ! Temperature part of v:
    ! v_pt = R2*v, R2: restriction matrix to temperature
    call VecGetArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(i) = ptr1(iv)
       else
          ptr2(i) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(v, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(v_pt, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! P1v_pt = P_2^{-1} v_pt, AMG
    call PCApply(pcamg_t, v_pt, P1v_pt, Ierr)
    CHKERRQ(Ierr)

    call VecGetArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecGetArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    do i=1, Nb
       iv = (i-1)*NbCompThermique + 2

       if(IsTprimNodeFracOwn(i) .eqv. .true.) then
          ptr2(iv) = ptr1(i)
       else
          ptr2(iv) = 0.d0
       end if
    end do

    call VecRestoreArrayReadF90(P1v_pt, ptr1, Ierr)
    CHKERRQ(Ierr)
    call VecRestoreArrayF90(P1v, ptr2, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = A*P1v
    call MatMult(A_mpi, P1v, AP1v, Ierr)
    CHKERRQ(Ierr)

    ! AP1v = v - AP1v
    call VecAYPX(AP1v, -1.d0, v, Ierr)
    CHKERRQ(Ierr)

    ! Pv = P_3^{-1} AP1v, ILU(0)
    call PCApply(pcilu0, AP1v, Pv, Ierr)
    CHKERRQ(Ierr)

    ! Pv = Pv + P1v
    call VecAXPY(Pv, 1.d0, P1v, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_cpramgPCApply_T_multiplicative

#endif


  ! Create Sm_mpi and x_mpi
  subroutine SolvePetsc_CreateSm

    PetscErrorCode :: Ierr

    ! create Sm_mpi
    call VecCreateMPI(ComPASS_COMM_WORLD, &
         NrowL, NrowG, &
         Sm_mpi, Ierr)
    CHKERRQ(Ierr)

    call VecSet(Sm_mpi, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(Sm_mpi, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(Sm_mpi, Ierr)
    CHKERRQ(Ierr)

   !  ! create x_mpi
   !  call VecDuplicate(Sm_mpi, x_mpi, Ierr)
   !  CHKERRQ(Ierr)

   !  call VecSet(x_mpi, 0.d0, Ierr)
   !  CHKERRQ(Ierr)
   !  call VecAssemblyBegin(x_mpi, Ierr)
   !  CHKERRQ(Ierr)
   !  call VecAssemblyEnd(x_mpi, Ierr)
   !  CHKERRQ(Ierr)

    ! create y_mpi
    call VecDuplicate(Sm_mpi, y_mpi, Ierr)
    CHKERRQ(Ierr)

    call VecSet(y_mpi, 0.d0, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyBegin(y_mpi, Ierr)
    CHKERRQ(Ierr)
    call VecAssemblyEnd(y_mpi, Ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_CreateSm


  ! Set Sm_mpi
  subroutine SolvePetsc_SetSm

    integer :: i, j
    integer :: start, blockstart
    double precision, dimension(:), pointer :: ptr
    PetscErrorCode :: Ierr

    ! set Sm to Sm_mpi
    call VecGetArrayF90(Sm_mpi, ptr, Ierr)
    CHKERRQ(Ierr)

    do i=1, NbNodeOwn+NbFracOwn

       do j=1, NbCompThermique
          ptr((i-1)*NbCompThermique+j) = Sm(j,i)
       end do
    end do

    start = (NbNodeOwn + NbFracOwn) * NbCompThermique
    blockstart = NbNodeOwn + NbFracOwn

    do i=1, NbWellInjOwn+NbWellProdOwn
       ptr(i+start) = Sm(1,i+blockstart)
    end do

    call VecRestoreArrayF90(Sm_mpi, ptr, Ierr)
    CHKERRQ(Ierr)


    ! call VecSetValues(Sm_mpi, JacA%Nb, Sm_ix, Sm, INSERT_VALUES, Ierr)

    ! call VecAssemblyBegin(Sm_mpi, Ierr)
    ! CHKERRQ(Ierr)
    ! call VecAssemblyEnd(Sm_mpi, Ierr)
    ! CHKERRQ(Ierr)

    ! call VecView(Sm_mpi, PETSC_VIEWER_STDOUT_WORLD, Ierr)

  end subroutine SolvePetsc_SetSm


  ! Create Ap and At
  subroutine SolvePetsc_cpramgCreateApAt

    integer, allocatable, dimension(:) :: &
         d_nnz, & ! number of nonzeros per row in diagonal portion
         o_nnz    ! number of nonzeros per row in the off-diagonal portion

    ! some tmp values, explained when setting values
    integer :: Nl_Fo, Nl_Fl, Nl_Fl_WIo, Nl_Fl_WIl, Nl_Fl_WIl_WPo, Nl_Fl_WIl_WPl

    ! tmp values
    integer :: i, j, lj
    PetscErrorCode :: Ierr

    ! number of nonzeros per row in diag or off-diag portion
    allocate(d_nnz(NBlockrowL))
    allocate(o_nnz(NBlockrowL))
    d_nnz(:) = 0
    o_nnz(:) = 0

    ! N: Node, F: Frac, WI: well inj, WP: well prod
    ! l: local, o: own
    Nl_Fo = NbNodeLocal + NbFracOwn
    Nl_Fl = NbNodeLocal + NbFracLocal
    Nl_Fl_WIo = NbNodeLocal + NbFracLocal + NbWellInjOwn
    Nl_Fl_WIl = NbNodeLocal + NbFracLocal + NbWellInjLocal
    Nl_Fl_WIl_WPo = NbNodeLocal + NbFracLocal + NbWellInjLocal + NbWellProdOwn
    Nl_Fl_WIl_WPl = NbNodeLocal + NbFracLocal + NbWellInjLocal + NbWellProdLocal

    do i=1, NBlockrowL

       do j=JacA%Pt(i)+1,JacA%Pt(i+1)

          lj = JacA%Num(j)

          if (lj<=NbNodeOwn) then ! node own
             d_nnz(i) = d_nnz(i) + 1
          else if (NbNodeOwn<lj .and. lj<=NbNodeLocal) then ! node ghost
             o_nnz(i) = o_nnz(i) + 1
          else if (NbNodeLocal<lj .and. lj<=Nl_Fo) then
             d_nnz(i) = d_nnz(i) + 1
          else if(Nl_Fo<lj .and. lj<=Nl_Fl) then
             o_nnz(i) = o_nnz(i) + 1
          else if(Nl_Fl<lj .and. lj<=Nl_Fl_WIo) then
             d_nnz(i) = d_nnz(i) + 1
          else if(Nl_Fl_WIo<lj .and. lj<=Nl_Fl_WIl) then
             o_nnz(i) = o_nnz(i) + 1
          else if(Nl_Fl_WIl<lj .and. lj<=Nl_Fl_WIl_WPo) then
             d_nnz(i) = d_nnz(i) + 1
          else if(Nl_Fl_WIl_WPo<lj .and. lj<=Nl_Fl_WIl_WPl) then
             o_nnz(i) = o_nnz(i) + 1
          end if
       end do
    end do

    ! create matrix Ap and At
    call MatCreateAIJ(ComPASS_COMM_WORLD, &
         NBlockrowL, NBlockcolL, &
         NBlockrowG, NBlockcolG, &
         0, d_nnz, 0, o_nnz, Ap, Ierr)
    CHKERRQ(Ierr)

! #ifdef _THERMIQUE_

!     call MatCreateAIJ(ComPASS_COMM_WORLD, &
!          NBlockrowL, NBlockcolL, &
!          NBlockrowG, NBlockcolG, &
!          0, d_nnz, 0, o_nnz, At, Ierr)
!     CHKERRQ(Ierr)
! #endif

    deallocate(d_nnz)
    deallocate(o_nnz)

  end subroutine SolvePetsc_CpramgCreateApAt

  subroutine SolvePetsc_dump_system(basename)
  
    character(len=*), intent(in) :: basename

    PetscViewer :: viewer
    PetscErrorCode :: Ierr

    call PetscViewerASCIIOpen(ComPASS_COMM_WORLD, trim(basename) // '_structure.dat', viewer, ierr)
    CHKERRQ(Ierr)
    call KSPView(ksp_mpi, viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerFlush(viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerDestroy(viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerASCIIOpen(ComPASS_COMM_WORLD, trim(basename) // '_A.dat', viewer, ierr)
    CHKERRQ(Ierr)
    call MatView(A_mpi, viewer, Ierr)
    CHKERRQ(Ierr)
    call PetscViewerFlush(viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerDestroy(viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerASCIIOpen(ComPASS_COMM_WORLD, trim(basename) // '_b.dat', viewer, ierr)
    CHKERRQ(Ierr)
    call VecView(Sm_mpi, viewer, Ierr)
    CHKERRQ(Ierr)
    call PetscViewerFlush(viewer, ierr)
    CHKERRQ(Ierr)
    call PetscViewerDestroy(viewer, ierr)
    CHKERRQ(Ierr)

  end subroutine SolvePetsc_dump_system

  function SolvePetsc_KspSolveIterationNumber() &
    result(NkspIter) &
    bind(C, name="SolvePetsc_KspSolveIterationNumber")
    
    integer(c_int) :: NkspIter
    PetscErrorCode :: Ierr

    call KSPGetIterationNumber(ksp_mpi, NkspIter, Ierr); CHKERRQ(Ierr)
    NkspIter = NkspIter + 1 ! CHECKME: why + 1?
       
  end function SolvePetsc_KspSolveIterationNumber

  subroutine SolvePetsc_Ksp_history(output)

    real(c_double), dimension(:), intent(inout) :: output
    integer :: i, n
    
    n = SolvePetsc_KspSolveIterationNumber()
    if(.not. allocated(kspHistory)) then
        call CommonMPI_abort("kspHistory should be allocated")
    end if
    if(n>size(output)) then
        call CommonMPI_abort("inconsistent sizes")
    end if
    do i=1, n
       output(i) = kspHistory(i)
    end do
    do i=n+1, size(output)
       output(i) = 0.d0
    end do

  end subroutine SolvePetsc_Ksp_history

  subroutine SolvePetsc_KspSolveIterations(cp, n) &
    bind(C, name="SolvePetsc_KspSolveIterations")
  
    type(c_ptr), intent(in), value :: cp
    integer(c_int), intent(in), value :: n
    real(c_double), pointer :: residuals(:)
    
    call c_f_pointer(cp, residuals, [n])
    call SolvePetsc_Ksp_history(residuals)
    
  end subroutine SolvePetsc_KspSolveIterations

  function SolvePetsc_KspSolve(x) result(reason)
    Vec, intent(inout) :: x
    integer(c_int) :: reason
    KSPConvergedReason :: native_reason ! this wraps a C enum
    PetscErrorCode :: Ierr

    call KSPSolve(ksp_mpi, Sm_mpi, x, Ierr)
    CHKERRQ(Ierr)

    call KSPGetConvergedReason(ksp_mpi, native_reason, Ierr)
    reason = native_reason ! FIXME: avoid conversion... use petsc4py / enums ?
    CHKERRQ(Ierr)

  end function SolvePetsc_KspSolve

  subroutine SolvePetsc_check_solution(x)
    Vec, intent(in) :: x
    PetscReal :: a
    PetscErrorCode :: Ierr

    call VecCopy(Sm_mpi, y_mpi, Ierr); CHKERRQ(Ierr)
    call VecScale(y_mpi, -1.d0, Ierr); CHKERRQ(Ierr)
    call MatMultAdd(A_mpi, x, y_mpi, y_mpi, Ierr); CHKERRQ(Ierr)

    write(*, *) 'linear solution check ||AX-b||' 
    call VecNorm(y_mpi, NORM_1, a, Ierr); CHKERRQ(Ierr)
    write(*, *) 'N1', a
    call VecNorm(y_mpi, NORM_2, a, Ierr); CHKERRQ(Ierr)
    write(*, *) 'N2', a
    call VecNorm(y_mpi, NORM_INFINITY, a, Ierr); CHKERRQ(Ierr)
    write(*, *) 'NI', a

  end subroutine SolvePetsc_check_solution 

  subroutine SolvePetsc_free

    PetscErrorCode :: Ierr

    ! Destroy
    call KSPDestroy(ksp_mpi, Ierr)
    CHKERRQ(Ierr)
    call MatDestroy(A_mpi, Ierr)
    CHKERRQ(Ierr)

    ! free ksp convergence history vector
    if(allocated(kspHistory)) then
       deallocate(kspHistory)
    end if

    ! free RowLToRowG ColLToColG
    deallocate(RowLToRowGBlock)
    deallocate(ColLToColGBlock)

    deallocate(RowLToRowG)
    deallocate(ColLToColG)

    call VecDestroy(Sm_mpi, Ierr)
    CHKERRQ(Ierr)
   !  call VecDestroy(x_mpi, Ierr)
   !  CHKERRQ(Ierr)

  end subroutine SolvePetsc_free


  ! Free
  subroutine SolvePetsc_cpramgFree

    PetscErrorCode :: Ierr

    call SolvePetsc_free

    ! Destroy
    call PCDestroy(pc_mpi, Ierr)
    CHKERRQ(Ierr)
    call KSPDestroy(ksp_mpi, Ierr)
    CHKERRQ(Ierr)
    call MatDestroy(A_mpi, Ierr)
    CHKERRQ(Ierr)

    call PCDestroy(pcamg_p, Ierr)
    CHKERRQ(Ierr)
    call PCDestroy(pcilu0, Ierr)
    CHKERRQ(Ierr)
    call MatDestroy(Ap, Ierr)
    CHKERRQ(Ierr)

    call VecDestroy(v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(P1v_pt, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(P1v, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(AP1v, Ierr)
    CHKERRQ(Ierr)

#ifdef _THERMIQUE_
    ! call PCDestroy(pcamg_t, Ierr)
    ! CHKERRQ(Ierr)
    ! call MatDestroy(At, Ierr)
    ! CHKERRQ(Ierr)

    call VecDestroy(AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(P2AP1v, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(P2AP1v_t, Ierr)
    CHKERRQ(Ierr)
    call VecDestroy(AP2AP1v, Ierr)
    CHKERRQ(Ierr)

    deallocate(IsTprimNodeFracOwn)
    deallocate(IsTprimNodeFracLocal)
#endif

  end subroutine SolvePetsc_cpramgFree


  ! row/col local to row/col global
  ! RowLtoRowG(i): gloabl row i in A_mpi where i is row in JacA
  ! ColLToColG(i): gloabl col i in A_mpi where i is col in JacA
  subroutine SolvePetsc_LtoG

    integer :: i, start, ndisp
    integer :: local_block_col_offset, global_col_offset(Ncpus)

    ! RowLtoRowG
    ! idea: add sum of numbers of own nodes/frac/well in the procs before commRank (rowstart)
    allocate( RowLToRowG(NBlockrowL) )

    ! node/frac
    do i=1, NbNodeOwn + NbFracOwn
       RowLToRowG(i) = rowstart(commRank+1) + (i-1)*NbCompThermique + 1
    end do

    ! well
    start = NbNodeOwn + NbFracOwn
    ndisp = rowstart(commRank+1) + (NbNodeOwn+NbFracOwn)*NbCompThermique

    do i=1, NbWellInjOwn + NbWellProdOwn
       RowLToRowG(i+start) = ndisp + i
    end do

    ! Colltocolg(i), where i is in (node local, frac local, well local)
    allocate( ColLToColG(NbNodeLocal+NbFracLocal + NbWellInjLocal +NbWellProdLocal) )
    ColLToColG(:) = 0

    local_block_col_offset = 0
    global_col_offset = colstart
    call SolvePetsc_LtoG_fill_ColLToColG(&
      local_block_col_offset, global_col_offset, NumNodebyProc, NbCompThermique, ColLToColG)
    local_block_col_offset = local_block_col_offset + NbNodeLocal
    global_col_offset = global_col_offset + NbCompThermique * NbNodeOwn_Ncpus
    if(NbFracLocal>0) then
      call SolvePetsc_LtoG_fill_ColLToColG(&
         local_block_col_offset, global_col_offset, NumFracbyProc, NbCompThermique, ColLToColG)
      local_block_col_offset = local_block_col_offset + NbFracLocal
    end if
    global_col_offset = global_col_offset + NbCompThermique * NbFracOwn_Ncpus
    if(NbWellInjLocal>0) then
      call SolvePetsc_LtoG_fill_ColLToColG(&
         local_block_col_offset, global_col_offset, NumWellInjbyProc, 1, ColLToColG)
      local_block_col_offset = local_block_col_offset + NbWellInjLocal
    endif
    global_col_offset = global_col_offset + NbWellInjOwn_Ncpus
    if(NbWellProdLocal>0) then
      call SolvePetsc_LtoG_fill_ColLToColG(&
         local_block_col_offset, global_col_offset, NumWellProdbyProc, 1, ColLToColG)
      ! local_block_col_offset = local_block_col_offset + NbWellProdLocal
    endif
    ! global_col_offset = global_col_offset + NbWellProdOwn_Ncpus

  end subroutine SolvePetsc_LtoG

  subroutine SolvePetsc_LtoG_fill_ColLToColG(local_block_offset, global_col_offset, dof_family, dof_size, LtoG_col_map)
    integer, intent(in) :: local_block_offset ! offset in terms of local blocks
    integer, dimension(Ncpus), intent(in) :: global_col_offset
    type(FamilyDOFIdCOC), intent(in) :: dof_family
    integer, intent(in) :: dof_size
    integer, dimension(:), intent(inout) :: LtoG_col_map

    integer :: i, proc, local_id

    do i=1, size(dof_family%ids)
      proc = dof_family%ids(i)%proc
      local_id = dof_family%ids(i)%local_id
      LtoG_col_map(local_block_offset + i) =  global_col_offset(proc+1) + (local_id - 1) * dof_size + 1
    end do

  end subroutine SolvePetsc_LtoG_fill_ColLToColG

  subroutine SolvePetsc_LtoGBlock

    integer :: i, local_block_col_offset, global_block_col_offset(Ncpus)

    ! RowLtoRowGBlock
    ! idea: add sum of numbers of own nodes/frac/well in the procs before commRank (rowstart)
    allocate( RowLToRowGBlock(NBlockrowL) )
    do i=1, NBlockrowL
       RowLToRowGBlock(i) = Blockrowstart(commRank+1) + i
    end do
    ! ColLtoColGBlock(i), where i is in (node local, frac local, well local)
    allocate( ColLToColGBlock(NbNodeLocal+NbFracLocal + NbWellInjLocal+NbWellProdLocal) )
    ColLToColGBlock(:) = 0

    local_block_col_offset = 0
    global_block_col_offset = Blockcolstart
    call SolvePetsc_LtoGBlock_fill_ColLToColGBlock( &
      local_block_col_offset, global_block_col_offset, NumNodebyProc, ColLToColGBlock)
    local_block_col_offset = local_block_col_offset + NbNodeLocal
    global_block_col_offset = global_block_col_offset + NbNodeOwn_Ncpus
    if(NbFracLocal>0) then
      call SolvePetsc_LtoGBlock_fill_ColLToColGBlock( &
         local_block_col_offset, global_block_col_offset, NumFracbyProc, ColLToColGBlock)
      local_block_col_offset = local_block_col_offset + NbFracLocal
      global_block_col_offset = global_block_col_offset + NbFracOwn_Ncpus
    end if
    if(NbWellInjLocal>0) then
      call SolvePetsc_LtoGBlock_fill_ColLToColGBlock( &
         local_block_col_offset, global_block_col_offset, NumWellInjbyProc, ColLToColGBlock)
      local_block_col_offset = local_block_col_offset + NbWellInjLocal
    end if
    global_block_col_offset = global_block_col_offset + NbWellInjOwn_Ncpus
    if(NbWellProdLocal>0) then
      call SolvePetsc_LtoGBlock_fill_ColLToColGBlock( &
         local_block_col_offset, global_block_col_offset, NumWellProdbyProc, ColLToColGBlock)
      ! local_block_col_offset = local_block_col_offset + NbWellProdLocal
    end if
    ! global_block_col_offset = global_block_col_offset + NbWellProdOwn_Ncpus

end subroutine SolvePetsc_LtoGBlock

!> Computes the offset in terms of dof_sizexdof_size small (dense) blocks
subroutine SolvePetsc_LtoGBlock_fill_ColLToColGBlock(local_col_offset, global_block_col_offset, dof_family, LtoG_block_col_map)
   integer, intent(in) :: local_col_offset
   integer, dimension(Ncpus), intent(in) :: global_block_col_offset
   type(FamilyDOFIdCOC), intent(in) :: dof_family
   integer, dimension(:), intent(inout) :: LtoG_block_col_map

   integer :: nb_proc, i, p, proc, local_id

   nb_proc = size(dof_family%offsets) - 1
   if(nb_proc<1) call CommonMPI_abort("Inconsistent number of procs")
   do p=1, nb_proc
      do i=dof_family%offsets(p)+1, dof_family%offsets(p+1)
         proc = dof_family%ids(i)%proc
         local_id = dof_family%ids(i)%local_id
         LtoG_block_col_map(local_col_offset + i) = global_block_col_offset(proc+1) + local_id
      end do
   end do

end subroutine SolvePetsc_LtoGBlock_fill_ColLToColGBlock

end module SolvePetsc

! FIXME: this is transitory
! This is out the module scope because of function names mangling
function compass_petsc_kspsolve(x) result(reason)

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petsc.h>
#endif

   use iso_c_binding, only: c_int
   use petsc
   use SolvePetsc, only: SolvePetsc_KspSolve

   implicit none

   Vec, intent(inout) :: x
   integer(c_int) :: reason

   reason = SolvePetsc_KspSolve(x)

end function compass_petsc_kspsolve

subroutine compass_check_solution(x)
#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petsc.h>
#endif

   use petsc
   use SolvePetsc, only: SolvePetsc_check_solution

   implicit none

   Vec, intent(in) :: x
   
   call SolvePetsc_check_solution(x)

end subroutine compass_check_solution
