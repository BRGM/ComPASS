!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecd

   use iso_c_binding, only: c_int, c_double, c_bool
   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use SchemeParameters, only: eps

   use DefModel, only: &
      NbComp, NbPhase, MCP, &
      NbIncTotalPrimMax, NbEqFermetureMax, NbIncTotalMax, &
      pschoice, psprim, pssecd, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      IndThermique, NbCompThermique

   use NumbyContext, only: &
      NbEqFermeture_ctx, NbIncTotalPrim_ctx, &
      NbIncPTC_ctx, NbIncTotal_ctx, NbCompCtilde_ctx

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, &
      IncCell, IncFrac, IncNode, &
      dPhasePressuredSNode, dPhasePressuredSFrac, dPhasePressuredSCell
   use MeshSchema, only: &
#ifdef _WITH_FREEFLOW_STRUCTURES_
      IsFreeflowNode, &
#endif
      NbCellLocal_Ncpus, NbFracLocal_Ncpus, NbNodeLocal_Ncpus

   use Newton, only: Newton_increments_pointers, Newton_increments, Newton_pointers_to_values
   use Thermodynamics, only: f_Fugacity
   use IncPrimSecdTypes, only: ControlVolumeInfo, IncPrimSecdTypes_collect_cv_info

   implicit none

   real(c_double), allocatable, dimension(:, :, :), target, public :: &
      dXssurdXpAll
   real(c_double), dimension(:, :, :), pointer, public :: &
      dXssurdXpCell, dXssurdXpFrac, dXssurdXpNode

   real(c_double), allocatable, dimension(:, :), target, public :: &
      SmdXsAll
   real(c_double), dimension(:, :), pointer, public :: &
      SmdXsCell, SmdXsFrac, SmdXsNode

   real(c_double), allocatable, dimension(:, :), target, public :: &
      SmFAll
   real(c_double), dimension(:, :), pointer, public :: &
      SmFCell, SmFFrac, SmFNode

   integer(c_int), allocatable, dimension(:, :), target, public :: &
      NumIncTotalPrimAll
   integer(c_int), dimension(:, :), pointer, public :: &
      NumIncTotalPrimCell, NumIncTotalPrimFrac, NumIncTotalPrimNode

   integer(c_int), allocatable, dimension(:, :), target, public :: &
      NumIncTotalSecondAll
   integer(c_int), dimension(:, :), pointer, public :: &
      NumIncTotalSecondCell, NumIncTotalSecondFrac, NumIncTotalSecondNode

   public :: &
      IncPrimSecd_allocate, &
      IncPrimSecd_free, &
      IncPrimSecd_compute, &
      IncPrimSecd_ps_cv, & ! choose first and second variables for each cv (called from IncPrimSecdFreeFlow)
      IncPrimSecd_dFsurdX_cv, & ! compute dF/dX for each cv (called from IncPrimSecdFreeFlow)
      IncPrimSecd_compPrim_nodes

   public :: &
      IncPrimSecd_compute_cv, & ! all operations for one cv
      IncPrimSecd_dXssurdXp_cv              ! compute dFs/dXp

contains

   !> \brief Main subroutine of this module,
   !! compute dXssurdXp, SmdXs, NumIncTotalPrim, NumIncTotalSecond
   !! loop of cell/frac/node
   subroutine IncPrimSecd_compute() &
      bind(C, name="IncPrimSecd_compute")

      !< cell
      call IncPrimSecd_compute_cv( &
         NbCellLocal_Ncpus(commRank + 1), &
         IncCell, dPhasePressuredSCell, &
         dXssurdXpCell, SmdXsCell, &
         SmFCell, &
         NumIncTotalPrimCell, NumIncTotalSecondCell)

      !< frac
      call IncPrimSecd_compute_cv( &
         NbFracLocal_Ncpus(commRank + 1), &
         IncFrac, dPhasePressuredSFrac, &
         dXssurdXpFrac, SmdXsFrac, &
         SmFFrac, &
         !
         NumIncTotalPrimFrac, NumIncTotalSecondFrac)

      !< node
      call IncPrimSecd_compute_cv( &
         NbNodeLocal_Ncpus(commRank + 1), &
         IncNode, dPhasePressuredSNode, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode, &
         !
         NumIncTotalPrimNode, NumIncTotalSecondNode &
#ifdef _WITH_FREEFLOW_STRUCTURES_
         , IsFreeflowNode &
#endif
         )

   end subroutine IncPrimSecd_compute

   !> \brief  all operations for a set of cv (cell/frac/node)
   subroutine IncPrimSecd_compute_cv( &
      NbIncLocal, &
      inc, dpadS, &
      dXssurdXp, SmdXs, SmF, &
      NumIncTotalPrimCV, NumIncTotalSecondCV, &
      skip_cv)
      integer, intent(in) :: NbIncLocal
      type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)
      real(c_double), intent(in) :: dpadS(NbPhase, NbIncLocal)
      integer, intent(out) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax, NbIncLocal), &
         NumIncTotalSecondCV(NbEqFermetureMax, NbIncLocal)
      real(c_double), intent(out) :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs(NbEqFermetureMax, NbIncLocal), &
         SmF(NbEqFermetureMax, NbIncLocal)
      logical(c_bool), optional, intent(in) :: skip_cv(:)

      integer :: k
      double precision :: &
         dFsurdX(NbIncTotalMax, NbEqFermetureMax) ! (col,row) index order
      type(ControlVolumeInfo) :: cv_info

#ifndef NDEBUG
      if (present(skip_cv)) then
         if (size(skip_cv) /= NbIncLocal) &
            call CommonMPI_abort("Inconsistent mask size.")
      endif
#endif
      do k = 1, NbIncLocal

#ifdef _WITH_FREEFLOW_STRUCTURES_
#ifndef NDEBUG
         if (inc(k)%ic <= 0) &
            call CommonMPI_abort("Inconsistent negative context.")
         if (present(skip_cv)) then
            if ((inc(k)%ic <= 2**NbPhase - 1) .and. skip_cv(k)) then
               write (*, *) "Context", inc(k)%ic, ": node is marked as freeflow."
               call CommonMPI_abort("Inconsistent context: node is marked as freeflow.")
            endif
            if ((inc(k)%ic > 2**NbPhase - 1) .and. .not. skip_cv(k)) then
               write (*, *) "Context", inc(k)%ic, ": node is not marked as freeflow."
               call CommonMPI_abort("Inconsistent context: node is not marked as freeflow.")
            endif
         endif
#endif
#endif
         if (present(skip_cv)) then
            if (skip_cv(k)) cycle
         endif

         ! init tmp values for each cv
         call IncPrimSecdTypes_collect_cv_info(inc(k)%ic, cv_info)

         !< compute dF/dX
         !< dFsurdX: (col, row) index order
         call IncPrimSecd_dFsurdX_cv(cv_info, inc(k), dpadS(:, k), dFsurdX, SmF(:, k))

         !< choose inconnues prim and secd
         call IncPrimSecd_ps_cv(cv_info, inc(k), dFsurdX, pschoice, &
                                NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k))

         !< compute dXssurdxp and SmdXs
         call IncPrimSecd_dXssurdXp_cv(cv_info, inc(k), dFsurdX, SmF(:, k), &
                                       NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                       dXssurdXp(:, :, k), SmdXs(:, k))

      enddo

   end subroutine IncPrimSecd_compute_cv

   !> \brief  Update prim/secd arrays of nodes
   subroutine IncPrimSecd_compPrim_nodes

      ! node
      call IncPrimSecd_compute_cv( &
         NbNodeLocal_Ncpus(commRank + 1), &
         IncNode, dPhasePressuredSNode, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode, &
         !
         NumIncTotalPrimNode, NumIncTotalSecondNode)

   end subroutine IncPrimSecd_compPrim_nodes

   !> \todo FIXME: reprendre notations de Xing et al. 2017
   !> \brief  compute dF/dX
   !! F = closure equations
   !!    * molar fractions sum to 1 for present phases
   !!    * thermodynamic equilibrium between phases if any (fugacities equality)
   !! Remark: Sg+Sl=1 is not a closure law, as well as Pphase=Pref+Pc(phase)
   !! compute dFsurdX for each control volume (Careful: the lines of the derivatives must coincide with the index of unknowns in DefModel.F90)
   !!      dFsurdX(1,:)                                            derivative reference Pressure
   !!      #ifdef _THERMIQUE_ dFsurdX(2,:)                         derivative Temperature
   !!      dFsurdX(2+IndThermique:NbEquilibre+IndThermique+1,:)    derivative Components
   !!      dFsurdX(NbIncPTC+1:NbIncPTC+NbPhasePresente+1, :)       derivative principal Saturations
   subroutine IncPrimSecd_dFsurdX_cv(cv_info, inc, dpadS, dFsurdX, SmF)
      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      real(c_double), intent(in) :: dpadS(NbPhase)
      real(c_double), intent(out) :: dFsurdX(NbIncTotalMax, NbEqFermetureMax)
      real(c_double), intent(out) :: SmF(NbEqFermetureMax)

      integer :: i, mi, iph, iph1, iph2, icp, j, jph, n, jph_n, numj
      real(c_double) :: &
         f1, df1dpa, df1dT, df1dC(NbComp), &
         f2, df2dpa, df2dT, df2dC(NbComp)

      dFsurdX = 0.d0  ! local to this file, cannot rely on previous computation in IncPrimSecd
      SmF = 0.d0

      ! --------------------------------------------------------------------------
      ! molar fractions sum to 1
      ! 1. F = sum_icp C_icp^iph(i) - 1, for i

      ! loop for rows associate with C_icp^iph(i) in dFsurdX
      do i = 1, cv_info%NbPhasePresente  ! row is i, col is j
         iph = cv_info%NumPhasePresente(i)

         do icp = 1, NbComp ! loop for cols
            if (MCP(icp, iph) == 1) then
               j = cv_info%NumIncComp2NumIncPTC(icp, iph)  ! num of C_icp^iph in IncPTC
               dFsurdX(j, i) = 1.d0  ! derivative wrt C
               SmF(i) = SmF(i) + inc%Comp(icp, iph)
            end if
         end do
         SmF(i) = SmF(i) - 1.d0
      enddo

      mi = cv_info%NbPhasePresente ! nb of closure eq already stored, mi+1 futur row of dFsurdX and SmF

      ! --------------------------------------------------------------------------
      ! thermodynamic equilibrium - fugacities equality
      ! 2. F = f_i^alpha - f_i^beta
#if defined ComPASS_SINGLE_PHASE
      ! FIXME: put this test in between NDEBUG preprocessors directives
      !        Once some CI tests are run in Debug mode
      if (cv_info%NbEqEquilibre /= 0) &
         call CommonMPI_abort("Model inconsistency...")
#else
      if (cv_info%NbEqEquilibre /= 0) then

         do i = 1, cv_info%NbEqEquilibre

            icp = cv_info%NumCompEqEquilibre(i) ! component i
            iph1 = cv_info%Num2PhasesEqEquilibre(1, i) ! phase alpha
            iph2 = cv_info%Num2PhasesEqEquilibre(2, i) ! phase beta

            ! fugacity and derivative
            call f_Fugacity(icp, iph1, inc%phase_pressure(iph1), inc%Temperature, inc%Comp(:, iph1), f1, df1dpa, df1dT, df1dC)
            call f_Fugacity(icp, iph2, inc%phase_pressure(iph2), inc%Temperature, inc%Comp(:, iph2), f2, df2dpa, df2dT, df2dC)

            ! derivative wrt reference pressure
            ! f(pa, ....) -> df/dP = dpa/dP * df/dpa = df/dpa because pa = P - Pc(S)
            dFsurdX(1, i + mi) = df1dpa - df2dpa

#ifdef _THERMIQUE_
            ! derivative wrt temperature
            dFsurdX(2, i + mi) = df1dT - df2dT
#endif

            ! derivative of df1dC wrt the components
            do j = 1, NbComp
               if (MCP(j, iph1) == 1) then ! phase iph1 contains component j
                  numj = cv_info%NumIncComp2NumIncPTC(j, iph1) ! num of C_j^iph1 in Inc
                  dFsurdX(numj, i + mi) = df1dC(j)
               end if
            end do

            ! derivative of (- df2dC ) wrt the Components
            do j = 1, NbComp
               if (MCP(j, iph2) == 1) then ! phase iph2 contains component j
                  numj = cv_info%NumIncComp2NumIncPTC(j, iph2) ! num of C_j^iph2 in Inc
                  dFsurdX(numj, i + mi) = -df2dC(j)
               end if
            end do

            ! SmF
            SmF(i + mi) = f1 - f2

            ! derivative wrt the primary saturations S_j; 1 <= j <= n-1
            ! Careful of the derivative of fn(pn, ....) wrt S_j due to
            !   the elimination of S_n using S_n = 1 - \Sum_{1 \leq j \leq n-1} S_j
            ! Remark: in (p_alpha, T_alpha, C_alpha), only p_alpha depends on S_alpha
            !   and dpadS = d(p_alpha)/d(S_alpha) with 1 <= alpha <= n
            !
            ! for 1 <= alpha <= n-1:
            !   f_alpha depends on (p_alpha, T_alpha, C_alpha) then
            !   d(f_alpha)/d(S_j) = d(p_alpha)/d(S_j) * d(f_alpha)/d(p_alpha)
            !   which is non-nul only when j=alpha (because p_alpha = p_ref + Pc(S_alpha))
            !
            ! fn(pn, Tn, Cn) -> dfn/dS_j = - dfn/dS_n = -dp_n/dS_n * dfn/dp_n
            n = cv_info%NbPhasePresente
            jph_n = cv_info%NumPhasePresente(n) ! phase of the eliminated saturation
            do j = 1, n - 1
               numj = j + cv_info%NbIncPTC
               jph = cv_info%NumPhasePresente(j)
#ifndef NDEBUG
               if (jph == jph_n) call CommonMPI_abort("IncPrimSecd_dFsurdX_cv inconsistent phase indexing.")
#endif
               ! d(f_iph1)/d(S_j) careful if iph1 = n (see above)
               if (iph1 == jph) dFsurdX(numj, i + mi) = dFsurdX(numj, i + mi) + dpadS(iph1)*df1dpa
               if (iph1 == jph_n) dFsurdX(numj, i + mi) = dFsurdX(numj, i + mi) - dpadS(iph1)*df1dpa
               ! - d(f_iph2)/d(S_j) careful if iph2 = n (see above)
               if (iph2 == jph) dFsurdX(numj, i + mi) = dFsurdX(numj, i + mi) - dpadS(iph2)*df2dpa
               if (iph2 == jph_n) dFsurdX(numj, i + mi) = dFsurdX(numj, i + mi) + dpadS(iph2)*df2dpa
            end do

         end do ! component i

      end if
#endif

   end subroutine IncPrimSecd_dFsurdX_cv

   !> \brief choose primary and secondary unknowns for each CV
   !! fill inc%Nb/NumIncTotalPrim/Secd
   subroutine IncPrimSecd_ps_cv(cv_info, inc, dFsurdX, pschoicecv, &
                                NumIncTotalPrimCV, NumIncTotalSecondCV)
      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      ! dFsurdX may be used depending on choice method
      ! (e.g. projection on closure equations)
      ! ??? pourrait-on faire mieux qu'une extraction des composantes ???
      double precision, intent(in) :: dFsurdX(NbIncTotalMax, NbEqFermetureMax)
      integer, intent(in) :: pschoicecv
      integer, intent(out) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(out) :: NumIncTotalSecondCV(NbEqFermetureMax)

      integer :: i, ic

      NumIncTotalPrimCV(:) = 0
      NumIncTotalSecondCV(:) = 0

      if (pschoicecv == 1) then ! manually

         ic = inc%ic

         ! prim variable
         do i = 1, NbIncTotalPrim_ctx(ic)
            NumIncTotalPrimCV(i) = psprim(i, ic)
         end do

         ! secd variable
         do i = 1, NbEqFermeture_ctx(ic)
            NumIncTotalSecondCV(i) = pssecd(i, ic)
         end do

      else if (pschoicecv == 2) then ! Glouton method
         ! call IncPrimSecd_IncSecondGluton(inc, dFsurdX, &
         ! NumIncTotalPrimCV, NumIncTotalSecondCV)
         call CommonMPI_abort("Glouton method is not implemented in IncPrimSecd")
      else if (pschoicecv == 3) then ! Gauss method
         call IncPrimSecd_IncSecondGauss(cv_info, inc, dFsurdX, &
                                         NumIncTotalPrimCV, NumIncTotalSecondCV)
      end if

   end subroutine IncPrimSecd_ps_cv

   !> \brief Compute dXssurdXp = dFsurdXs**(-1) * dFsurdXp
   !! and SmdXs = dFsurdXs**(-1) * SmF
   subroutine IncPrimSecd_dXssurdXp_cv(cv_info, inc, dFsurdX, SmF, &
                                       NumIncTotalPrimCV, NumIncTotalSecondCV, dXssurdXp, SmdXs)

      ! inputs
      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      double precision, intent(in) :: & ! (col, row) index order
         dFsurdX(NbIncTotalMax, NbEqFermetureMax), &
         SmF(NbEqFermetureMax)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax), &
         NumIncTotalSecondCV(NbEqFermetureMax)

      ! output
      double precision, intent(out) :: & ! (col, row) index order
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax), &
         SmdXs(NbEqFermetureMax)

      ! tmp
      double precision :: & ! (row,col) index order, lapack
         dFsurdX_prim(NbEqFermetureMax, NbIncTotalPrimMax), &
         dFsurdX_secd(NbEqFermetureMax, NbEqFermetureMax)

      ! parameters for lapack
      integer :: ipiv(NbEqFermetureMax), blas_info, Ierr, errcode
      integer :: i, j, iph

      ! from dFsurdX, take out the cols of prim and secd variable
      do j = 1, cv_info%NbIncTotalPrim
         do i = 1, cv_info%NbEqFermeture
            dFsurdX_prim(i, j) = dFsurdX(NumIncTotalPrimCV(j), i)
         enddo
      enddo

      do j = 1, cv_info%NbEqFermeture ! = Nb of secd unknowns
         do i = 1, cv_info%NbEqFermeture
            dFsurdX_secd(i, j) = dFsurdX(NumIncTotalSecondCV(j), i)
         enddo
      enddo

      ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
      call dgetrf(cv_info%NbEqFermeture, cv_info%NbEqFermeture, &
                  dFsurdX_secd, NbEqFermetureMax, ipiv, blas_info)

      if (blas_info /= 0) then
         write (0, '(A,I0)') "dgetrf error in dXssurdxp, info = ", blas_info

         write (*, *) ' inc ic ', inc%ic
         write (*, *) ' inc P ', inc%Pression
         write (*, *) ' inc T ', inc%Temperature
         write (*, *) ' Sat ', inc%Saturation
         DO i = 1, cv_info%NbPhasePresente
            iph = cv_info%NumPhasePresente(i)
            write (*, *) ' phase  ', iph
            write (*, *) ' C_i  ', inc%Comp(:, iph)
         ENDDO
         write (*, *)

         do i = 1, cv_info%NbEqFermeture ! = Nb of secd unknowns
            do j = 1, cv_info%NbEqFermeture
               write (*, *) ' dFsurdXs ', i, j, dFsurdX(NumIncTotalSecondCV(j), i)
            enddo
            write (*, *)
         enddo

         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      call dgetrs('N', cv_info%NbEqFermeture, cv_info%NbIncTotalPrim, &
                  dFsurdX_secd, NbEqFermetureMax, &
                  ipiv, dFsurdX_prim, NbEqFermetureMax, blas_info)
      if (blas_info /= 0) then
         write (0, '(A,I0)') "dgetrs error in dXssurdxp, info = ", blas_info
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      do j = 1, cv_info%NbIncTotalPrim
         do i = 1, cv_info%NbEqFermeture
            dXssurdXp(j, i) = dFsurdX_prim(i, j)
         enddo
      enddo

      ! SmdXs = dFsurdXs**(-1) * SmF
      call dgetrs('N', cv_info%NbEqFermeture, 1, &
                  dFsurdX_secd, NbEqFermetureMax, &
                  ipiv, SmF, NbEqFermetureMax, blas_info)
      if (blas_info /= 0) then
         write (0, '(A,I0)') "dgetrs error in SmdXs, info = ", blas_info
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      do i = 1, cv_info%NbEqFermeture
         SmdXs(i) = SmF(i)
      end do

   end subroutine IncPrimSecd_dXssurdXp_cv

   !> \brief  Choice of primary and secondary unknowns
   !! from the matrix dFsurdX with the Gauss algorithm
   subroutine IncPrimSecd_IncSecondGauss(cv_info, inc, dFsurdX, &
                                         NumIncTotalPrimCV, NumIncTotalSecondCV)

      ! input
      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      double precision, intent(in) :: dFsurdX(NbIncTotalMax, NbEqFermetureMax)

      ! output
      integer, intent(out) :: NumIncTotalPrimCV(NbIncTotalPrimMax)
      integer, intent(out) :: NumIncTotalSecondCV(NbEqFermetureMax)

      ! tmps
      double precision :: &
         BB(NbIncTotalMax, NbEqFermetureMax) ! copy of dFsurdX

      double precision :: pivotmax
      integer :: npivot, pivot(2, NbEqFermetureMax*NbIncTotalMax) !> \todo FIXME: check the size
      logical :: pivot_P, pivot_T

      integer :: & ! E^{eq}, E^{inc}
         NbSetInc, &
         NbSetEq, &
         NumSetInc(NbIncTotalMax), &
         NumSetEq(NbEqFermetureMax)

      logical :: &
         NumIncTotalPrim_idx(NbIncTotalMax)

      integer :: is, i, j, numi, numj, n
      integer :: i1, j1, icp, iph, icp1, iph1, k

      ! 1. NumIncTotalSecondcv
      !    1.1 choose secd in (P,T,C), Gauss
      !    1.2 choose secd in S, the first is secd
      ! 2. NumIncTotalPrimcv

      ! 1.1 choose secd in (P,T,C), Gauss

      ! if(commRank==0) then
      !    do i=1, NbEqFermeture
      !       do j=1, NbIncPTC
      !          print*, dFsurdX(j,i)
      !       end do
      !       print*, ""
      !    end do
      ! end if

      ! init set of inconnus and equations
      NbSetInc = cv_info%NbIncPTC
      NbSetEq = cv_info%NbEqFermeture

      do j = 1, NbSetInc
         NumSetInc(j) = j
      end do

      do i = 1, NbSetEq
         NumSetEq(i) = i
      end do

      BB(:, :) = dFsurdX(:, :) ! (col, row) index order

      do is = 1, cv_info%NbEqFermeture ! = nb of secd inconnues

         ! max element of abs(BB(:,:))
         pivotmax = -1.d0
         do numi = 1, NbSetEq
            i = NumSetEq(numi)

            do numj = 1, NbSetInc
               j = NumSetInc(numj)

               if (abs(BB(j, i)) > pivotmax) then
                  pivotmax = abs(BB(j, i))
               end if
            end do
         end do

         ! set of element (BB) that takes the maximum value
         npivot = 0
         pivot_P = .false.
         pivot_T = .false.

         do numi = 1, NbSetEq
            i = NumSetEq(numi)

            do numj = 1, NbSetInc
               j = NumSetInc(numj)

               if (abs(abs(BB(j, i)) - pivotmax) < eps) then

                  npivot = npivot + 1

                  pivot(1, npivot) = i ! num eq in BB, not in subset of BB
                  pivot(2, npivot) = j ! num inc in BB, not in subset of BB

                  ! check if j is P or T
                  if (j == 1) then
                     pivot_P = .true.
                  else if (j == 2) then
                     pivot_T = .true.
                  end if
               end if

            end do
         end do

         ! if(commRank==1 .and. is==2) then
         !    print*, pivotmax, npivot
         !    do j=1, npivot
         !       print*, pivot(1,j), pivot(2,j)
         !    end do
         !    !print*, ""
         ! end if

         ! choose (i1,j1) is in set pivot(:,),
         ! (i1,j1) is num of BB, not subset of BB

         ! sinon on prend T si elle est dans l'ensemble des pivots max
         if (pivot_T .eqv. .true.) then

            do k = 1, npivot
               if (pivot(2, k) == 2) then
                  i1 = pivot(1, k)
                  j1 = 2
               end if
            end do

            NumIncTotalSecondCV(is) = 2 ! j1=2
         else
            ! sinon on prend la plus grande composition dans la phase
            ! avec la plus petite saturation

            i1 = pivot(1, 1) ! num of Eq
            j1 = pivot(2, 1) ! num of Inc
            icp1 = cv_info%NumIncPTC2NumIncComp_comp(j1)
            iph1 = cv_info%NumIncPTC2NumIncComp_phase(j1)

            ! if(is==2 .and. commRank==1) then
            !    print*,"init", i1, j1
            ! end if

            do k = 2, npivot

               i = pivot(1, k)
               j = pivot(2, k)
               icp = cv_info%NumIncPTC2NumIncComp_comp(j)  ! icp in C_{icp}^iph
               iph = cv_info%NumIncPTC2NumIncComp_phase(j) ! iph in C_{icp}^iph

               !              if(commRank==1 .and. is==2) then
               ! !                print*, i1,j1
               !                 print*, & !i, j ,inc%Saturation(iph), inc%Saturation(iph1), &
               !                      icp, iph, j!, inc%Comp(icp,iph), inc%Comp(icp1, iph1)
               !                 print*, ""
               !              end if

               ! update (i1, j1)
               if (inc%Saturation(iph) < inc%Saturation(iph1)) then ! if < (Saturation)
                  i1 = i
                  j1 = j
                  icp1 = icp
                  iph1 = iph

               else if ((abs(inc%Saturation(iph) - inc%Saturation(iph1)) < eps)) then ! if = (Saturation) and ...
                  if (inc%Comp(icp, iph) > inc%Comp(icp1, iph1)) then !
                     i1 = i
                     j1 = j
                     icp1 = icp
                     iph1 = iph
                  end if
               end if

            end do

            NumIncTotalSecondCV(is) = j1
         end if ! end for choosing (i1, j1), NumIncTotalSecondcv(is)=j1

         ! update sets: NumSetInc and NumSetEq, remove (i1, j1)
         n = 0
         do i = 1, NbSetEq
            numi = NumSetEq(i)

            if (numi .ne. i1) then
               n = n + 1
               NumSetEq(n) = numi
            end if
         end do
         NbSetEq = n

         n = 0
         do j = 1, NbSetInc
            numj = NumSetInc(j)

            if (numj .ne. j1) then
               n = n + 1
               NumSetInc(n) = numj
            end if
         end do
         NbSetInc = n

         ! schur complement
         do numi = 1, NbSetEq
            i = NumSetEq(numi)

            do numj = 1, NbSetInc
               j = NumSetInc(numj)

               !print*, i1, j1, numi, numj

               BB(j, i) = BB(j, i) &
                          - BB(j1, i)*BB(j, i1)/BB(j1, i1)
            end do
         end do

      end do ! end loop of is

      ! 2. NumIncTotalPrimcv = {1,2,...,NbIncPTC}/NumIncTotalSecondcv
      !                       + S^alpha, alpha=1:Nphase-1
      NumIncTotalPrim_idx(:) = .true.
      do j = 1, cv_info%NbEqFermeture
         NumIncTotalPrim_idx(NumIncTotalSecondCV(j)) = .false. ! not prim
      end do

      n = 0
      do j = 1, cv_info%NbIncPTC
         if (NumIncTotalPrim_idx(j) .eqv. .true.) then ! prim
            n = n + 1
            NumIncTotalPrimCV(n) = j
         end if
      end do

      ! last S is secd
      ! if there is only one phase, phase is secd
      do j = cv_info%NbIncPTC + 1, cv_info%NbIncTotal - 1  ! NbIncTotal is not possible here, if there are more unknowns than P T C S
         n = n + 1
         NumIncTotalPrimCV(n) = j
      end do

   end subroutine IncPrimSecd_IncSecondGauss

   !> \brief  Link with the C code
   subroutine IncPrimSecd_PrimToSecd_C(increments_pointers) &
      bind(C, name="IncPrimSecd_PrimToSecd")

      type(Newton_increments_pointers), intent(in), value :: increments_pointers
      type(Newton_increments) :: increments

      call Newton_pointers_to_values(increments_pointers, increments)
      call IncPrimSecd_PrimToSecd( &
         increments%nodes, increments%fractures, increments%cells &
         )

   end subroutine IncPrimSecd_PrimToSecd_C

   !> \brief Compute secd values using prim values
   !! secd = SmdX - dXssurdXp * prim
   !! for node and frac and cell
   ! v for variation ?
   subroutine IncPrimSecd_PrimToSecd( &
      vnode, vfrac, vcell)

      real(c_double), dimension(:, :), intent(inout) :: &
         vnode, vfrac, vcell

      ! node
      call IncPrimSecd_PrimToSecd_cv( &
         IncNode, NbNodeLocal_Ncpus(commRank + 1), &
         dXssurdXpNode, SmdXsNode, &
         NumIncTotalPrimNode, NumIncTotalSecondNode, &
         vnode)

      ! frac
      call IncPrimSecd_PrimToSecd_cv( &
         IncFrac, NbFracLocal_Ncpus(commRank + 1), &
         dXssurdXpFrac, SmdXsFrac, &
         NumIncTotalPrimFrac, NumIncTotalSecondFrac, &
         vfrac)

      ! cell
      call IncPrimSecd_PrimToSecd_cv( &
         IncCell, NbCellLocal_Ncpus(commRank + 1), &
         dXssurdXpCell, SmdXsCell, &
         NumIncTotalPrimCell, NumIncTotalSecondCell, &
         vcell)

   end subroutine IncPrimSecd_PrimToSecd

   !> \brief Compute secd values using prim values
   !! for a set of cv (node or frac or cell)
   subroutine IncPrimSecd_PrimToSecd_cv( &
      inc, NbIncLocal, &
      dXssurdXp, SmdXs, &
      NumIncTotalPrimCV, NumIncTotalSecondCV, &
      var_inc)

      ! input
      type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)

      integer, intent(in) :: NbIncLocal

      double precision, intent(in) :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs(NbEqFermetureMax, NbIncLocal)

      integer, intent(in) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax, NbIncLocal), &
         NumIncTotalSecondCV(NbEqFermetureMax, NbIncLocal)

      ! inout: need dvar_prim to fill dvar (both prim and secd)
      real(c_double), dimension(:, :), intent(inout) :: var_inc

      integer :: k, ic, i, iph
      integer :: &
         NbPhasePresente, NbEqFermeture, NbNodeLocal, &
         NbIncPTC, NbIncTotal, NbIncPTCPrim, NbIncTotalPrim

      double precision :: &
         xp(NbCompThermique), &
         xs(NbEqFermetureMax)

      do k = 1, NbIncLocal

         ic = inc(k)%ic
         NbNodeLocal = NbNodeLocal_Ncpus(commRank + 1)
         NbPhasePresente = NbPhasePresente_ctx(ic)
         NbEqFermeture = NbEqFermeture_ctx(ic)
         NbIncPTC = NbIncPTC_ctx(ic)
         NbIncTotal = NbIncTotal_ctx(ic)
         NbIncTotalPrim = NbIncTotalPrim_ctx(ic)
         NbIncPTCPrim = NbIncPTC - NbEqFermeture

         xp(1:NbCompThermique) = var_inc(1:NbCompThermique, k)
         xs(1:NbEqFermeture) = SmdXs(1:NbEqFermeture, k)

         ! http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
         ! y := alpha*A*x + beta*y
         ! Compute xs = xs - dXssurdXp * xp, where xs=SmdXs and xp=var_prim
         call dgemv('T', NbIncTotalPrim, NbEqFermeture, &
                    -1.d0, dXssurdXp(:, :, k), NbIncTotalPrimMax, &
                    xp(:), 1, -1.d0, xs(:), 1)

         !-----------------------------------------------------
         ! update var_inc with the primary/secondary increments
         !-----------------------------------------------------

         var_inc(:, k) = 0.d0

         ! copy prim P,T,C,S
         do i = 1, NbIncTotalPrim
            var_inc(NumIncTotalPrimCV(i, k), k) = xp(i)
         end do

         ! copy secd P,T,C
         do i = 1, NbEqFermeture
            var_inc(NumIncTotalSecondCV(i, k), k) = xs(i)
         end do

         ! fill secd S
         ! if NbPhasePresente=1,
         !    then this phase is eliminated
         ! else last saturation is eliminated, the others are prim (in reservoir dof)
         !    eliminated S = - sum_{S^alpha is prim} S^alpha
         ! iph is last present phase in vector (P,T,C,S,n_i)
         iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente, ic)
         var_inc(iph, k) = 0.d0 !> \todo FIXME: is useful, because wrong numerotation of FreeFlow sat
         do i = 1, NbPhasePresente - 1
            var_inc(iph, k) = var_inc(iph, k) - xp(NbIncPTCPrim + i)
         end do

#ifdef _WITH_FREEFLOW_STRUCTURES_
         ! FIXME: CRITICAL: we should check against a freeflow tag
         if (ic > 2**NbPhase - 1 .and. NbPhasePresente > 1) then ! FIXME: loop over freeflow dof only, avoid reservoir node
            ! Correct the value for the last saturation
            ! Remark : it is possible that no saturation belongs to the unknowns (when liquid outflow)
            iph = NbIncPTC + NumPhasePresente_ctx(NbPhasePresente, ic)
            var_inc(iph, k) = 0.d0
            ! look for saturations which are primary unknowns
            !    eliminated S = - sum_{S^alpha is prim} S^alpha
            do i = NbIncPTCPrim + 1, NbIncPTCPrim + NbPhasePresente - 1
               if (NumIncTotalPrimCV(i, k) > NbIncPTC .and. NumIncTotalPrimCV(i, k) < NbIncPTC + NbPhasePresente) then ! test if it is a saturation
                  var_inc(iph, k) = var_inc(iph, k) - xp(i)
               endif
            enddo
         endif

         if (ic >= 2**NbPhase) then ! FIXME: loop over freeflow dof only, avoid reservoir node
            !> \todo FIXME: change the numbering in DefModel to remove this part
            !! Correct the value of var_inc because the eliminated saturation is inserted in the index
            !! whereas it is not counted in the numbering in DefModel
            !! necessary to coincide with IncCVReservoir_NewtonIncrement_reservoir
            !! copy prim
            do i = 1, NbIncTotalPrim
               if (NumIncTotalPrimCV(i, k) >= NbIncPTC + NbPhasePresente) var_inc(NumIncTotalPrimCV(i, k) + 1, k) = xp(i)
            end do

            ! copy secd
            do i = 1, NbEqFermeture
               if (NumIncTotalSecondCV(i, k) >= NbIncPTC + NbPhasePresente) var_inc(NumIncTotalSecondCV(i, k) + 1, k) = xs(i)
            end do
         endif

#endif

         ! term prim n_k(X_j^n), components which are present only in absent phase(s)
         ! n_k(X_j^n) are not part of NbIncTotal
         ! copy prim n_i
         do i = 1, NbCompCtilde_ctx(ic) ! =NbCompThermique-NbIncTotalPrim
            var_inc(NbIncTotal + i, k) = xp(NbIncTotalPrim + i)
         end do

      end do ! loop over local inc

   end subroutine IncPrimSecd_PrimToSecd_cv

   ! allocate
   subroutine IncPrimSecd_allocate

      integer :: nbCell, nbFrac, nbNode, n
      integer :: begin_node, end_node
      integer :: begin_frac, end_frac
      integer :: begin_cell, end_cell

      nbCell = NbCellLocal_Ncpus(commRank + 1)
      nbFrac = NbFracLocal_Ncpus(commRank + 1)
      nbNode = NbNodeLocal_Ncpus(commRank + 1)
      n = nbNode + nbFrac + nbCell
      begin_node = 1
      end_node = nbNode
      begin_frac = end_node + 1
      end_frac = end_node + nbFrac
      begin_cell = end_frac + 1
      end_cell = end_frac + nbCell

      allocate (dXssurdXpAll(NbIncTotalPrimMax, NbEqFermetureMax, n))
      dXssurdXpNode => dXssurdXpAll(:, :, begin_node:end_node)
      dXssurdXpFrac => dXssurdXpAll(:, :, begin_frac:end_frac)
      dXssurdXpCell => dXssurdXpAll(:, :, begin_cell:end_cell)

      allocate (SmdXsAll(NbEqFermetureMax, n))
      SmdXsNode => SmdXsAll(:, begin_node:end_node)
      SmdXsFrac => SmdXsAll(:, begin_frac:end_frac)
      SmdXsCell => SmdXsAll(:, begin_cell:end_cell)

      allocate (SmFAll(NbEqFermetureMax, n))
      SmFNode => SmFAll(:, begin_node:end_node)
      SmFFrac => SmFAll(:, begin_frac:end_frac)
      SmFCell => SmFAll(:, begin_cell:end_cell)

      allocate (NumIncTotalPrimAll(NbIncTotalPrimMax, n))
      NumIncTotalPrimNode => NumIncTotalPrimAll(:, begin_node:end_node)
      NumIncTotalPrimFrac => NumIncTotalPrimAll(:, begin_frac:end_frac)
      NumIncTotalPrimCell => NumIncTotalPrimAll(:, begin_cell:end_cell)

      allocate (NumIncTotalSecondAll(NbEqFermetureMax, n))
      NumIncTotalSecondNode => NumIncTotalSecondAll(:, begin_node:end_node)
      NumIncTotalSecondFrac => NumIncTotalSecondAll(:, begin_frac:end_frac)
      NumIncTotalSecondCell => NumIncTotalSecondAll(:, begin_cell:end_cell)

   end subroutine IncPrimSecd_allocate

   ! free
   subroutine IncPrimSecd_free

      nullify (dXssurdXpCell)
      nullify (dXssurdXpFrac)
      nullify (dXssurdXpNode)
      deallocate (dXssurdXpAll)

      nullify (SmdXsCell)
      nullify (SmdXsFrac)
      nullify (SmdXsNode)
      deallocate (SmdXsall)

      nullify (SmFCell)
      nullify (SmFFrac)
      nullify (SmFNode)
      deallocate (SmFAll)

      nullify (NumIncTotalPrimCell)
      nullify (NumIncTotalPrimFrac)
      nullify (NumIncTotalPrimNode)
      deallocate (NumIncTotalPrimAll)

      nullify (NumIncTotalSecondCell)
      nullify (NumIncTotalSecondFrac)
      nullify (NumIncTotalSecondNode)
      deallocate (NumIncTotalSecondAll)

   end subroutine IncPrimSecd_free

end module IncPrimSecd
