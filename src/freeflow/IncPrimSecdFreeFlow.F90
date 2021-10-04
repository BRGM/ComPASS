!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module IncPrimSecdFreeFlow

   ! FIXME: add only: check that all members are necessary
   use iso_c_binding, only: c_int, c_double
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort
   use DefModel, only: &
      NbComp, NbPhase, MCP, GAS_PHASE, &
      NbIncTotalPrimMax, NbEqFermetureMax, NbIncTotalMax, &
      pschoice, psprim, pssecd, &
      NumPhasePresente_ctx, NbPhasePresente_ctx, &
      NbIncPTCMax, IndThermique, NbCompThermique, NbEqEquilibreMax
   use NumbyContext, only: &
      NbEqFermeture_ctx, NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, &
      NumIncComp2NumIncPTC_ctx, NbIncTotalPrim_ctx, NumCompCtilde_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
      NbEqEquilibre_ctx, NbIncPTC_ctx, NbIncTotal_ctx, NbCompCtilde_ctx
   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, &
      IncNode, dPhasePressuredSNode
   use IncPrimSecd, only: IncPrimSecd_ps_cv, IncPrimSecd_dFsurdX_cv, &
                          dXssurdXpCell, dXssurdXpFrac, dXssurdXpNode, &
                          SmdXsCell, SmdXsFrac, SmdXsNode, &
                          SmFCell, SmFFrac, SmFNode, &
                          NumIncTotalPrimCell, NumIncTotalPrimFrac, NumIncTotalPrimNode, &
                          NumIncTotalSecondCell, NumIncTotalSecondFrac, NumIncTotalSecondNode
   use Physics, only: atm_pressure
   use MeshSchema, only: &
#ifdef _WIP_FREEFLOW_STRUCTURES_
      IdFFNodeLocal, &
#endif
      NbCellLocal_Ncpus, NbFracLocal_Ncpus, NbNodeLocal_Ncpus, &
      NodeByWellInjLocal
   use IncPrimSecdTypes, only: ControlVolumeInfo, IncPrimSecdTypes_collect_cv_info

   implicit none

   public :: &
      IncPrimSecdFreeFlow_compute !, &
   !   IncPrimSecdFreeFlow_compPrim_nodes

   private :: &
      IncPrimSecdFreeFlow_compute_cv, & ! all operations for one cv (called with nodes only)
      IncPrimSecdFreeFlow_dFsurdX_cv, & ! compute dF/dX for each cv (called with nodes only)
      IncPrimSecdFreeFlow_dXssurdXp_cv  ! compute dFs/dXp (called with nodes only)

contains

   ! main subroutine of this module
   ! NumIncTotalPrim, NumIncTotalSecond filled in IncPrimSecd
   ! compute dXssurdXp, SmdXs
   ! of the FreeFlow d.o.f. (nodes only)
   subroutine IncPrimSecdFreeFlow_compute() &
      bind(C, name="IncPrimSecdFreeFlow_compute")

      ! FreeFlow BC nodes
      call IncPrimSecdFreeFlow_compute_cv( &
         NbNodeLocal_Ncpus(commRank + 1), &
         IncNode, dPhasePressuredSNode, &
         dXssurdXpNode, SmdXsNode, &
         SmFNode, &
         NumIncTotalPrimNode, NumIncTotalSecondNode)

   end subroutine IncPrimSecdFreeFlow_compute

   ! all operations for a set of cv (called with nodes only)
   subroutine IncPrimSecdFreeFlow_compute_cv( &
      NbIncLocal, &
      inc, dpadS, &
      dXssurdXp, SmdXs, SmF, &
      NumIncTotalPrimCV, NumIncTotalSecondCV)
      integer, intent(in) :: NbIncLocal
      type(TYPE_IncCVReservoir), intent(in) :: inc(NbIncLocal)
      real(c_double), intent(in) :: pa(NbPhase, NbIncLocal)
      real(c_double), intent(in) :: dpadS(NbPhase, NbIncLocal)
      integer, intent(out) :: &
         NumIncTotalPrimCV(NbIncTotalPrimMax, NbIncLocal), &
         NumIncTotalSecondCV(NbEqFermetureMax, NbIncLocal)
      real(c_double), intent(out) :: &
         dXssurdXp(NbIncTotalPrimMax, NbEqFermetureMax, NbIncLocal), & ! (col,row) index order
         SmdXs(NbEqFermetureMax, NbIncLocal), &
         SmF(NbEqFermetureMax, NbIncLocal)

      integer :: k
      double precision :: &
         dFsurdX(NbIncTotalMax, NbEqFermetureMax) ! (col,row) index order
      type(ControlVolumeInfo) :: cv_info

      do k = 1, NbIncLocal

         ! done only for ff dof
         if (.not. IdFFNodeLocal(k)) cycle ! loop over freeflow dof only, avoid reservoir node

         ! init tmp values for each cv
         call IncPrimSecdTypes_collect_cv_info(inc(k)%ic, cv_info)

         ! compute dF/dX
         ! dFsurdX: (col, row) index order
         call IncPrimSecdFreeFlow_dFsurdX_cv(cv_info, inc(k), dpadS(:, k), dFsurdX, SmF(:, k))

         ! FIXME: si je peux faire IncPrimSecd_compute sur les dof FF, enlever ce call et remettre protected à NumIncTotal...
         ! (mais cela m'étonnerait car le numb d'eq de fermeture ne correspond pas au nb d'inconnus secds dans ce cas)
         ! choose prim and sced unknowns
         call IncPrimSecd_ps_cv(cv_info, inc(k), dFsurdX, pschoice, &
                                NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k))

         ! compute dXssurdxp
         call IncPrimSecdFreeFlow_dXssurdXp_cv(cv_info, inc(k), dFsurdX, SmF(:, k), &
                                               NumIncTotalPrimCV(:, k), NumIncTotalSecondCV(:, k), &
                                               dXssurdXp(:, :, k), SmdXs(:, k))

      enddo

   end subroutine IncPrimSecdFreeFlow_compute_cv

!    !> Update prim/secd arrays of nodes (same as IncPrimSecdFreeFlow_compute)
!    ! FIXME NumIncTotalPrimNode, NumIncTotalSecondNode are updated in IncPrimSecd.F90, is it necessary to update dXssurdXpNode, SmdXsNode ?
!    subroutine IncPrimSecdFreeFlow_compPrim_nodes

!       ! node
!       call IncPrimSecdFreeFlow_compute_cv( &
!          NbNodeLocal_Ncpus(commRank + 1), &
!          IncNode, NodeDarcyRocktypesLocal, &
!          dXssurdXpNode, SmdXsNode, &
!          SmFNode, &
!          !
!          NumIncTotalPrimNode, NumIncTotalSecondNode)

!    end subroutine IncPrimSecdFreeFlow_compPrim_nodes

   ! FIXME: reprendre notations de Xing et al. 2017
   ! F = closure equations
   !    * molar fractions sum to 1 for present phases
   !    * thermodynamic equilibrium between phases if any (fugacities equality)
   !    * P^g = P^atm
   ! Remark: Sg+Sl=1 is not a closure law, as well as Pphase=Pref+Pc(phase) and q^l_atm * S^g = 0, it is forced in the implementation
   ! compute dFsurdX for each control volume (Careful: the lines of the derivatives must coincide with the index of unknowns in DefModel.F90)
   !      dFsurdX(1,:)                                            derivative Pressure
   !      #ifdef _THERMIQUE_ dFsurdX(2,:)                         derivative Temperature
   !      dFsurdX(2+IndThermique:NbEquilibre+IndThermique+1,:)    derivative Components
   !      dFsurdX(NbIncPTC+1:NbIncPTC+NbPhasePresente+1, :)       derivative principal Saturations
   !      dFsurdX(, :)                                            derivative freeflow flowrate(s)
   subroutine IncPrimSecdFreeFlow_dFsurdX_cv(cv_info, inc, dpadS, dFsurdX, SmF)
      type(ControlVolumeInfo), intent(in) :: cv_info
      type(TYPE_IncCVReservoir), intent(in) :: inc
      real(c_double), intent(in) :: dpadS(NbPhase)
      real(c_double), intent(out) :: dFsurdX(NbIncTotalMax, NbEqFermetureMax)
      real(c_double), intent(out) :: SmF(NbEqFermetureMax)

      integer :: mi, j, jph, n, jph_n, numj

      call IncPrimSecd_dFsurdX_cv(cv_info, inc, dpadS, dFsurdX, SmF)

      ! --------------------------------------------------------------------------
      ! P^g - P^atm = 0     ie      Pref + Pc(GAS_PHASE) - P^atm = 0

      mi = cv_info%NbPhasePresente + cv_info%NbEqEquilibre + 1 ! row of dFsurdX and SmF

      ! derivative wrt reference pressure
      ! pa = Pref - Pc(S) -> dpa/dS = -dPc/dS
      ! f(pa, ....) -> df/dP = dpa/dP * df/dpa = df/dpa
      dFsurdX(1, mi) = 1.d0

      ! derivative wrt primary saturations with contribution of secondary saturation
      ! because sum(saturations)=1 is eliminated
      ! pa = Pref - Pc(S) -> dpa/dS = -dPc/dS
      ! f(pa, ....) -> df/dS = -dPc/dS * df/dpa = dpa/dS * df/dpa
      ! for last saturation S_n = 1 - \Sum_{1 \leq k \leq n-1} S_k
      ! pa_n = Pref - Pc_n(1 - \Sum_{1 \leq k \leq n-1} S_k) -> dpa_n / dS_k = dPc_n/dS_n
      n = cv_info%NbPhasePresente
      jph_n = cv_info%NumPhasePresente(n)
      do j = 1, n - 1
         numj = j + cv_info%NbIncPTC
         jph = cv_info%NumPhasePresente(j)
#ifndef NDEBUG
         if (jph == jph_n) call CommonMPI_abort("IncPrimSecdFreeflow_dFsurdX_cv inconsistent phase indexing.")
#endif
         if (GAS_PHASE == jph) dFsurdX(numj, mi) = dpadS(GAS_PHASE)
         if (GAS_PHASE == jph_scd) dFsurdX(numj, mi) = -dpadS(GAS_PHASE)
      end do

      SmF(mi) = pa(GAS_PHASE) - atm_pressure

   end subroutine IncPrimSecdFreeFlow_dFsurdX_cv

   ! dXssurdXp = dFsurdXs**(-1) * dFsurdXp
   subroutine IncPrimSecdFreeFlow_dXssurdXp_cv(cv_info, inc, dFsurdX, SmF, &
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
      end if ! info

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

   end subroutine IncPrimSecdFreeFlow_dXssurdXp_cv

end module IncPrimSecdFreeFlow
