!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

#include "gravity_average_rho.def"

module Flux

   use CommonMPI, only: commRank
   use DefModel, only: &
      NbPhase, NbComp, NumPhasePresente_ctx, NbPhasePresente_ctx, phase_can_be_present

   use IncCVReservoir, only: &
      IncNode, IncCell, IncFrac, &
      NbCellLocal_Ncpus, NbFracLocal_Ncpus

   use IncCVReservoirTypes, only: TYPE_IncCVReservoir

   use LoisThermoHydro, only: &
      DensiteMassiqueNode, DensiteMassiqueCell, DensiteMassiqueFrac
   use VAGFrac, only: &
      TkLocal_Darcy, TkFracLocal_Darcy, TkLocal_Fourier, TkFracLocal_Fourier
   use Physics, only: gravity

   use MeshSchema, only: &
      NbFracCellMax, NbNodeCellMax, NbNodeFaceMax, &
      NodebyCellLocal, FracbyCellLocal, NodebyFaceLocal, &
      FracToFaceLocal, FaceToFracLocal, XNodeLocal, XFaceLocal, XCellLocal

   use Thermodynamics, only: f_DensiteMassique

   implicit none

   ! flux Darcy V_{k,s}^alpha, k is cell/frac
   double precision, allocatable, dimension(:, :, :), public :: &
      FluxDarcyKI, & !< Darcy flux from cell K to dof I (may be node or frac)
      FluxDarcyFI    !< Darcy flux from frac F to dof I (may be cell or frac)

   ! flux Fourier
   double precision, allocatable, dimension(:, :), protected :: &
      FluxFourierKI, &  !< Fourier flux from cell K to dof I (may be node or frac)
      FluxFourierFI     !< Fourier flux from frac F to dof I (may be cell or frac)

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY
   double precision, public, parameter :: epsilon_avrho = 1.d-6
#endif

   public :: &
      Flux_allocate, &
      Flux_free, &
      Flux_DarcyFlux_Cell, &
      Flux_DarcyFlux_Frac, &
      Flux_FourierFlux_Cell, &
      Flux_FourierFlux_Frac

contains

   !> \brief Allocate FluxDarcy and FluxFourier
   subroutine Flux_allocate

      ! flux
      allocate (FluxDarcyKI &
                (NbPhase, NbNodeCellMax + NbFracCellMax, NbCellLocal_Ncpus(commRank + 1)))
      allocate (FluxDarcyFI &
                (NbPhase, NbNodeFaceMax, NbFracLocal_Ncpus(commRank + 1)))

#ifdef _THERMIQUE_

      allocate (FluxFourierKI &
                (NbNodeCellMax + NbFracCellMax, NbCellLocal_Ncpus(commRank + 1)))
      allocate (FluxFourierFI &
                (NbNodeFaceMax, NbFracLocal_Ncpus(commRank + 1)))

#endif

   end subroutine Flux_allocate

   !> \brief Deallocate FluxDarcy and FluxFourier
   subroutine Flux_free

      deallocate (FluxDarcyKI)
      deallocate (FluxDarcyFI)

#ifdef _THERMIQUE_

      deallocate (FluxFourierKI)
      deallocate (FluxFourierFI)

#endif

   end subroutine Flux_free

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES

   subroutine Flux_compute_density_gravity_term(X1, rho1, X2, rho2, rho)
      type(TYPE_IncCVReservoir), intent(in) :: X1
      double precision, intent(in) :: rho1(NbPhase)
      type(TYPE_IncCVReservoir), intent(in) :: X2
      double precision, intent(in) :: rho2(NbPhase)
      double precision, intent(out) :: rho(NbPhase)

      integer :: k, phik, n(NbPhase)

      rho = 0.d0 ! should be ok by Fortran standard (intent(out))
      n = 0

      do k = 1, NbPhasePresente_ctx(X1%ic) ! set of present phases: Q_1
         phik = NumPhasePresente_ctx(k, X1%ic)
         rho(phik) = rho(phik) + rho1(phik)
         n(phik) = n(phik) + 1
      end do

      do k = 1, NbPhasePresente_ctx(X2%ic) ! set of present phases: Q_2
         phik = NumPhasePresente_ctx(k, X2%ic)
         rho(phik) = rho(phik) + rho2(phik)
         n(phik) = n(phik) + 1
      end do

      rho = rho/max(n, 1)

   end subroutine Flux_compute_density_gravity_term

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_EXTRAPOLATE

   subroutine Flux_compute_density_gravity_term(X1, rho1, X2, rho2, rho)
      type(TYPE_IncCVReservoir), intent(in) :: X1
      double precision, intent(in) :: rho1(NbPhase)
      type(TYPE_IncCVReservoir), intent(in) :: X2
      double precision, intent(in) :: rho2(NbPhase)
      double precision, intent(out) :: rho(NbPhase)

      integer :: k
      double precision :: rhoext, dP, dT, dC(NbComp)

      rho = 0.d0

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               rho(k) = 0.5*(rho1(k) + rho2(k))
            else
               ! FIXME: how to use phase pressure?
               call f_DensiteMassique(k, X2%Pression, X2%Temperature, X2%Comp, rhoext, dP, dT, dC)
               rho(k) = 0.5*(rho1(k) + rhoext)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) then
               ! FIXME: how to use phase pressure?
               call f_DensiteMassique(k, X1%Pression, X1%Temperature, X1%Comp, rhoext, dP, dT, dC)
               rho(k) = 0.5*(rhoext + rho2(k))
            endif
         endif
      end do

   end subroutine Flux_compute_density_gravity_term

! ComPASS_GRAVITY_AVERAGE_RHO_EXTRAPOLATE
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_MIN

   subroutine Flux_compute_density_gravity_term(X1, rho1, X2, rho2, rho)
      type(TYPE_IncCVReservoir), intent(in) :: X1
      double precision, intent(in) :: rho1(NbPhase)
      type(TYPE_IncCVReservoir), intent(in) :: X2
      double precision, intent(in) :: rho2(NbPhase)
      double precision, intent(out) :: rho(NbPhase)

      integer :: k

      rho = 0.d0 ! should be ok by Fortran standard (intent(out))

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               rho(k) = min(rho1(k), rho2(k))
            else
               rho(k) = rho1(k)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) rho(k) = rho2(k)
         endif
      end do

   end subroutine Flux_compute_density_gravity_term

! ComPASS_GRAVITY_AVERAGE_RHO_MIN
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY

   subroutine Flux_compute_density_gravity_term(X1, rho1, X2, rho2, rho)
      type(TYPE_IncCVReservoir), intent(in) :: X1
      double precision, intent(in) :: rho1(NbPhase)
      type(TYPE_IncCVReservoir), intent(in) :: X2
      double precision, intent(in) :: rho2(NbPhase)
      double precision, intent(out) :: rho(NbPhase)

      integer :: k
      double precision :: Stot

      rho = 0.d0 ! should be ok by Fortran standard (intent(out))

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               Stot = X1%Saturation(k) + X2%Saturation(k)
               if (Stot > epsilon_avrho) then
                  rho(k) = (X1%Saturation(k)*rho1(k) + X2%Saturation(k)*rho2(k))/Stot
               else
                  rho(k) = 0.5d0*(rho1(k) + rho2(k))
               end if
            else
               rho(k) = rho1(k)
            end if
         else
            if (phase_can_be_present(k, X2%ic)) rho(k) = rho2(k)
         endif
      end do

   end subroutine Flux_compute_density_gravity_term

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS_SINGULARITY
#endif

#ifdef ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS

   subroutine Flux_compute_density_gravity_term(X1, rho1, X2, rho2, rho)
      type(TYPE_IncCVReservoir), intent(in) :: X1
      double precision, intent(in) :: rho1(NbPhase)
      type(TYPE_IncCVReservoir), intent(in) :: X2
      double precision, intent(in) :: rho2(NbPhase)
      double precision, intent(out) :: rho(NbPhase)

      integer :: k
      double precision :: chi

      rho = 0.d0 ! should be ok by Fortran standard (intent(out))

      do k = 1, NbPhase
         if (phase_can_be_present(k, X1%ic)) then
            if (phase_can_be_present(k, X2%ic)) then
               if (X1%Saturation(k) < X2%Saturation(k)) then
                  if (X2%Saturation(k) > 0.d0) then
                     chi = X1%Saturation(k)/X2%Saturation(k)
                     rho(k) = (chi*rho1(k) + rho2(k))/(1 + chi)
                  endif
               else
                  if (X1%Saturation(k) > 0.d0) then
                     chi = X2%Saturation(k)/X1%Saturation(k)
                     rho(k) = (rho1(k) + chi*rho2(k))/(1 + chi)
                  endif
               endif
            else
               rho(k) = rho1(k)
            endif
         else
            if (phase_can_be_present(k, X2%ic)) rho(k) = rho2(k)
         endif
      end do

   end subroutine Flux_compute_density_gravity_term

! ComPASS_GRAVITY_AVERAGE_RHO_PRESENT_PHASES_CONTINUOUS
#endif

   !> \brief Structure of this subroutine:                             <br>
   !! loop of cell k                             <br>
   !!   a. loop of node i of cell k                             <br>
   !!       1. compute rho_ki_alpha                             <br>
   !!          loops of Q_k and Q_i                             <br>
   !!
   !!       2. loop of node j                             <br>
   !!          loops of Q_k and Q_i                             <br>
   !!
   !!       3. loop of frac j                             <br>
   !!          loops of Q_k and Q_i                             <br>
   !!
   !!   b. loop of frac i of cell k                             <br>
   !!       1. compute rho_ki_alpha                             <br>
   !!          loops of Q_k and Q_i                             <br>
   !!
   !!       2. loop of node j                             <br>
   !!          loops of Q_k and Q_i                             <br>
   !!
   !!       3. loop of frac j                             <br>
   !!          loops of Q_k and Q_i
   subroutine Flux_DarcyFlux_Cell() &
      bind(C, name="Flux_DarcyFlux_Cell")

      integer :: k, i, j, fj, fi
      integer :: numi, numj, alpha
      integer :: NbNodeCell, NbFracCell

      double precision :: rho_ki_alpha(NbPhase)
      double precision :: Tkij, dpkj, zkj

      FluxDarcyKI = 0.d0

      ! FluxDarcyKI
      do k = 1, NbCellLocal_Ncpus(commRank + 1) ! loop of cell

         ! number of nodes/fracs in cell k
         NbNodeCell = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         NbFracCell = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)

         ! i is node. Two parts:
         !   1. compute rho_ki^alpha
         !   2. loop of j
         do i = 1, NbNodeCell

            numi = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + i) ! num of node i

            call Flux_compute_density_gravity_term( &
               IncCell(k), DensiteMassiqueCell(:, k), &
               IncNode(numi), DensiteMassiqueNode(:, numi), &
               rho_ki_alpha)

            ! i is node, j is node
            do j = 1, NbNodeCell

               numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + j) ! num of node j

               Tkij = TkLocal_Darcy(k)%Pt(i, j) ! a_{k,s}^s'
               zkj = gravity*(XCellLocal(3, k) - XNodeLocal(3, numj)) ! g*(z_k - z_s')

               do alpha = 1, NbPhase
                  if (phase_can_be_present(alpha, IncCell(k)%ic) .or. &
                      phase_can_be_present(alpha, IncNode(numi)%ic)) then
                     dpkj = IncCell(k)%phase_pressure(alpha) - IncNode(numj)%phase_pressure(alpha)
                     FluxDarcyKI(alpha, i, k) = FluxDarcyKI(alpha, i, k) &
                                                + Tkij*(dpkj + rho_ki_alpha(alpha)*zkj)
                  end if
               end do

            end do

            ! i is node, j is frac
            do j = 1, NbFracCell

               numj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + j)
               fj = FracToFaceLocal(numj) ! fj is face number, numj is frac number

               Tkij = TkLocal_Darcy(k)%Pt(i, j + NbNodeCell) ! a_{k,s}^s'
               zkj = gravity*(XCellLocal(3, k) - XFaceLocal(3, fj)) ! g*(z_k - z_s')

               do alpha = 1, NbPhase
                  if (phase_can_be_present(alpha, IncCell(k)%ic) .or. &
                      phase_can_be_present(alpha, IncNode(numi)%ic)) then
                     dpkj = IncCell(k)%phase_pressure(alpha) - IncFrac(numj)%phase_pressure(alpha)
                     FluxDarcyKI(alpha, i, k) = FluxDarcyKI(alpha, i, k) &
                                                + Tkij*(dpkj + rho_ki_alpha(alpha)*zkj)
                  end if
               end do

            end do ! end of j
         end do

         ! i is frac. Two parts:
         !   1. compute rho_ki^alpha
         !   2. loop of j
         do i = 1, NbFracCell

            numi = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + i) ! numi is frac number
            fi = FracToFaceLocal(numi) ! fi is face number

            call Flux_compute_density_gravity_term( &
               IncCell(k), DensiteMassiqueCell(:, k), &
               IncFrac(numi), DensiteMassiqueFrac(:, numi), &
               rho_ki_alpha)

            ! i is frac, j is node
            do j = 1, NbNodeCell

               numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + j)

               Tkij = TkLocal_Darcy(k)%Pt(i + NbNodeCell, j)            ! a_{k,s}^s'
               zkj = gravity*(XCellLocal(3, k) - XNodeLocal(3, numj)) ! g*(z_k - z_s')

               do alpha = 1, NbPhase
                  if (phase_can_be_present(alpha, IncCell(k)%ic) .or. &
                      phase_can_be_present(alpha, IncFrac(numi)%ic)) then
                     dpkj = IncCell(k)%phase_pressure(alpha) - IncNode(numj)%phase_pressure(alpha)
                     FluxDarcyKI(alpha, i + NbNodeCell, k) = FluxDarcyKI(alpha, i + NbNodeCell, k) &
                                                             + Tkij*(dpkj + rho_ki_alpha(alpha)*zkj)
                  end if
               end do

            end do ! end of j

            ! i is frac, j is frac
            do j = 1, NbFracCell

               numj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + j)
               fj = FracToFaceLocal(numj) ! fj is face number, numj is frac number

               Tkij = TkLocal_Darcy(k)%Pt(i + NbNodeCell, j + NbNodeCell) ! a_{k,s}^s'
               zkj = gravity*(XCellLocal(3, k) - XFaceLocal(3, fj))   ! g*(z_k - z_s')

               do alpha = 1, NbPhase
                  if (phase_can_be_present(alpha, IncCell(k)%ic) .or. &
                      phase_can_be_present(alpha, IncFrac(numi)%ic)) then
                     dpkj = IncCell(k)%phase_pressure(alpha) - IncFrac(numj)%phase_pressure(alpha)
                     FluxDarcyKI(alpha, i + NbNodeCell, k) = FluxDarcyKI(alpha, i + NbNodeCell, k) &
                                                             + Tkij*(dpkj + rho_ki_alpha(alpha)*zkj)
                  end if
               end do

            end do ! fin loop of j

         end do ! fin loop of i in cell k

      end do ! fin loop cell k

   end subroutine Flux_DarcyFlux_Cell

   ! Structure of this subroutine:
   ! loop of frac k
   !   a. loop of node i of frac k
   !       1. compute rho_ki_alpha
   !          loops of Q_k and Q_i

   !       2. loop of node j
   !          loops of Q_k and Q_i
   subroutine Flux_DarcyFlux_Frac() &
      bind(C, name="Flux_DarcyFlux_Frac")

      integer :: k, fk, i, j, numi, numj, tmp_compt(NbPhase)
      integer :: alpha, nph_i, numph_i, nph_k, numph_k
      integer :: NbNodeFrac

      double precision :: dpkj, Tkij, zkj

      double precision :: rho_ki_alpha(NbPhase)

      FluxDarcyFI(:, :, :) = 0

      ! loop of frac
      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         fk = FracToFaceLocal(k) ! fk is face number

         ! number of nodes in a frac
         NbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)

         do i = 1, NbNodeFrac

            numi = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + i) ! node number of i

            ! compute rho_ki^alpha: loop of Q_k and loop of Q_i
            rho_ki_alpha(:) = 0.d0
            tmp_compt(:) = 0

            do nph_k = 1, NbPhasePresente_ctx(IncFrac(k)%ic) ! phases present: Q_k
               numph_k = NumPhasePresente_ctx(nph_k, IncFrac(k)%ic)
               rho_ki_alpha(numph_k) = rho_ki_alpha(numph_k) + DensiteMassiqueFrac(numph_k, k)
               tmp_compt(numph_k) = tmp_compt(numph_k) + 1
            end do

            do nph_i = 1, NbPhasePresente_ctx(IncNode(numi)%ic) ! phases present: Q_i
               numph_i = NumPhasePresente_ctx(nph_i, IncNode(numi)%ic)
               rho_ki_alpha(numph_i) = rho_ki_alpha(numph_i) + DensiteMassiqueNode(numph_i, numi) ! Attention: numph_k used for densitemassique
               tmp_compt(numph_i) = tmp_compt(numph_i) + 1
            end do

            rho_ki_alpha(:) = rho_ki_alpha(:)/max(tmp_compt(:), 1)

            do j = 1, NbNodeFrac
               numj = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + j)

               Tkij = TkFracLocal_Darcy(k)%Pt(i, j) ! a_{k,s}^s'
               zkj = gravity*(XFaceLocal(3, fk) - XNodeLocal(3, numj)) ! g*(z_k - z_s')

               do alpha = 1, NbPhase
                  if (phase_can_be_present(alpha, IncFrac(k)%ic) .or. &
                      phase_can_be_present(alpha, IncNode(numi)%ic)) then
                     dpkj = IncFrac(k)%phase_pressure(alpha) - IncNode(numj)%phase_pressure(alpha)
                     FluxDarcyFI(alpha, i, k) = FluxDarcyFI(alpha, i, k) &
                                                + Tkij*(dpkj + rho_ki_alpha(alpha)*zkj)
                  end if
               end do

            end do
         end do

      end do ! end loop of frac k

   end subroutine Flux_DarcyFlux_Frac

   subroutine Flux_FourierFlux_Cell() &
      bind(C, name="Flux_FourierFlux_Cell")

      integer :: k, i, j
      integer :: numj
      integer :: NbNodeCell, NbFracCell

      FluxFourierKI(:, :) = 0.d0

      ! FluxFourierKI
      do k = 1, NbCellLocal_Ncpus(commRank + 1) ! loop of cell

         ! number of nodes/fracs in cell k
         NbNodeCell = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         NbFracCell = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)

         ! i is node
         do i = 1, NbNodeCell

            ! i is node, j is node
            do j = 1, NbNodeCell

               numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + j)

               FluxFourierKI(i, k) = FluxFourierKI(i, k) &
                                     + TkLocal_Fourier(k)%Pt(i, j)* &
                                     (IncCell(k)%Temperature - IncNode(numj)%Temperature)
            end do

            ! i is node, j is frac
            do j = 1, NbFracCell

               numj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + j)

               FluxFourierKI(i, k) = FluxFourierKI(i, k) &
                                     + TkLocal_Fourier(k)%Pt(i, j + NbNodeCell) &
                                     *(IncCell(k)%Temperature - IncFrac(numj)%Temperature)
            end do
         end do

         ! i is frac
         do i = 1, NbFracCell

            ! i is frac, j is node
            do j = 1, NbNodeCell

               numj = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + j)

               FluxFourierKI(i + NbNodeCell, k) = FluxFourierKI(i + NbNodeCell, k) &
                                                  + TkLocal_Fourier(k)%Pt(i + NbNodeCell, j) &
                                                  *(IncCell(k)%Temperature - IncNode(numj)%Temperature)
            end do

            ! i is frac, j is frac
            do j = 1, NbFracCell

               numj = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + j)

               FluxFourierKI(i + NbNodeCell, k) = FluxFourierKI(i + NbNodeCell, k) &
                                                  + TkLocal_Fourier(k)%Pt(i + NbNodeCell, j + NbNodeCell) &
                                                  *(IncCell(k)%Temperature - IncFrac(numj)%Temperature)
            end do
         end do

      end do ! fin loop cell: k

   end subroutine Flux_FourierFlux_Cell

   subroutine Flux_FourierFlux_Frac() &
      bind(C, name="Flux_FourierFlux_Frac")

      integer :: k, fk, i, j, numj
      integer :: NbNodeFrac

      FluxFourierFI(:, :) = 0

      ! loop of frac
      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         fk = FracToFaceLocal(k) ! fk is face number

         ! number of nodes in a frac
         NbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)

         do i = 1, NbNodeFrac

            do j = 1, NbNodeFrac
               numj = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + j)

               FluxFourierFI(i, k) = FluxFourierFI(i, k) &
                                     + TkFracLocal_Fourier(k)%pt(i, j) &
                                     *(IncFrac(k)%Temperature - IncNode(numj)%Temperature)
            end do
         end do

      end do ! end loop of frac k

   end subroutine Flux_FourierFlux_Frac

end module Flux
