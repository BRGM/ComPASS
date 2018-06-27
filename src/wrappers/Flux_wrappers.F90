!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module FluxWrappers

   use CommonTypesWrapper
   use DefModel
   use IncCVReservoir
   use LoisThermoHydro
   use Flux

   implicit none

   private :: &
      ks_mass_fluxes, &
      average_cell_mass_fluxes, &
      average_fracture_mass_fluxes

   public :: &
      retrieve_mass_fluxes

contains

   subroutine ks_mass_fluxes(k, Xk, Inck, Mk, sr, Xs, Incs, Ms, Fks)

      integer, intent(in) :: k, sr ! sr is s index relative to k
      real(c_double), intent(in) :: Xk(3), Xs(3)
      type(TYPE_IncCVReservoir), intent(in) :: Inck, Incs
      real(c_double), intent(in) :: Mk(:, :), Ms(:, :)
      real(c_double), intent(out) :: Fks(3, NbComp)

      integer :: m, mph, icp

      real(c_double) :: flux_ks(NbComp), DarcyFlux, XkXs(3)

      flux_ks(:) = 0.d0
      do m = 1, NbPhasePresente_ctx(Inck%ic) ! Q_k
         mph = NumPhasePresente_ctx(m, Inck%ic)
         DarcyFlux = FluxDarcyKI(mph, sr, k)
         if (DarcyFlux >= 0.d0) then ! K_{k,s}^{alpha}=k
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i
                  flux_ks(icp) = flux_ks(icp) + Mk(icp, m)*DarcyFlux
               end if
            end do
         end if
      end do
      do m = 1, NbPhasePresente_ctx(Incs%ic) ! Q_s
         mph = NumPhasePresente_ctx(m, Incs%ic)
         DarcyFlux = FluxDarcyKI(mph, sr, k)
         if (DarcyFlux < 0.d0) then ! K_{k,s}^{alpha}=s
            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! \cap P_i
                  flux_ks(icp) = flux_ks(icp) + Ms(icp, m)*DarcyFlux
               end if
            end do
         end if
      end do
      XkXs = Xs - Xk
      do icp = 1, NbComp
         Fks(:, icp) = Fks(:, icp) + flux_ks(icp)*XkXs
      end do

   end subroutine ks_mass_fluxes

   subroutine average_cell_mass_fluxes(mass_fluxes)

      real(c_double), intent(out) :: mass_fluxes(:, :, :)

      integer :: k, s, fs, nums, NbNodeCell, NbFracCell
      real(c_double) :: Xk(3)

      integer :: errcode, Ierr

      if (.not. (size(mass_fluxes, 1) == 3 .and. &
                 size(mass_fluxes, 2) == NbComp .and. &
                 size(mass_fluxes, 3) == NbCellLocal_Ncpus(commRank + 1))) then
         write (*, *) 'Output cell mass fluxes should have size:', &
            (/3, NbComp, NbCellLocal_Ncpus(commRank + 1)/)
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      mass_fluxes(:, :, :) = 0.d0

      do k = 1, NbCellLocal_Ncpus(commRank + 1)
         Xk(:) = XCellLocal(:, k)
         NbNodeCell = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         NbFracCell = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)
         do s = 1, NbNodeCell
            nums = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + s)
            call ks_mass_fluxes( &
               k, Xk, IncCell(k), DensitemolaireKrViscoCompCell(:, :, k), &
               s, XNodeLocal(:, nums), IncNode(nums), DensitemolaireKrViscoCompNode(:, :, nums), &
               mass_fluxes(:, :, k))
         end do
         do s = 1, NbFracCell
            fs = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + s) ! fs is face number
            nums = FaceToFracLocal(fs) ! nums is frac number
            call ks_mass_fluxes( &
               k, Xk, IncCell(k), DensitemolaireKrViscoCompCell(:, :, k), &
               s + NbNodeCell, XFaceLocal(:, fs), IncFrac(nums), DensitemolaireKrViscoCompFrac(:, :, nums), &
               mass_fluxes(:, :, k))
         end do
      end do

   end subroutine average_cell_mass_fluxes

   subroutine average_fracture_mass_fluxes(mass_fluxes)

      real(c_double), intent(out) :: mass_fluxes(:, :, :)

      integer :: k, s, fk, nums, NbNodeFrac

      integer :: errcode, Ierr

      if (.not. (size(mass_fluxes, 1) == 3 .and. &
                 size(mass_fluxes, 2) == NbComp .and. &
                 size(mass_fluxes, 3) == NbFracLocal_Ncpus(commRank + 1))) then
         write (*, *) 'Output fracture mass fluxes should have size:', &
            (/3, NbComp, NbFracLocal_Ncpus(commRank + 1)/)
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      mass_fluxes(:, :, :) = 0.d0

      do k = 1, NbFracLocal_Ncpus(commRank + 1)
         fk = FracToFaceLocal(k) ! fk is face number
         NbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)
         do s = 1, NbNodeFrac
            nums = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + s)
            call ks_mass_fluxes( &
               k, XFaceLocal(:, fk), IncFrac(k), DensitemolaireKrViscoCompFrac(:, :, k), &
               s, XNodeLocal(:, nums), IncNode(nums), DensitemolaireKrViscoCompNode(:, :, nums), &
               mass_fluxes(:, :, k))
         end do
      end do

   end subroutine average_fracture_mass_fluxes

   subroutine retrieve_mass_fluxes(Fc, Ff) &
      bind(C, name="retrieve_mass_fluxes")

      type(cpp_narray_wrapper), intent(in) :: Fc, Ff

      real(c_double), pointer :: Fcp(:, :, :), Ffp(:, :, :)

      call bind_3array(Fc, Fcp)
      call average_cell_mass_fluxes(Fcp)
      call bind_3array(Ff, Ffp)
      call average_fracture_mass_fluxes(Ffp)

   end subroutine retrieve_mass_fluxes

end module FluxWrappers
