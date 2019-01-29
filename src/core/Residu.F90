!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Residu

   use, intrinsic :: iso_c_binding 

   use DefModel
   use IncCVReservoir
   use IncCVWells
   use NeumannContribution
   use LoisThermoHydro
   use Flux
   use InteroperabilityStructures

   implicit none

   ! Residu cell, fracture face, node
   real(c_double), pointer :: &
      ResiduCell(:, :), &
      ResiduFrac(:, :), &
      ResiduNode(:, :)

   ! Residu injection well and production well
   real(c_double), pointer :: &
      ResiduWellInj(:), &
      ResiduWellProd(:)

   ! AccVol of time step n-1
   real(c_double), pointer :: &
      AccVolCell_1(:, :), &
      AccVolFrac_1(:, :), &
      AccVolNode_1(:, :)

   type, bind(C) :: CTVector
      real(c_double) :: values(NbCompThermique)
   end type CTVector

   public :: &
      Residu_reset_history, &
      Residu_compute, &
      Residu_RelativeNorm, &
      Residu_associate_pointers

   private :: &
      Residu_AccVol, &
      Residu_add_flux_contributions, &
      Residu_reset_Dirichlet_nodes, &
      Residu_add_Neumann_contributions, &
      Residu_add_flux_contributions_wells, &
      Residu_clear_absent_components_accumulation

contains

   subroutine Residu_associate_pointers( &
      nodes, fractures, cells, injectors, producers, &
      nodes_accumulation, fractures_accumulation, cells_accumulation &
   ) &
   bind(C, name="Residu_associate_pointers")
   
      type(cpp_array_wrapper), intent(in), value :: &
         nodes, cells, fractures, injectors, producers, &
         nodes_accumulation, cells_accumulation, fractures_accumulation
   
      if(nodes%n/=NbNodeLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent nodes residual size')
      if(nodes_accumulation%n/=NbNodeLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent accumulation nodes size')
      if(cells%n/=NbCellLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent cells residual size')
      if(cells_accumulation%n/=NbCellLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent residual accumulation cells size')
      if(fractures%n/=NbFracLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent residual fractures size')
      if(fractures_accumulation%n/=NbFracLocal_Ncpus(commRank + 1) * NbCompThermique) &
         call CommonMPI_abort('inconsistent residual accumulation fractures size')
      if(injectors%n/=NbWellInjLocal_Ncpus(commRank + 1)) &
         call CommonMPI_abort('inconsistent injectors residual size')
      if(producers%n/=NbWellProdLocal_Ncpus(commRank + 1)) &
         call CommonMPI_abort('inconsistent producers residual size')
   
      call c_f_pointer(nodes%p, ResiduNode, [NbCompThermique, NbNodeLocal_Ncpus(commRank + 1)])
      call c_f_pointer(cells%p, ResiduCell, [NbCompThermique, NbCellLocal_Ncpus(commRank + 1)])
      call c_f_pointer(fractures%p, ResiduFrac, [NbCompThermique, NbFracLocal_Ncpus(commRank + 1)])
      call c_f_pointer(injectors%p, ResiduWellInj, [NbWellInjLocal_Ncpus(commRank + 1)])
      call c_f_pointer(producers%p, ResiduWellProd, [NbWellProdLocal_Ncpus(commRank + 1)])
      call c_f_pointer(nodes_accumulation%p, AccVolNode_1, [NbCompThermique, NbNodeLocal_Ncpus(commRank + 1)])
      call c_f_pointer(cells_accumulation%p, AccVolCell_1, [NbCompThermique, NbCellLocal_Ncpus(commRank + 1)])
      call c_f_pointer(fractures_accumulation%p, AccVolFrac_1, [NbCompThermique, NbFracLocal_Ncpus(commRank + 1)])
   
   end subroutine Residu_associate_pointers

   ! Newton is initialized with the unknown at time step n-1
   subroutine Residu_reset_history() &
      bind(C, name="Residu_reset_history")

      integer :: i, k

      call Residu_AccVol

      do k = 1, NbCellLocal_Ncpus(commRank + 1) ! cell
         do i = 1, NbCompThermique
            AccVolCell_1(i, k) = IncCell(k)%AccVol(i)
         end do
      end do

      do k = 1, NbFracLocal_Ncpus(commRank + 1) ! frac
         do i = 1, NbCompThermique
            AccVolFrac_1(i, k) = IncFrac(k)%AccVol(i)
         end do
      end do

      do k = 1, NbNodeLocal_Ncpus(commRank + 1) ! node
         do i = 1, NbCompThermique
            AccVolNode_1(i, k) = IncNode(k)%AccVol(i)
         end do
      end do

   end subroutine Residu_reset_history

   subroutine Residu_compute(Delta_t) &
      bind(C, name="Residu_compute")

      real(c_double), intent(in), value :: Delta_t

      integer :: i, k

      ! init Residu as zero
      ResiduCell(:, :) = 0.d0
      ResiduFrac(:, :) = 0.d0
      ResiduNode(:, :) = 0.d0

      ResiduWellInj(:) = 0.d0
      ResiduWellProd(:) = 0.d0

      call Residu_AccVol

      ! Residu from terms of accumulation
      do k = 1, NbCellLocal_Ncpus(commRank + 1) ! cell
         do i = 1, NbCompThermique
            ResiduCell(i, k) = &
               (IncCell(k)%AccVol(i) - AccVolCell_1(i, k))/Delta_t
         end do
      end do

      do k = 1, NbFracLocal_Ncpus(commRank + 1) ! Frac
         do i = 1, NbCompThermique
            ResiduFrac(i, k) = &
               (IncFrac(k)%AccVol(i) - AccVolFrac_1(i, k))/Delta_t
         end do
      end do

      do k = 1, NbNodeLocal_Ncpus(commRank + 1) ! node
         do i = 1, NbCompThermique
            ResiduNode(i, k) = &
               (IncNode(k)%AccVol(i) - AccVolNode_1(i, k))/Delta_t
         end do
      end do

#ifdef _THERMIQUE_
      ! Thermal source
      ResiduCell(NbComp + 1, :) = ResiduCell(NbComp + 1, :) - CellThermalSourceVol
      ResiduFrac(NbComp + 1, :) = ResiduFrac(NbComp + 1, :) - FracThermalSourceVol
      ResiduNode(NbComp + 1, :) = ResiduNode(NbComp + 1, :) - NodeThermalSourceVol
#endif

      ! Residu for dir node will be reset as 0 in the end of this subroutine
      ! We do not consider dirichlet nodes during computing the residu such that
      ! the code is clearer

      ! Residu from conservation composants
      call Residu_add_flux_contributions

      call Residu_reset_Dirichlet_nodes

      call Residu_add_Neumann_contributions

   end subroutine Residu_compute

   subroutine Residu_reset_Dirichlet_nodes

      integer :: k

      ! Residu for dir node should be reset as 0 for computing norm of Residu
#ifdef _THERMIQUE_
      do k = 1, NbNodeLocal_Ncpus(commRank + 1) ! node

         if (IdNodeLocal(k)%P == "d") then
            ResiduNode(1, k) = 0.d0
            ResiduNode(3:NbCompThermique, k) = 0.d0
         end if

         if (IdNodeLocal(k)%T == "d") then
            ResiduNode(2, k) = 0.d0
         end if
      end do
#else
      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         if (IdNodeLocal(k)%P == "d") then
            ResiduNode(1:NbCompThermique, k) = 0.d0
         end if
      end do
#endif

   end subroutine Residu_reset_Dirichlet_nodes

   subroutine Residu_add_Neumann_contributions

      integer :: k

      do k = 1, NbNodeLocal_Ncpus(commRank + 1) ! node
         ResiduNode(1:NbComp, k) = ResiduNode(1:NbComp, k) - NodeNeumannBC(k)%molar_flux
#ifdef _THERMIQUE_
         ResiduNode(NbCompThermique, k) = ResiduNode(NbCompThermique, k) - NodeNeumannBC(k)%heat_flux
#endif
      end do

   end subroutine Residu_add_Neumann_contributions

   ! FIXME: Could be simpler if we multiply accumulations by a mask with 0 and 1
   subroutine Residu_clear_absent_components_accumulation(ic, accumulations)

      integer, intent(in) :: ic ! context
      real(c_double), intent(inout) :: accumulations(NbCompThermique) ! unknowns
      integer :: i, icp
      real(c_double) :: copy(NbComp) ! unknowns

      do i = 1, NbCompCtilde_ctx(ic)
         icp = NumCompCtilde_ctx(i, ic)
         copy(i) = accumulations(icp)
      end do

      ! CHECKME: is size differnent from NbCompThermique: inc(:) = 0.d0 would be enough
      accumulations(1:NbCompThermique) = 0.d0

      do i = 1, NbCompCtilde_ctx(ic)
         icp = NumCompCtilde_ctx(i, ic)
         accumulations(icp) = copy(i)
      end do

   end subroutine Residu_clear_absent_components_accumulation

   subroutine Residu_AccVol

      integer :: k, m, mph, i, ic, icp
      real(c_double) :: tmpAccVol(NbComp)

      ! Cells
      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         ic = IncCell(k)%ic

         call Residu_clear_absent_components_accumulation(ic, IncCell(k)%AccVol)

         do m = 1, NbPhasePresente_ctx(ic) ! Q_k
            mph = NumPhasePresente_ctx(m, ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  IncCell(k)%AccVol(icp) = IncCell(k)%AccVol(icp) &
                                           + PoroVolDarcyCell(k)*DensitemolaireSatCompCell(icp, m, k)
               end if
            end do

#ifdef _THERMIQUE_

            IncCell(k)%AccVol(NbComp + 1) = IncCell(k)%AccVol(NbComp + 1) &
                                            + PoroVolFourierCell(k)*DensitemolaireEnergieInterneSatCell(m, k)
#endif
         end do ! end of phase m

#ifdef _THERMIQUE_

         IncCell(k)%AccVol(NbComp + 1) = IncCell(k)%AccVol(NbComp + 1) &
                                         + Poro_1volFourierCell(k)*CpRoche*IncCell(k)%Temperature
#endif

      end do ! end of cell k

      ! Fractures
      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         ic = IncFrac(k)%ic

         call Residu_clear_absent_components_accumulation(ic, IncFrac(k)%AccVol)

         do m = 1, NbPhasePresente_ctx(ic) ! Q_k
            mph = NumPhasePresente_ctx(m, ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  IncFrac(k)%AccVol(icp) = IncFrac(k)%AccVol(icp) &
                                           + PoroVolDarcyFrac(k)*DensitemolaireSatCompFrac(icp, m, k)
               end if
            end do

#ifdef _THERMIQUE_

            IncFrac(k)%AccVol(NbComp + 1) = IncFrac(k)%AccVol(NbComp + 1) &
                                            + PoroVolFourierFrac(k)*DensitemolaireEnergieInterneSatFrac(m, k)
#endif

         end do ! end of phase m

#ifdef _THERMIQUE_

         IncFrac(k)%AccVol(NbComp + 1) = IncFrac(k)%AccVol(NbComp + 1) &
                                         + Poro_1volFourierFrac(k)*CpRoche*IncFrac(k)%Temperature
#endif

      end do ! end of frac k

      ! Nodes
      do k = 1, NbNodeLocal_Ncpus(commRank + 1)

         ic = IncNode(k)%ic

         call Residu_clear_absent_components_accumulation(ic, IncNode(k)%AccVol)

         do m = 1, NbPhasePresente_ctx(IncNode(k)%ic) ! Q_k
            mph = NumPhasePresente_ctx(m, IncNode(k)%ic)

            do icp = 1, NbComp
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                  IncNode(k)%AccVol(icp) = IncNode(k)%AccVol(icp) &
                                           + PoroVolDarcyNode(k)*DensitemolaireSatCompNode(icp, m, k)
               end if
            end do

#ifdef _THERMIQUE_

            IncNode(k)%AccVol(NbComp + 1) = IncNode(k)%AccVol(NbComp + 1) &
                                            + PoroVolFourierNode(k)*DensitemolaireEnergieInterneSatNode(m, k)

#endif

         end do ! end of phase m

#ifdef _THERMIQUE_

         IncNode(k)%AccVol(NbComp + 1) = IncNode(k)%AccVol(NbComp + 1) &
                                         + Poro_1volFourierNode(k)*CpRoche*IncNode(k)%Temperature

#endif

      end do ! end of node k

   end subroutine Residu_AccVol

   subroutine Residu_add_flux_contributions

      integer :: k, s, fs, fk, nums, m, mph, icp
      integer :: NbNodeCell, NbFracCell, NbNodeFrac

      double precision :: Flux_ks(NbComp), FluxT_ks, DarcyFlux
      double precision :: Pws, Tws, Ps, Ts, WIDws, WIFws, qw, Ps_Pws

      do k = 1, NbCellLocal_Ncpus(commRank + 1)

         ! number of nodes/fracs in cell k
         NbNodeCell = NodebyCellLocal%Pt(k + 1) - NodebyCellLocal%Pt(k)
         NbFracCell = FracbyCellLocal%Pt(k + 1) - FracbyCellLocal%Pt(k)

         ! s is node in cell k
         do s = 1, NbNodeCell
            nums = NodebyCellLocal%Num(NodebyCellLocal%Pt(k) + s)

            Flux_ks(:) = 0.d0
            FluxT_ks = 0.d0

            ! alpha in Q_{K_{k,s}^{alpha}} =
            !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
            ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
            ! OPTIMIZE: could directly be an array with NumPhasePresente?
            do m = 1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
               mph = NumPhasePresente_ctx(m, IncCell(k)%ic)

               DarcyFlux = FluxDarcyKI(mph, s, k)

               if (DarcyFlux >= 0.d0) then ! K_{k,s}^{alpha}=k

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompCell(icp, m, k)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieCell(m, k) &
                             *DarcyFlux

#endif

               end if
            end do

            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

               DarcyFlux = FluxDarcyKI(mph, s, k)

               if (DarcyFlux < 0.d0) then ! K_{k,s}^{alpha}=s

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompNode(icp, m, nums)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m, nums) &
                             *DarcyFlux

#endif

               end if
            end do

            ResiduCell(1:NbComp, k) = ResiduCell(1:NbComp, k) + Flux_ks(:)
            ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) - Flux_ks(:)

#ifdef _THERMIQUE_

            ResiduCell(NbComp + 1, k) = ResiduCell(NbComp + 1, k) &
                                        + FluxT_ks &
                                        + FluxFourierKI(s, k)

            ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) &
                                           - FluxT_ks &
                                           - FluxFourierKI(s, k)
#endif

         end do ! end of node s in cell k

         ! s is frac in cell k
         do s = 1, NbFracCell

            fs = FracbyCellLocal%Num(FracbyCellLocal%Pt(k) + s)
            nums = FaceToFracLocal(fs) ! fs is face number, nums is frac number

            Flux_ks(:) = 0.d0
            FluxT_ks = 0.d0

            ! alpha in Q_{K_{k,s}^{alpha}} =
            !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
            ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
            do m = 1, NbPhasePresente_ctx(IncCell(k)%ic) ! Q_k
               mph = NumPhasePresente_ctx(m, IncCell(k)%ic)

               DarcyFlux = FluxDarcyKI(mph, s + NbNodeCell, k)

               if (DarcyFlux >= 0.d0) then ! K_{k,s}^{alpha}=k

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompCell(icp, m, k)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieCell(m, k) &
                             *DarcyFlux
#endif

               end if
            end do

            do m = 1, NbPhasePresente_ctx(IncFrac(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncFrac(nums)%ic)

               DarcyFlux = FluxDarcyKI(mph, s + NbNodeCell, k)

               if (DarcyFlux < 0.d0) then ! K_{k,s}^{alpha}=s

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompFrac(icp, m, nums)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieFrac(m, nums) &
                             *DarcyFlux
#endif

               end if
            end do

            ResiduCell(1:NbComp, k) = ResiduCell(1:NbComp, k) + Flux_ks(:)
            ResiduFrac(1:NbComp, nums) = ResiduFrac(1:NbComp, nums) - Flux_ks(:)

#ifdef _THERMIQUE_

            ResiduCell(NbComp + 1, k) = ResiduCell(NbComp + 1, k) &
                                        + FluxT_ks &
                                        + FluxFourierKI(s + NbNodeCell, k)

            ResiduFrac(NbComp + 1, nums) = ResiduFrac(NbComp + 1, nums) &
                                           - FluxT_ks &
                                           - FluxFourierKI(s + NbNodeCell, k)
#endif

         end do ! end of frac s in cell k

      end do ! end of cell k

      do k = 1, NbFracLocal_Ncpus(commRank + 1)

         fk = FracToFaceLocal(k) ! fk is face number

         ! number of nodes in a frac
         NbNodeFrac = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)

         do s = 1, NbNodeFrac
            nums = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(fk) + s)

            Flux_ks(:) = 0.d0
            FluxT_ks = 0.d0

            ! alpha in Q_{K_{k,s}^{alpha}} =
            !      ( {alpha | alpha in Q_k} \cap {alpha | K_{k,s}^{alpha}=k})
            ! \cup ( {alpha | alpha in Q_s} \cap {alpha | K_{k,s}^{alpha}=s})
            do m = 1, NbPhasePresente_ctx(IncFrac(k)%ic) ! Q_k
               mph = NumPhasePresente_ctx(m, IncFrac(k)%ic)

               if (FluxDarcyFI(mph, s, k) >= 0.d0) then ! K_{k,s}^{alpha}=k

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompFrac(icp, m, k)*FluxDarcyFI(mph, s, k)

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieFrac(m, k) &
                             *FluxDarcyFI(mph, s, k)
#endif
               end if
            end do

            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

               if (FluxDarcyFI(mph, s, k) < 0.d0) then ! K_{k,s}^{alpha}=s

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i

                        Flux_ks(icp) = Flux_ks(icp) &
                                       + DensitemolaireKrViscoCompNode(icp, m, nums)*FluxDarcyFI(mph, s, k)

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m, nums) &
                             *FluxDarcyFI(mph, s, k)
#endif
               end if
            end do

            ResiduFrac(1:NbComp, k) = ResiduFrac(1:NbComp, k) + Flux_ks(:)
            ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) - Flux_ks(:)

#ifdef _THERMIQUE_

            ResiduFrac(NbComp + 1, k) = ResiduFrac(NbComp + 1, k) &
                                        + FluxT_ks &
                                        + FluxFourierFI(s, k)

            ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) &
                                           - FluxT_ks &
                                           - FluxFourierFI(s, k)
#endif

         end do ! end of node s in frac k

      end do ! end of frac k

      call Residu_add_flux_contributions_wells

   end subroutine Residu_add_flux_contributions

   subroutine Residu_add_flux_contributions_wells

      integer :: k, s, fs, fk, nums, m, mph, icp

      double precision :: Flux_ks(NbComp), FluxT_ks
      double precision :: Pws, Tws, Ps, Ts, WIDws, WIFws, qw, Ps_Pws

      ! Injection well
      ! It is possible that one ghost well contains some own nodes,
      ! thus, when assembly of flux, we need to take into account
      ! not only own well, but also ghost well
      do k = 1, NbWellInjLocal_Ncpus(commRank + 1) ! injection well k

         qw = 0.d0

         ! nodes of well k
         do s = NodebyWellInjLocal%Pt(k) + 1, NodebyWellInjLocal%Pt(k + 1)
            nums = NodebyWellInjLocal%Num(s)

            Pws = PerfoWellInj(s)%Pression ! P_{w,s}
            Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
            Ps = IncNode(nums)%Pression ! P_s
            Ts = IncNode(nums)%Temperature ! T_s

            WIDws = NodeDatabyWellInjLocal%Val(s)%WID
            WIFws = NodeDatabyWellInjLocal%Val(s)%WIF

            Flux_ks(:) = 0.d0
            FluxT_ks = 0.d0

            ! Flux q_{w,s,i} and q_{w,s,e}
            Ps_Pws = Ps - Pws

            if (Ps_Pws < 0.d0) then ! < 0

               do icp = 1, NbComp
                  Flux_ks(icp) = DensitemolaireKrViscoCompWellInj(icp, s)*WIDws*Ps_Pws
               end do
#ifdef _THERMIQUE_
               FluxT_ks = DensitemolaireKrViscoEnthalpieWellInj(s)*WIDws*Ps_Pws
#endif
            end if

            ! node equation
            ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) + Flux_ks(:)

#ifdef _THERMIQUE_
            ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) + FluxT_ks ! + WIFws * (Ts-Tws)
#endif
            ! qw
            do icp = 1, NbComp
               qw = qw + Flux_ks(icp)
            end do
         end do

         ! inj well equation
         if (DataWellInjLocal(k)%IndWell == 'p') then
            ResiduWellInj(k) = DataWellInjLocal(k)%PressionMax - IncPressionWellInj(k)
         else if (DataWellInjLocal(k)%IndWell == 'f') then
            ResiduWellInj(k) = qw - DataWellInjLocal(k)%ImposedFlowrate
         else
            call CommonMPI_abort("Injection well index should be 'p' (pressure mode) or 'f' (flowrate mode)!")
         end if

      end do ! end of node s in injection well k

      ! Production well
      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)

         qw = 0.d0

         ! nodes of well k
         do s = NodebyWellProdLocal%Pt(k) + 1, NodebyWellProdLocal%Pt(k + 1)
            nums = NodebyWellProdLocal%Num(s)

            Pws = PerfoWellProd(s)%Pression ! P_{w,s}
            Ps = IncNode(nums)%Pression ! P_s
            Ps_Pws = Ps - Pws
            WIDws = NodeDatabyWellProdLocal%Val(s)%WID

            ! print*, "ps_pw: ", Ps_Pws

            Flux_ks(:) = 0.d0
            FluxT_ks = 0.d0

            do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
               mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

               if (Ps_Pws > 0.d0) then

                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i
                        Flux_ks(icp) = Flux_ks(icp) + DensitemolaireKrViscoCompNode(icp, m, nums)*WIDws*Ps_Pws
                     end if
                  end do

#ifdef _THERMIQUE_
                  FluxT_ks = FluxT_ks + DensitemolaireKrViscoEnthalpieNode(m, nums)*WIDws*Ps_Pws
#endif
               end if
            end do

            ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) + Flux_ks(:)

#ifdef _THERMIQUE_
            ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) + FluxT_ks
#endif
            ! qw
            do icp = 1, NbComp
               qw = qw + Flux_ks(icp)
            end do
         end do

         ! prod well equation
         if (DataWellProdLocal(k)%IndWell == 'p') then
            ResiduWellProd(k) = IncPressionWellProd(k) - DataWellProdLocal(k)%PressionMin
         else if (DataWellProdLocal(k)%IndWell == 'f') then
            ResiduWellProd(k) = DataWellProdLocal(k)%ImposedFlowrate - qw
         else
            call CommonMPI_abort("Production well index should be 'p' (pressure mode) or 'f' (flowrate mode)!")
         end if

      end do

   end subroutine Residu_add_flux_contributions_wells

   subroutine Residu_RelativeNorm_local_conservation(ResConvLocal) &
      bind(C, name="Residu_RelativeNorm_local_conservation")

      type(CTVector), intent(out) :: ResConvLocal
      integer :: k

      ResConvLocal%values(:) = 0.d0

      do k = 1, NbNodeOwn_Ncpus(commRank + 1) ! node
         ResConvLocal%values(:) = ResConvLocal%values(:) + abs(ResiduNode(:, k))
      enddo

      do k = 1, NbFracOwn_Ncpus(commRank + 1) ! frac
         ResConvLocal%values(:) = ResConvLocal%values(:) + abs(ResiduFrac(:, k))
      enddo

      do k = 1, NbCellOwn_Ncpus(commRank + 1) ! cell
         ResConvLocal%values(:) = ResConvLocal%values(:) + abs(ResiduCell(:, k))
      enddo

      do k = 1, NbWellInjOwn_Ncpus(commRank + 1) ! well inj
         ResConvLocal%values(1) = ResConvLocal%values(1) + abs(ResiduWellInj(k))
      enddo

      do k = 1, NbWellProdOwn_Ncpus(commRank + 1) ! well prod
         ResConvLocal%values(1) = ResConvLocal%values(1) + abs(ResiduWellProd(k))
      enddo

   end subroutine Residu_RelativeNorm_local_conservation

   function Residu_RelativeNorm_local_closure() &
      result(ResClosLocal) &
      bind(C, name="Residu_RelativeNorm_local_closure")

      real(c_double) :: ResClosLocal
      integer :: i, k, ic, start

      ResClosLocal = 0.d0

      do k = 1, NbNodeOwn_Ncpus(commRank + 1) ! node
         ic = IncNode(k)%ic

         start = NbPhasePresente_ctx(ic)
         do i = 1, NbEqEquilibre_ctx(ic)
            ResClosLocal = ResClosLocal + abs(SmFNode(start + i, k))
         enddo
      enddo

      do k = 1, NbFracOwn_Ncpus(commRank + 1) ! frac
         ic = IncFrac(k)%ic

         start = NbPhasePresente_ctx(ic)
         do i = 1, NbEqEquilibre_ctx(ic)
            ResClosLocal = ResClosLocal + abs(SmFFrac(start + i, k))
         enddo
      enddo

      do k = 1, NbCellOwn_Ncpus(commRank + 1) ! cell
         ic = IncCell(k)%ic

         start = NbPhasePresente_ctx(ic)
         do i = 1, NbEqEquilibre_ctx(ic)
            ResClosLocal = ResClosLocal + abs(SmFCell(start + i, k))
         enddo
      enddo

   end function Residu_RelativeNorm_local_closure

   subroutine Residu_RelativeNorm_reference_conservation(Delta_t, ResConvRefLocal) &
      bind(C, name="Residu_RelativeNorm_reference_conservation")

      real(c_double), intent(in), value :: Delta_t
      type(CTVector), intent(out) :: ResConvRefLocal
      integer :: k

      ResConvRefLocal%values(:) = 1.d0 ! CHECKME: why ?

      do k = 1, NbNodeOwn_Ncpus(commRank + 1)
         ResConvRefLocal%values(:) = ResConvRefLocal%values(:) + &
                                     abs(IncNode(k)%AccVol(:))/Delta_t/1000.d0 ! CHECKME: why ? magic number!
      end do
      do k = 1, NbFracOwn_Ncpus(commRank + 1)
         ResConvRefLocal%values(:) = ResConvRefLocal%values(:) + &
                                     abs(IncFrac(k)%AccVol(:))/Delta_t/1000.d0 ! CHECKME: why ? magic number!
      end do
      do k = 1, NbCellOwn_Ncpus(commRank + 1)
         ResConvRefLocal%values(:) = ResConvRefLocal%values(:) + &
                                     abs(IncCell(k)%AccVol(:))/Delta_t/1000.d0 ! CHECKME: why ? magic number!
      end do

   end subroutine Residu_RelativeNorm_reference_conservation

   subroutine Residu_RelativeNorm_initial_conservation(Delta_t, ResConvLocal, ResConvInit) &
      bind(C, name="Residu_RelativeNorm_initial_conservation")

      real(c_double), intent(in), value :: Delta_t
      type(CTVector), intent(in) :: ResConvLocal
      type(CTVector), intent(out) :: ResConvInit
      type(CTVector) :: ResConvRefLocal
      type(CTVector) :: ResConvInitLocal
      integer :: i, Ierr

      call Residu_RelativeNorm_reference_conservation(Delta_t, ResConvRefLocal)
      do i = 1, NbCompThermique
         ResConvInitLocal%values(i) = max(ResConvLocal%values(i), ResConvRefLocal%values(i))
      end do

      call MPI_Allreduce( &
         ResConvInitLocal%values, ResConvInit%values, NbCompThermique, &
         MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

   end subroutine Residu_RelativeNorm_initial_conservation

   function Residu_RelativeNorm_initial_closure(ResClosLocal) &
      result(ResClosInit) &
      bind(C, name="Residu_RelativeNorm_initial_closure")

      real(c_double), intent(in), value :: ResClosLocal
      real(c_double) :: ResClosInit
      integer :: Ierr
      real(c_double) :: ResClosInitLocal

      ResClosInitLocal = 1.d0
      ResClosInitLocal = max(ResClosLocal, ResClosInitLocal)

      call MPI_Allreduce(ResClosInitLocal, ResClosInit, 1, &
                         MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

   end function Residu_RelativeNorm_initial_closure

   function Residu_compute_relative_norm( &
      ResConvInit, ResConvLocal, ResClosInit, ResClosLocal) &
      result(relative_norm) &
      bind(C, name="Residu_compute_relative_norm")

      type(CTVector), intent(in) :: ResConvInit, ResConvLocal
      real(c_double), intent(in), value :: ResClosInit, ResClosLocal
      real(c_double) :: relative_norm
      integer :: i, Ierr
      real(c_double) :: local

      local = ResClosLocal/ResClosInit
      do i = 1, NbCompThermique
         local = max(ResConvLocal%values(i)/ResConvInit%values(i), local)
      enddo
      call MPI_Allreduce( &
         local, relative_norm, 1, &
         MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)

   end function Residu_compute_relative_norm

   subroutine Residu_RelativeNorm(Iter, Delta_t, &
                                  ResNormRel, ResConvInit, ResClosInit)

      integer(c_int), intent(in), value :: Iter ! Newton iteration: Iter
      real(c_double), intent(in), value :: Delta_t
      real(c_double), intent(out) :: ResNormRel
      ! Residu relative norm from first Newton iteration used for all Newton iterations
      type(CTVector), intent(inout) :: ResConvInit
      real(c_double), intent(inout) :: ResClosInit

      type(CTVector) :: ResConvLocal
      real(c_double) :: ResClosLocal

      call Residu_RelativeNorm_local_conservation(ResConvLocal)
      if (Iter == 1) &
         call Residu_RelativeNorm_initial_conservation( &
         Delta_t, ResConvLocal, ResConvInit)

      ResClosLocal = Residu_RelativeNorm_local_closure()
      if (Iter == 1) &
         ResClosInit = Residu_RelativeNorm_initial_closure(ResClosLocal)

      ResNormRel = Residu_compute_relative_norm( &
                   ResConvInit, ResConvLocal, ResClosInit, ResClosLocal)

   end subroutine Residu_RelativeNorm

end module Residu
