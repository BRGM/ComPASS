!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module Residu

   use, intrinsic :: iso_c_binding, only: c_double, c_f_pointer

   use CommonMPI, only: commRank, CommonMPI_abort
   use InteroperabilityStructures, only: cpp_array_wrapper

   use DefModel, only: &
      NbPhase, NbComp, NbCompThermique, IndThermique, &
      LIQUID_PHASE, MCP, &
      NbEqEquilibreMax, NbIncTotalPrimMax, &
      NbIncPTCMax, NbIncTotalMax, NbEqFermetureMax, &
      NbPhasePresente_ctx, NumPhasePresente_ctx

   use LoisThermoHydro, only: &
#ifdef _WIP_FREEFLOW_STRUCTURES_
      FreeFlowMolarFlowrateCompNode, &
      FreeFlowHmCompNode, &
      FreeFlowHTTemperatureNetRadiationNode, &
      FreeFlowMolarFlowrateEnthalpieNode, &
      AtmEnthalpieNode, &
#endif
      DensitemolaireKrViscoCompWellInj, &
      DensitemolaireKrViscoEnthalpieWellInj, &
      DensitemolaireKrViscoCompNode, &
      DensitemolaireKrViscoCompCell, &
      DensitemolaireKrViscoCompFrac, &
      DensitemolaireKrViscoEnthalpieNode, &
      DensitemolaireKrViscoEnthalpieCell, &
      DensitemolaireKrViscoEnthalpieFrac, &
      DensitemolaireSatComp, &
      DensitemolaireEnergieInterneSat

   use Physics, only: CpRoche, atm_comp, rain_flux
   use IncPrimSecd, only: SmFNode, SmFCell, SmFFrac
   use Flux, only: FluxDarcyKI, FluxFourierKI, FluxDarcyFI, FluxFourierFI

   use NeumannContribution, only: NodeNeumannBC

   use NumbyContext, only: &
      NbEqEquilibre_ctx, NumCompCtilde_ctx, NbCompCtilde_ctx
   use IncCVReservoir, only: &
      NbCellLocal_Ncpus, NbNodeLocal_Ncpus, NbFracLocal_Ncpus, &
      IncAll, IncNode, IncCell, IncFrac
   use IncCVWells, only: &
      PerfoWellInj, DataWellInjLocal, &
      NbWellInjLocal_Ncpus, &
      NodebyWellProdLocal, PerfoWellProd, IncPressionWellInj, &
      IncPressionWellProd, NodeDatabyWellProdLocal
   use VAGFrac, only: &
      ThermalSourceVol, &
      PoroVolDarcy, &
      PoroVolFourier, &
      Poro_1VolFourier
   use MeshSchema, only: &
      DataWellProdLocal, NodebyCellLocal, FracbyCellLocal, &
      FaceToFracLocal, FracToFaceLocal, NodebyFaceLocal, &
#ifdef _WIP_FREEFLOW_STRUCTURES_
      IdFFNodeLocal, &
#endif
      IdNodeLocal, NbNodeOwn_Ncpus, NbCellOwn_Ncpus, NbFracOwn_Ncpus, &
      NodebyWellInjLocal, NodeDatabyWellInjLocal, NbWellProdLocal_Ncpus, &
      SurfFreeFlowLocal

   implicit none

   ! Define a residual for all values (DOFFamily array)
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
   ! FIXME: explicit name like AccVolCell_previous or AccVolCell_0 would be better
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
      Residu_AccVol, &
      Residu_associate_pointers

   private :: &
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

   !> \brief Clear ResiduCell/Frac/Node/Well
   subroutine Residu_clear_residuals()

      ResiduCell(:, :) = 0.d0
      ResiduFrac(:, :) = 0.d0
      ResiduNode(:, :) = 0.d0
      ResiduWellInj(:) = 0.d0
      ResiduWellProd(:) = 0.d0

   end subroutine Residu_clear_residuals


   !> \brief Newton is initialized with the unknown at time step n-1
   !! What about Ctilde ?????? shoud be 0 at the beginning of the time step
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

   !> \brief Compute the residu and store it into inc%AccVol
   subroutine Residu_compute(Delta_t) &
      bind(C, name="Residu_compute")

      real(c_double), intent(in), value :: Delta_t

      integer :: i, k

      !> Clear ResiduCell/Frac/...
      call Residu_clear_residuals

      call Residu_AccVol

      !> Initialise Residu with terms of accumulation
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
      ResiduCell(NbComp + 1, :) = ResiduCell(NbComp + 1, :) - ThermalSourceVol%cells
      ResiduFrac(NbComp + 1, :) = ResiduFrac(NbComp + 1, :) - ThermalSourceVol%fractures
      ResiduNode(NbComp + 1, :) = ResiduNode(NbComp + 1, :) - ThermalSourceVol%nodes
#endif

      ! Residu for dir node will be reset as 0 in the end of this subroutine
      ! We do not consider dirichlet nodes during computing the residu such that
      ! the code is clearer

      call Residu_add_flux_contributions

      call Residu_reset_Dirichlet_nodes

      call Residu_add_Neumann_contributions

   end subroutine Residu_compute

   !> \brief Clear the residu of Dirichlet nodes 
   !! (distinction between P and T)
   subroutine Residu_reset_Dirichlet_nodes

      integer :: k

      ! Residu for dirichlet node should be reset to 0 for computing norm of Residu
      do k = 1, NbNodeLocal_Ncpus(commRank + 1)
         if (IdNodeLocal(k)%P == "d") ResiduNode(1:NbComp, k) = 0.d0
#ifdef _THERMIQUE_
         if (IdNodeLocal(k)%T == "d") ResiduNode(NbCompThermique, k) = 0.d0
#endif
      end do

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

   !> \brief Reset residu for components which are not Ctilde
   !> \todo FIXME: Could be simpler if we multiply accumulations by a mask with 0 and 1 (MCP...)
   subroutine Residu_clear_absent_components_accumulation(ic, accumulations)

      integer, intent(in) :: ic ! context
      real(c_double), intent(inout) :: accumulations(NbCompThermique) ! unknowns
      integer :: i, icp
      real(c_double) :: copy(NbComp) ! unknowns

      do i = 1, NbCompCtilde_ctx(ic)
         icp = NumCompCtilde_ctx(i, ic)
         copy(i) = accumulations(icp)
      end do

      ! CHECKME: is size different from NbCompThermique: inc(:) = 0.d0 would be enough
      accumulations(1:NbCompThermique) = 0.d0

      do i = 1, NbCompCtilde_ctx(ic)
         icp = NumCompCtilde_ctx(i, ic)
         accumulations(icp) = copy(i)
      end do

   end subroutine Residu_clear_absent_components_accumulation

   !> \brief Compute AccVol for the present phase and component
   !! Keep previous value for Ctilde
   subroutine Residu_AccVol() &
      bind(C, name="Residu_update_accumulation")

      integer :: k, m, mph, ic, icp

      ! Loop over all degrees of freedom (nodes, fractures, cells)
      do k = 1, size(IncAll)
         ic = IncAll(k)%ic
         call Residu_clear_absent_components_accumulation(ic, IncAll(k)%AccVol)
         do m = 1, NbPhasePresente_ctx(ic) ! Q_k
            mph = NumPhasePresente_ctx(m, ic)
            do icp = 1, NbComp
               ! FIXME: Could be simpler if we multiply accumulations by a mask with 0 and 1 (MCP...)
               if (MCP(icp, mph) == 1) then ! Q_k \cap P_i
                  IncAll(k)%AccVol(icp) = IncAll(k)%AccVol(icp) &
                                           + PoroVolDarcy%values(k)*DensiteMolaireSatComp%values(icp, m, k)
               end if
            end do
#ifdef _THERMIQUE_
            IncAll(k)%AccVol(NbComp + 1) = IncAll(k)%AccVol(NbComp + 1) &
                                            + PoroVolFourier%values(k)*DensiteMolaireEnergieInterneSat%values(m, k)
#endif
         end do ! end of phase m
#ifdef _THERMIQUE_
         IncAll(k)%AccVol(NbComp + 1) = IncAll(k)%AccVol(NbComp + 1) &
                                         + Poro_1VolFourier%values(k)*CpRoche*IncAll(k)%Temperature
#endif
      end do

   end subroutine Residu_AccVol

   subroutine Residu_add_flux_contributions

      integer :: k, s, fs, fk, nums, m, mph, icp
      integer :: NbNodeCell, NbFracCell, NbNodeFrac

      double precision :: Flux_ks(NbComp), FluxT_ks, DarcyFlux

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
                                       + DensiteMolaireKrViscoCompCell(icp, m, k)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieCell(m, k) &
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
                                       + DensiteMolaireKrViscoCompNode(icp, m, nums)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieNode(m, nums) &
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
                                       + DensiteMolaireKrViscoCompCell(icp, m, k)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieCell(m, k) &
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
                                       + DensiteMolaireKrViscoCompFrac(icp, m, nums)*DarcyFlux

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieFrac(m, nums) &
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
                                       + DensiteMolaireKrViscoCompFrac(icp, m, k)*FluxDarcyFI(mph, s, k)

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieFrac(m, k) &
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
                                       + DensiteMolaireKrViscoCompNode(icp, m, nums)*FluxDarcyFI(mph, s, k)

                     end if
                  end do

#ifdef _THERMIQUE_

                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieNode(m, nums) &
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

#ifdef _WIP_FREEFLOW_STRUCTURES_
      call Residu_add_flux_contributions_FF_node
#endif

   end subroutine Residu_add_flux_contributions

   subroutine Residu_add_flux_contributions_wells

      integer :: k, s, nums, m, mph, icp

      double precision :: Flux_ks(NbComp), FluxT_ks
      double precision :: Pws, Tws, Ps, Ts, WIDws, WIFws, qw, qe, Ps_Pws
      logical  something_is_produced,  something_is_injected

      ! Injection well
      ! It is possible that one ghost well contains some own nodes,
      ! thus, when assembly of flux, we need to take into account
      ! not only own well, but also ghost well
      do k = 1, NbWellInjLocal_Ncpus(commRank + 1) ! injection well k

         ! Check if the well is closed
         if (DataWellInjLocal(k)%IndWell == 'c') then
            ResiduWellInj(k) = 0.d0
            cycle ! FIXME: if Fourier contribion is considered we should take it into account
         end if

         qw = 0.d0
         qe = 0.d0
         something_is_injected = .false.

         ! nodes of well k
         do s = NodebyWellInjLocal%Pt(k) + 1, NodebyWellInjLocal%Pt(k + 1)
            nums = NodebyWellInjLocal%Num(s)

            Pws = PerfoWellInj(s)%Pression ! P_{w,s}
            Tws = PerfoWellInj(s)%Temperature ! T_{w,s}
            Ps = IncNode(nums)%Pression ! P_s DCQ-TODO: Add capilarry pressure for 2-phase
            Ts = IncNode(nums)%Temperature ! T_s

            WIDws = NodeDatabyWellInjLocal%Val(s)%WID
            WIFws = NodeDatabyWellInjLocal%Val(s)%WIF

            ! Flux q_{w,s,i} and q_{w,s,e}
            Ps_Pws = Ps - Pws
            if (Ps_Pws < 0.d0)  then
               something_is_injected = .true.
               Flux_ks = DensiteMolaireKrViscoCompWellInj(:, s)*WIDws*Ps_Pws
#ifdef _THERMIQUE_
               FluxT_ks = DensiteMolaireKrViscoEnthalpieWellInj(s)*WIDws*Ps_Pws
#endif
               ! node equation
               ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) + Flux_ks

#ifdef _THERMIQUE_
               ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) + FluxT_ks ! + WIFws * (Ts-Tws)
#endif
               qw = qw + sum(Flux_ks)
               qe = qe + FluxT_ks
            end if
         end do ! end of node s in injection well k

         DataWellInjLocal(k)%actual_mass_flowrate = qw
         DataWellInjLocal(k)%actual_energy_flowrate = qe

         ! inj well equation
         if (DataWellInjLocal(k)%IndWell == 'p') then
            ResiduWellInj(k) = DataWellInjLocal(k)%PressionMax - IncPressionWellInj(k)
         else if (DataWellInjLocal(k)%IndWell == 'f') then
            if (something_is_injected) then
               ResiduWellInj(k) = qw - DataWellInjLocal(k)%ImposedFlowrate
            else
               ResiduWellInj(k) = 0.d0
            end if
         else
            call CommonMPI_abort("Injection well index should be 'p' (pressure mode) or 'f' (flowrate mode)!")
         end if

      end do ! end of injection well k

      ! Production well
      do k = 1, NbWellProdLocal_Ncpus(commRank + 1)

         ! Check if the well is closed
         if (DataWellProdLocal(k)%IndWell == 'c') then
            ResiduWellProd(k) = 0.d0
            cycle ! FIXME: if Fourier contribion is considered we should take it into account
         end if

         qw = 0.d0
         qe = 0.d0
         something_is_produced = .false.

         ! nodes of well k
         do s = NodebyWellProdLocal%Pt(k) + 1, NodebyWellProdLocal%Pt(k + 1)
            nums = NodebyWellProdLocal%Num(s)
            Pws = PerfoWellProd(s)%Pression ! P_{w,s}
            Ps = IncNode(nums)%Pression ! P_s
            Ps_Pws = Ps - Pws
            WIDws = NodeDatabyWellProdLocal%Val(s)%WID
            if (Ps_Pws > 0.d0) then
               something_is_produced = .true.
               Flux_ks(:) = 0.d0
               FluxT_ks = 0.d0
               do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
                  mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)
                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! \cap P_i
                        Flux_ks(icp) = Flux_ks(icp) + DensiteMolaireKrViscoCompNode(icp, m, nums)*WIDws*Ps_Pws
                     end if
                  end do
#ifdef _THERMIQUE_
                  FluxT_ks = FluxT_ks + DensiteMolaireKrViscoEnthalpieNode(m, nums)*WIDws*Ps_Pws
#endif
               end do
               ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) + Flux_ks(:)
#ifdef _THERMIQUE_
               ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) + FluxT_ks
#endif
               qw = qw + sum(Flux_ks)
               qe = qe + FluxT_ks
            end if
         end do ! end of node s in production well k

         DataWellProdLocal(k)%actual_mass_flowrate = qw
         DataWellProdLocal(k)%actual_energy_flowrate = qe

         ! prod well equation
         if (DataWellProdLocal(k)%IndWell == 'p') then
            ResiduWellProd(k) = IncPressionWellProd(k) - DataWellProdLocal(k)%PressionMin
         else if (DataWellProdLocal(k)%IndWell == 'f') then
            if (something_is_produced) then
               ResiduWellProd(k) = DataWellProdLocal(k)%ImposedFlowrate - qw
            else
               ResiduWellProd(k) = 0.d0
            end if
         else
            call CommonMPI_abort("Production well index should be 'p' (pressure mode) or 'f' (flowrate mode)!")
         end if

      end do ! end of production well k

   end subroutine Residu_add_flux_contributions_wells

#ifdef _WIP_FREEFLOW_STRUCTURES_
   subroutine Residu_add_flux_contributions_FF_node

    integer :: nums, m, mph, icp
    double precision :: Flux_FreeFlow(NbComp), FluxT_FreeFlow, FreeFlowMolarFlowrate

    do nums=1, NbNodeOwn_Ncpus(commRank+1)

      if(IdFFNodeLocal(nums)) then ! loop over freeflow dof only, avoid reservoir node

         Flux_FreeFlow(:) = 0.d0
         FluxT_FreeFlow = 0.d0

         do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic)
            mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

            FreeFlowMolarFlowrate = IncNode(nums)%FreeFlow_flowrate(mph)

            if(FreeFlowMolarFlowrate>=0.d0) then

               do icp=1, NbComp
                  if(MCP(icp,mph)==1) then ! \cap P_i

                     Flux_FreeFlow(icp) = Flux_FreeFlow(icp) + SurfFreeFlowLocal(nums) * ( &
                             FreeFlowMolarFlowrateCompNode(icp,m,nums) &
                           + FreeFlowHmCompNode(icp,m,nums) &
                           + atm_comp(icp,mph) * rain_flux(mph) )! rain source term (rain_flux(gas)=0)
                  end if
               end do ! end of icp
#ifdef _THERMIQUE_
               FluxT_FreeFlow = FluxT_FreeFlow + SurfFreeFlowLocal(nums) * ( &
                       FreeFlowMolarFlowrateEnthalpieNode(m,nums) &
                     + AtmEnthalpieNode(m,nums) * rain_flux(mph) ) ! rain source term (rain_flux(gas)=0)
#endif

            else ! FreeFlowMolarFlowrate<0.d0
            ! liq phase never enters in this loop because always FreeFlow_flowrate(liq)>=0.d0
         
               do icp=1, NbComp
                  if(MCP(icp,mph)==1) then ! \cap P_i
         
                     Flux_FreeFlow(icp) = Flux_FreeFlow(icp) + SurfFreeFlowLocal(nums) * ( &
                                FreeFlowMolarFlowrate * atm_comp(icp,mph) &
                              + FreeFlowHmCompNode(icp,m,nums) &
                              + atm_comp(icp,mph) * rain_flux(mph) ) ! rain source term (rain_flux(gas)=0)
                  end if
               end do ! end of icp
#ifdef _THERMIQUE_
               FluxT_FreeFlow = FluxT_FreeFlow + SurfFreeFlowLocal(nums) * ( &
                       FreeFlowMolarFlowrate * AtmEnthalpieNode(m,nums) &
                     + AtmEnthalpieNode(m,nums) * rain_flux(mph) ) ! rain source term (rain_flux(gas)=0)
#endif
            endif ! sign of FreeFlowMolarFlowrate
         enddo ! m

#ifdef _THERMIQUE_
         FluxT_FreeFlow = FluxT_FreeFlow &
               + SurfFreeFlowLocal(nums) * FreeFlowHTTemperatureNetRadiationNode(nums)
         ResiduNode(NbComp + 1, nums) = ResiduNode(NbComp + 1, nums) + FluxT_FreeFlow
#endif

         ResiduNode(1:NbComp, nums) = ResiduNode(1:NbComp, nums) + Flux_FreeFlow(1:NbComp)

      endif ! avoid reservoir node

   enddo ! node nums

   end subroutine Residu_add_flux_contributions_FF_node
#endif

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

end module Residu
