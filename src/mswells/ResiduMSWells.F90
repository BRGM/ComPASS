module ResiduMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use DefModel
   use mpi
   use CommonMPI
   use CommonType
   use MeshSchema
   use IncCVReservoir
   use NumbyContext
   use LoisThermoHydro
   use LoisThermoHydroMSWells
   use IncCVMSWells
   use VSHydroMSWells
   use MeshSchemaMSWells
   use MSWellsData
#else
   use iso_c_binding, only: c_double, c_bool, c_int, c_f_pointer

   use DefModel, only: &
      NbPhase, NbComp, NbContexte, NbCompThermique, MCP
#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      GAS_PHASE, LIQUID_PHASE
#endif

   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use CommonType, only: &
      CSRArray2dble, CSR, CommonType_deallocCSR

   use MeshSchema, only: &
      XNodeLocal, &
      IdNodeLocal, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus, &
      NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbMSWellNodeLocal_Ncpus

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, NumPhasePresente_ctx, NbPhasePresente_ctx

   use NumbyContext, only: &
      NbCompCtilde_ctx, NumCompCtilde_ctx, NbIncPTC_ctx, &
      NbIncTotalPrim_ctx, NbIncTotalPrim_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   use IncCVMSWells, only: &
      IncMSWell

   use LoisThermoHydro, only: &
      ContextInfo, &
      LoisThermoHydro_init_cv

   use LoisThermoHydroMSWells, only: &
      LoisThermoHydroMSWell, TYPE_LoisThermoHydro

   use VSHydroMSWells, only: &
      VSHydroMSWell, &
      VSHydroHeadMSWell

   use MeshSchemaMSWells, only: &
      NbMSWellLocal, NbMSWellNodeOwn, &
      NbMSWellNodeLocal, &
      IncIdxNodebyMSWellLocal

   use InteroperabilityStructures, only: cpp_array_wrapper
   use MSWellsData, only: MSWellDataNodebyMSWell

#endif
   implicit none

   ! Residu for injection  and production wells
   real(c_double), pointer, dimension(:, :), public:: &
      ResiduMSWell      !Indexed with the new local-mswell-indexing

   ! AccMole of time
   real(c_double), pointer, dimension(:, :), public:: &
      AccVolMSWell_1   !Indexed with the new local-mswell-indexing

   ! Qw by MSWell to manage BDCs at the head
   real(c_double), pointer, dimension(:), public:: &
      QwByMSWell       !Value per mswell producer

   logical :: initialize_acc_terms !Should be true when the Newton Algorithm starts

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Some convenient stuff

#ifdef ComPASS_DIPHASIC_CONTEXT
   !This implementation works  only for two phases: LIQUID and GAS
   integer, parameter, private :: MSWellNumPhase(NbPhase) = (/GAS_PHASE, LIQUID_PHASE/)
#endif

   public :: &
      ResiduMSWells_allocate, &
      ResiduMSWells_free, &
      ResiduMSWells_compute, &
      ResiduMSWells_reset_history, &
      ResiduMSWells_associate_pointers, &
      ResiduMSWells_AccVol, &
      ResiduMSWells_compute_prod

contains

   subroutine ResiduMSWells_associate_pointers( &
      mswell_nodes, mswell_nodes_accumulation) &
      bind(C, name="ResiduMSWells_associate_pointers")

      type(cpp_array_wrapper), intent(in), value :: &
         mswell_nodes, mswell_nodes_accumulation

      if (mswell_nodes%n /= NbMSWellNodeLocal_Ncpus(commRank + 1)*NbCompThermique) &
         call CommonMPI_abort('inconsistent mswell_nodes residual size')
      if (mswell_nodes_accumulation%n /= NbMSWellNodeLocal_Ncpus(commRank + 1)*NbCompThermique) &
         call CommonMPI_abort('inconsistent accumulation mswell_nodes size')

      call c_f_pointer(mswell_nodes%p, ResiduMSWell, [NbCompThermique, NbMSWellNodeLocal_Ncpus(commRank + 1)])
      call c_f_pointer(mswell_nodes_accumulation%p, AccVolMSWell_1, [NbCompThermique, NbMSWellNodeLocal_Ncpus(commRank + 1)])

   end subroutine ResiduMSWells_associate_pointers

   subroutine ResiduMSWells_allocate

      allocate (QwByMSWell(NbMSWellLocal_Ncpus(commRank + 1)))
      initialize_acc_terms = .true.
      !initialize_acc_terms = .false. !Useful when we  reinitialize from a file

   end subroutine ResiduMSWells_allocate

   !> \brief Reset residu for components which are not Ctilde
   !> \todo FIXME: Could be simpler if we multiply accumulations by a mask with 0 and 1 (MCP...)
   !> TODO: This function is the same as the one in Residu, put it somewhere else
   subroutine ResiduMSWells_clear_absent_components_accumulation(ic, accumulations)

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

   end subroutine ResiduMSWells_clear_absent_components_accumulation

   !> \brief Compute the residu and store it into inc%AccVol
   subroutine ResiduMSWells_compute(Delta_t) &
      bind(C, name="ResiduMSWells_compute")

      real(c_double), intent(in), value :: Delta_t

      !Clear residuals and aux data
      ResiduMSWell(:, :) = 0.d0
      QwByMSWell(:) = 0.d0

      !Compute Acc Term
      call ResiduMSWells_AccVol

      !Compute mswells producer stuff
      call ResiduMSWells_compute_prod(Delta_t)

   end subroutine ResiduMSWells_compute

   !> \brief Init AccVol_1 with inc%AccVol
   !> What about Ctilde ?????? shoud be 0 at the beginning of the time step
   subroutine ResiduMSWells_reset_history() &
      bind(C, name="ResiduMSWells_reset_history")

      integer :: k, s, num_s, unk_idx_s, icp, mph, m
      type(TYPE_IncCVReservoir) :: coats_s
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! This blocl is useful when we reinit from a file, if not comment this block
      !!if (initialize_acc_terms) then
      !!   initialize_acc_terms = .false.
      !!   call ResiduMSWells_AccVol !Update accumulation term, useful when init ResiduMSWell from a file
      !!endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call ResiduMSWells_AccVol !Update accumulation term

      !For each producer
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !We begin with the accumulation terms
      do k = 1, NbMSWellLocal

         !For all nodes at the well, from queue to head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) !numbering in mesh

            if (IdNodeLocal(num_s)%Proc == "g") cycle

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx

            coats_s = IncMSWell(s)%coats
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! loop of component

            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i
                     !Use new indexing
                     AccVolMSWell_1(icp, unk_idx_s) = coats_s%AccVol(icp)

                  end if !MCP

               end do !icp
#ifdef _THERMIQUE_
               !Use new indexing
               AccVolMSWell_1(NbComp + 1, unk_idx_s) = coats_s%AccVol(NbComp + 1)
#endif
            end do !m
         end do ! end of loop s

      end do! end of loop NbMSWellLocal

   end subroutine ResiduMSWells_reset_history

   !> \brief Compute AccVol for the present phase and component
   !! Keep previous value for Ctilde
   subroutine ResiduMSWells_AccVol() &
      bind(C, name="ResiduMSWells_update_accumulation")

      !tmp
      integer :: k, s, num_s, unk_idx_s, icp, mph, m
      type(TYPE_IncCVReservoir) :: coats_s
      type(TYPE_LoisThermoHydro) lthydro_s

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Note: Not mandatory to put the ResiduMSWell equal to zero for ghost nodes,
      !       but it is very useful for debugging in parallel

      !For each open mswell
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !We begin with the accumulation terms for all mswells
      do k = 1, NbMSWellLocal

         if (DataMSWellLocal(k)%IndWell == 'c') cycle

         !For all nodes at the well, from queue to head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) !numbering in mesh

            if (IdNodeLocal(num_s)%Proc == "g") cycle !ghost node

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx

            lthydro_s = LoisThermoHydroMSWell(s)

            call ResiduMSWells_clear_absent_components_accumulation(IncMSWell(s)%coats%ic, IncMSWell(s)%coats%AccVol)
            coats_s = IncMSWell(s)%coats
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! loop of component

            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                     coats_s%AccVol(icp) = coats_s%AccVol(icp) &
                                           + (lthydro_s%DensiteMolaire(mph) &
                                              *coats_s%Saturation(mph) &
                                              *coats_s%Comp(icp, mph) &
                                              )*MSWellDataNodebyMSWell(s)%vol

                  end if !MCP

               end do !icp
#ifdef _THERMIQUE_

               coats_s%AccVol(NbComp + 1) = coats_s%AccVol(NbComp + 1) &
                                            + (lthydro_s%DensiteMolaire(mph) &
                                               *coats_s%Saturation(mph) &
                                               *lthydro_s%Enthalpie(mph) &
                                               )*MSWellDataNodebyMSWell(s)%vol

#endif
            end do !m

            IncMSWell(s)%coats = coats_s !save state

         end do ! end of loop s
      end do ! end of loop k

   end subroutine ResiduMSWells_AccVol

   subroutine ResiduMSWells_compute_prod(Delta_t)

      real(c_double), intent(in), value :: Delta_t
      !tmp
      integer :: k, s, sp, num_s, num_sp, unk_idx_s, i, icp, mph, m
      integer :: upwind_idx(NbPhase), unk_idx_sp, uidx, well_head_idx
      type(TYPE_IncCVReservoir) :: coats_s, coats_sp, coats_up
      type(TYPE_LoisThermoHydro) lthydro_s
      double precision :: ecross_sec, xk_s, xk_sp, esize, flux_mph, flux_mph_icp, flux_mph_e
      double precision :: VsbyPhase(NbPhase), e_thcond, fluxF
      logical ::  phase_in_up_context
#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort("Multi-segmented wells are only implemented for diphasic physics!")

#else

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Note: Not mandatory to put the ResiduMSWell equal to zero for ghost nodes,
      !       but it is very useful for debugging in parallel

      !We begin with the accumulation terms for all mswells
      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         !For all nodes at the well, from queue to head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) !numbering in mesh

            if (IdNodeLocal(num_s)%Proc == "g") cycle !ghost node

            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx

            coats_s = IncMSWell(s)%coats
            lthydro_s = LoisThermoHydroMSWell(s)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! loop of component

            do icp = 1, NbComp
               ResiduMSWell(icp, unk_idx_s) = ResiduMSWell(icp, unk_idx_s) &
                                              + (coats_s%AccVol(icp) - AccVolMSWell_1(icp, unk_idx_s))/Delta_t

            end do !icp
#ifdef _THERMIQUE_

            ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                  + (coats_s%AccVol(NbComp + 1) - AccVolMSWell_1(NbComp + 1, unk_idx_s))/Delta_t
#endif

         end do ! end of loop s
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Qw for BC at the well head node
         s = NodebyMSWellLocal%Pt(k + 1) !head idx at well

         unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx
         coats_s = IncMSWell(s)%coats

         do icp = 1, NbComp

            QwByMSWell(k) = QwByMSWell(k) &
                            - (coats_s%AccVol(icp) - AccVolMSWell_1(icp, unk_idx_s))/Delta_t
         end do

      end do! end of loop NbMSWellLocal

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Assemble by edges
      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         ! Save the well_head_idx
         well_head_idx = NodebyMSWellLocal%Pt(k + 1)
         !For all nodes at the well except the head
         !s is the well node numbering

         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            sp = NodeDatabyMSWellLocal%Val(s)%PtParent ! parent of s (well numbering)
            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the parent in the reservoir

            xk_s = XNodeLocal(1, num_s)    ! x-coordinate of node s
            xk_sp = XNodeLocal(1, num_sp)  ! x-coordinate of parent of s
            esize = (xk_s - xk_sp)**2

            xk_s = XNodeLocal(2, num_s)    ! y-coordinate of node s
            xk_sp = XNodeLocal(2, num_sp)   ! y-coordinate of parent of s
            esize = esize + (xk_s - xk_sp)**2

            xk_s = XNodeLocal(3, num_s)    ! z-coordinate of node s
            xk_sp = XNodeLocal(3, num_sp)   ! z-coordinate of parent of s

            esize = esize + (xk_s - xk_sp)**2
            esize = dsqrt(esize) !Edge-Size

            ecross_sec = MSWellDataNodebyMSWell(s)%edata%cross_section

            coats_s = IncMSWell(s)%coats
            coats_sp = IncMSWell(sp)%coats
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Prepare indices to be used for the Residual
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown  idx for son node
            unk_idx_sp = IncIdxNodebyMSWellLocal%Num(sp)!unkown idx for parent node

#ifdef _THERMIQUE_
            e_thcond = MSWellDataNodebyMSWell(s)%thcond(1) !TODO: Get the the avg thermal conductivity using Sats
#endif

            !We get first the upwind index for each phase
            do m = 1, NbPhase
               mph = MSWellNumPhase(m)

               !Recall: the son has the edge data
               VsbyPhase(mph) = VSHydroMSWell(s)%edata%SVel(mph)

               if (VsbyPhase(mph) .ge. 0.d0) then
                  upwind_idx(mph) = s
               else
                  upwind_idx(mph) = sp !thus up_idx is the idx at the local mesh
               endif

            end do !m

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Molar fluxes
            do m = 1, NbPhase

               mph = MSWellNumPhase(m)
               uidx = upwind_idx(mph)

               coats_up = IncMSWell(uidx)%coats !Upwind

               !Check if  mph is in the current context of coats_up
               phase_in_up_context = .false.
               do i = 1, NbPhasePresente_ctx(coats_up%ic)
                  if (mph == NumPhasePresente_ctx(i, coats_s%ic)) then
                     phase_in_up_context = .true.
                     exit
                  endif
               end do

               if (.not. phase_in_up_context) cycle

               flux_mph = VSHydroMSWell(s)%edata%Svel(mph) &
                          *ecross_sec &
                          *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph)

               do icp = 1, NbComp

                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                     flux_mph_icp = flux_mph &
                                    *IncMSWell(uidx)%coats%Comp(icp, mph)

                     if (IdNodeLocal(num_s)%Proc == "o") then
                        ResiduMSWell(icp, unk_idx_s) = ResiduMSWell(icp, unk_idx_s) &
                                                       + flux_mph_icp
                     end if

                     if (IdNodeLocal(num_sp)%Proc == "o") then
                        ResiduMSWell(icp, unk_idx_sp) = ResiduMSWell(icp, unk_idx_sp) &
                                                        - flux_mph_icp
                     end if

                     !Qw for BC at the well head node
                     if (sp == well_head_idx) then
                        QwByMSWell(k) = QwByMSWell(k) + flux_mph_icp
                     end if

                  end if
               end do !icp

#ifdef _THERMIQUE_

               flux_mph_e = flux_mph &
                            *LoisThermoHydroMSWell(uidx)%Enthalpie(mph)
               if (IdNodeLocal(num_s)%Proc == "o") then

                  ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                        + flux_mph_e
               endif
               if (IdNodeLocal(num_sp)%Proc == "o") then

                  ResiduMSWell(NbComp + 1, unk_idx_sp) = ResiduMSWell(NbComp + 1, unk_idx_sp) &
                                                         - flux_mph_e
               endif
#endif

            end do !m
#ifdef _THERMIQUE_
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Fourier flux
            fluxF = (coats_s%Temperature - coats_sp%Temperature) &
                    *e_thcond*ecross_sec/esize

            if (IdNodeLocal(num_s)%Proc == "o") then

               ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                     + fluxF

            end if
            if (IdNodeLocal(num_sp)%Proc == "o") then

               ResiduMSWell(NbComp + 1, unk_idx_sp) = ResiduMSWell(NbComp + 1, unk_idx_sp) &
                                                      - fluxF
            endif
#endif

         end do !s

      end do !k
#endif
   end subroutine ResiduMSWells_compute_prod

   subroutine ResiduMSWells_free

      deallocate (QwByMSWell)

   end subroutine ResiduMSWells_free

end module ResiduMSWells
