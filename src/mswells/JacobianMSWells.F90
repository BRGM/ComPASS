module JacobianMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use DefModel
   use mpi
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort
   use CommonType
   use InteroperabilityStructures
   use MeshSchema
   use IncCVReservoir
   use NumbyContext
   use LoisThermoHydro
   use LoisThermoHydroMSWells
   use IncCVMSWells
   use VSHydroMSWells
   use MeshSchemaMSWells
   use ResiduMSWells
   use LeafMSWells
   use MSWellsData
#else
   use iso_c_binding, only: c_double, c_bool, c_int

   use DefModel, only: &
      NbPhase, NbComp, NbContexte, NbCompThermique, NbIncTotalPrimMax, MCP, &
      aligmat
#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      GAS_PHASE, LIQUID_PHASE, DIPHASIC_CONTEXT
#endif

   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use CommonType, only: &
      CSRArray2dble, CSR, CommonType_deallocCSR

   use InteroperabilityStructures, only: &
      csr_block_matrix_wrapper, retrieve_csr_block_matrix, &
      cpp_array_wrapper_dim2, retrieve_double_array_dim2

   use MeshSchema, only: &
      XNodeLocal, &
      IdNodeLocal, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus, &
      NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbMSWellNodeOwn_Ncpus, NbMSWellNodeLocal_Ncpus

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, NumPhasePresente_ctx, NbPhasePresente_ctx

   use NumbyContext, only: &
      NbCompCtilde_ctx, NumCompCtilde_ctx, NbIncPTC_ctx, &
      NbIncTotalPrim_ctx, NbIncTotalPrim_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx

   use IncCVMSWells, only: &
      IncMSWell, IncCVMSWells_compute_coupling

   use LoisThermoHydro, only: &
      ContextInfo, &
      LoisThermoHydro_init_cv

   use LoisThermoHydroMSWells, only: &
      LoisThermoHydroMSWell, &
      TYPE_LoisThermoHydro

   use VSHydroMSWells, only: &
      VSHydroMSWell, &
      VSHydroHeadMSWell, &
      VSHydroMSWells_compute_bdcs_prod_well, &
      TYPE_VSHydroMSWells, &
      TYPE_VSHydroMSWellsHead

   use MeshSchemaMSWells, only: &
      NbMSWellLocal, &
      NbMSWellNodeOwn, &
      NbMSWellNodeLocal, &
      IncIdxNodebyMSWellLocal

   use ResiduMSWells, only: &
      ResiduMSWell, &
      QwByMSWell

   use LeafMSWells, only: &
      LeafNodebyMSWellLocal, LeafDataMSWell
   use MSWellsData, only: &
      MSWellDataNodebyMSWell

#endif
   implicit none

   type(CSRArray2dble), public, target :: JacA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set this to one to print mswells monitoring info
! See function JacobianMSWells_print_IP_info_to_file
#define _DEBUG_JAC_IP_MSWELLS_ 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _DEBUG_JAC_MSWELLS_ 0

#if  _DEBUG_JAC_MSWELLS_ ==1
   type(CSRArray2dble), public, target :: JacANoAlig
#endif

   ! JacA is a sparse matrix.
   ! And for a node in a mswell, we need to know the colnums which are linked
   ! to the row of the node and to the row of the parent in order to fill JacA
   ! The following vectors are used for this purpose
   integer, allocatable, dimension(:), private :: csrS, csrSP

   ! right hand side
   real(c_double), pointer, dimension(:, :), public :: &
      Sm

   ! right hand side  for the well head node
   real(c_double), pointer, dimension(:), public :: &
      SmQw

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Some convenient stuff
#ifdef ComPASS_DIPHASIC_CONTEXT
   !This implementation works  only for two phases: LIQUID and GAS
   integer, parameter, private :: MSWellNumPhase(NbPhase) = (/GAS_PHASE, LIQUID_PHASE/)
#endif

   public :: &
      JacobianMSWells_StrucJacA, &    !< init non-zero structure of Jacobian, Acc and Residuals
      JacobianMSWells_JacA_Sm, &    !< Jacobian and right hand side
      JacobianMSWells_free, &
      JacobianMSWells_print_LA_info_to_file, &
      JacobianMSWells_print_IP_info_to_file

   private :: &
      JacobianMSWells_JacA_Sm_prod, &
      JacobianMSWells_SetBDCs_prod, &
      JacobianMSWells_JacA_Sm_init_from_residual, &
      JacobianMSWells_Alignment_man, &
      JacobianMSWells_Alignment_man_row
contains

   subroutine retrieve_jacobian_mswells(A) &
      bind(C, name="retrieve_jacobian_mswells")
      type(csr_block_matrix_wrapper), intent(inout) :: A

      call retrieve_csr_block_matrix(JacA, A)

   end subroutine retrieve_jacobian_mswells

   subroutine retrieve_right_hand_side_mswells(a) &
      bind(C, name="retrieve_right_hand_side_mswells")
      type(cpp_array_wrapper_dim2), intent(inout) :: a

      call retrieve_double_array_dim2(Sm, a)

   end subroutine retrieve_right_hand_side_mswells

   subroutine JacobianMSWells_StrucJacA

      integer :: s, k, num_s, num_s_p, unk_idx_s, unk_idx_s_p
      integer :: Nz, ncol, j, i
      !Maximum number of  connections per node inside a well, change if needed
      integer, parameter ::NbNzbyLineMax = 100

      integer, dimension(:), allocatable :: &
         NbNnzbyline ! number of non zeros of each line
      integer, dimension(:, :), allocatable :: &
         NumNzbyLine ! columns of non zeros for each line

      !TODO: Safety check--- Make sure we have only two phases: GAS and LIQUID

      !Let us set the new  unknowns ordering for producers to be used for the Jacobian

      allocate (csrS(NbMSWellNodeLocal))
      allocate (csrSP(NbMSWellNodeLocal))
      csrS(:) = 0
      csrSP(:) = 0

      !Start creating the Jacobian structure
      JacA%Nb = NbMSWellNodeOwn !number of mswell node own
      allocate (JacA%Pt(JacA%Nb + 1))
      allocate (NbNnzbyLine(JacA%Nb))
      allocate (NumNzbyLine(JacA%Nb, NbNzbyLineMax))

      !** We first set the diagonal terms **
      do k = 1, NbMSWellLocal
         ! looping from  queue to  the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s)
            !Unkown index
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)

            if (IdNodeLocal(num_s)%Proc == "o") then ! node own
               NbNnzbyLine(unk_idx_s) = 1 !with itself
               !Fill the column
               NumNzbyLine(unk_idx_s, NbNnzbyLine(unk_idx_s)) = unk_idx_s
            end if

         end do
      end do

      !**Coupling with other well nodes**
      !  Doing the loop for each edge
      do k = 1, NbMSWellLocal
         ! looping from  queue to the node before the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            num_s = NodebyMSWellLocal%Num(s)
            num_s_p = NodeDatabyMSWellLocal%Val(s)%Parent ! parent  num in the local mesh
            !Unkown indices of node s and its parent
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)
            unk_idx_s_p = IncIdxNodebyMSWellLocal%Num(NodeDatabyMSWellLocal%Val(s)%PtParent)

            if (IdNodeLocal(num_s)%Proc == "o") then ! node own
               NbNnzbyLine(unk_idx_s) = NbNnzbyLine(unk_idx_s) + 1 !with its parent
               !Fill the column, TODO: Check array limit
               NumNzbyLine(unk_idx_s, NbNnzbyLine(unk_idx_s)) = unk_idx_s_p

            end if

            if (IdNodeLocal(num_s_p)%Proc == "o") then ! node own
               NbNnzbyLine(unk_idx_s_p) = NbNnzbyLine(unk_idx_s_p) + 1 !with its son
               !Fill the column, TODO: Check array limit
               NumNzbyLine(unk_idx_s_p, NbNnzbyLine(unk_idx_s_p)) = unk_idx_s

            end if

         end do

      end do

      Nz = 0
      JacA%Pt(1) = 0
      do i = 1, JacA%Nb
         Nz = Nz + NbNnzbyLine(i)
         JacA%Pt(i + 1) = Nz
      end do

      allocate (JacA%Num(Nz))
      JacA%Num(:) = 0
      ncol = 1
      !fill the columns in the JacA
      do s = 1, NbMSWellNodeOwn
         unk_idx_s = s
         do j = 1, NbNnzbyLine(s)
            JacA%Num(ncol) = NumNzbyLine(unk_idx_s, j)
            ncol = ncol + 1
         end do
      end do

      ! allocate JacA%Val
      allocate (JacA%Val(NbCompThermique, NbCompThermique, Nz)) ! number of non zero
#if  _DEBUG_JAC_MSWELLS_ == 1
      allocate (JacANoAlig%Val(NbCompThermique, NbCompThermique, Nz)) ! number of non zero
#endif
      ! allocate Sm and SmQw
      if (JacA%Nb .eq. 0) then
         allocate (Sm(NbCompThermique, 1)) !InteroperabilityStructures.F90 needs to have a dummy Sm
      else
         allocate (Sm(NbCompThermique, JacA%Nb))
      endif

      allocate (SmQw(NbMSWellLocal_Ncpus(commRank + 1)))

      deallocate (NbNnzbyLine)
      deallocate (NumNzbyLine)

      if (commRank == 0) then !Master proc
#if  _DEBUG_JAC_MSWELLS_ == 1
         open (unit=110, file='reswell_vec.val', status='unknown', action='write')
         write (110, *) ""
         close (110)

         open (unit=110, file='sm_vec.val', status='unknown', action='write')
         write (110, *) ""
         close (110)

         open (unit=110, file='jacobian.val', status='unknown', action='write')
         write (110, *) ""
         close (110)

         open (unit=110, file='jacobian_ali.val', status='unknown', action='write')
         write (110, *) ""
         close (110)
#endif
#if _DEBUG_JAC_IP_MSWELLS_ == 1

         open (unit=110, file='qwwell.val', status='unknown', action='write')
         write (110, *) ""
         close (110)

         open (unit=110, file='pwwell.val', status='unknown', action='write')
         write (110, *) ""
         close (110)

         open (unit=110, file='QPwell.val', status='unknown', action='write')
         write (110, *) ""
         close (110)
#endif
      endif
   end subroutine JacobianMSWells_StrucJacA

   !> \brief Fill Sm with the residual values
   !!
   !!   If the node is dir, it is set in JacobianMSWells_SetBDCs_prod
   !!
   subroutine JacobianMSWells_JacA_Sm_init_from_residual()
      integer :: k, row_s, s, unk_idx_s, num_s

      do k = 1, NbMSWellLocal
         ! looping from  queue to the  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s)
            !Unkown indices of node s
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)

            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            row_s = unk_idx_s! row of s
            Sm(:, row_s) = -ResiduMSWell(:, row_s)
         end do
      end do

   end subroutine JacobianMSWells_JacA_Sm_init_from_residual

   !> \brief fill Jacobian and right hand side: main subroutine
   subroutine JacobianMSWells_JacA_Sm(Delta_t, alignment_flag) &
      bind(C, name="JacobianMSWells_ComputeJacSm")

      real(c_double), intent(in), value :: Delta_t
      logical(c_bool), intent(in), value:: alignment_flag

      JacA%Val(:, :, :) = 0.d0
      Sm(:, :) = 0.d0
      SmQw(:) = 0.d0

      call JacobianMSWells_JacA_Sm_init_from_residual

      call JacobianMSWells_JacA_Sm_prod(Delta_t)

      if (alignment_flag) then
         ! Alignment of Jacobian
#if  _DEBUG_JAC_MSWELLS_ == 1
         JacANoAlig%Val(:, :, :) = JacA%Val(:, :, :)
#endif
         call JacobianMSWells_Alignment_man !call manual method
      endif

   end subroutine JacobianMSWells_JacA_Sm

   !> \brief fill Jacobian and right hand side for mswell producers
   subroutine JacobianMSWells_JacA_Sm_prod(Delta_t)

      double precision, intent(in) :: Delta_t
      !tmp
      integer :: k, s, sp, num_s, num_sp, unk_idx_s, Nz, i, j, icp, mph, m
      integer :: upwind_idx(NbPhase), unk_idx_sp, uidx, well_head_idx
      integer :: Nz_ss, Nz_ssp, Nz_sps, Nz_spsp, Nz_sup, Nz_spup
      type(TYPE_IncCVReservoir) :: coats_s, coats_sp, coats_up
      type(TYPE_LoisThermoHydro) :: lthydro_s
      double precision :: A_s, A_sp, A_up, ecross_sec, xk_s, xk_sp, esize
      double precision :: Sm_s, Sm_sp, Sm_s_sp, Sm_up
      double precision :: VsbyPhase(NbPhase), e_thcond
      logical :: is_node_s_own, is_node_sp_own, phase_in_up_context

#ifndef ComPASS_DIPHASIC_CONTEXT
      if (NbMSWellLocal > 0) then
         call CommonMPI_abort( &
            "In function JacobianMSWells_JacA_Sm_prod: Multi-segmented wells are only implemented for diphasic physics!")
      endif
#else

      !For each producer

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !We begin assembling the diagonal terms
      do k = 1, NbMSWellLocal
         !Only mswells producers
         if (DataMSWellLocal(k)%Well_type == 'i') cycle

         !Save the well_head_idx
         well_head_idx = NodebyMSWellLocal%Pt(k + 1)

         !For all nodes at the well, from queue to head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s) !numbering in mesh
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx
            sp = NodeDatabyMSWellLocal%Val(s)%PtParent ! parent of s (well numbering)

            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            !Nz: Diagonal idx in Jac (Recall we built the diagonals terms first for each row)
            !TODO: we can search for Nz manually to make it more robust
            Nz = JacA%Pt(unk_idx_s) + 1

            coats_s = IncMSWell(s)%coats
            lthydro_s = LoisThermoHydroMSWell(s)

            if (DataMSWellLocal(k)%IndWell == 'c') then  !mswell is closed
               !Put the Identity matrix
               do i = 1, NbComp
                  JacA%Val(i, i, Nz) = 1.d0
               end do

#ifdef _THERMIQUE_
               JacA%Val(NbComp + 1, NbComp + 1, Nz) = 1.d0

#endif
               cycle !continue s-loop
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! loop of components in Ctilde, n_k is an unknown independent

            do i = 1, NbCompCtilde_ctx(coats_s%ic)
               icp = NumCompCtilde_ctx(i, coats_s%ic)

               j = NbIncTotalPrim_ctx(coats_s%ic) + i
               JacA%val(j, icp, Nz) = MSWellDataNodebyMSWell(s)%vol/Delta_t
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! loop of component, n_k is not an unknown independent

            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                     do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                        A_s = &
                           ( &
                           lthydro_s%divDensiteMolaire(j, mph) &
                           *coats_s%Saturation(mph) &
                           *coats_s%Comp(icp, mph) &
                           + lthydro_s%DensiteMolaire(mph) &
                           *lthydro_s%divSaturation(j, mph) &
                           *coats_s%Comp(icp, mph) &
                           + lthydro_s%DensiteMolaire(mph) &
                           *coats_s%Saturation(mph) &
                           *lthydro_s%divComp(j, icp, mph) &
                           )*MSWellDataNodebyMSWell(s)%vol

                        JacA%Val(j, icp, Nz) = JacA%Val(j, icp, Nz) + A_s/Delta_t

                     end do !j

                     Sm_s = &
                        ( &
                        lthydro_s%SmDensiteMolaire(mph) &
                        *coats_s%Saturation(mph) &
                        *coats_s%Comp(icp, mph) &
                        + lthydro_s%DensiteMolaire(mph) &
                        *coats_s%Saturation(mph) &
                        *lthydro_s%SmComp(icp, mph) &
                        )*MSWellDataNodebyMSWell(s)%vol

                     Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) + Sm_s/Delta_t

                     !For BC at well head node
                     if (sp == well_head_idx) then
                        SmQw(k) = SmQw(k) - Sm_s/Delta_t
                     endif

                  end if !MCP

               end do !icp
#ifdef _THERMIQUE_

               do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                  A_s = &
                     (lthydro_s%divDensiteMolaire(j, mph) &
                      *coats_s%Saturation(mph) &
                      *lthydro_s%Enthalpie(mph) &
                      + lthydro_s%DensiteMolaire(mph) &
                      *lthydro_s%divSaturation(j, mph) &
                      *lthydro_s%Enthalpie(mph) &
                      + lthydro_s%DensiteMolaire(mph) &
                      *coats_s%Saturation(mph) &
                      *lthydro_s%divEnthalpie(j, mph) &
                      )*MSWellDataNodebyMSWell(s)%vol

                  JacA%Val(j, NbComp + 1, Nz) = JacA%Val(j, NbComp + 1, Nz) + A_s/Delta_t
               end do !j

               Sm_s = &
                  (lthydro_s%SmDensiteMolaire(mph) &
                   *coats_s%Saturation(mph) &
                   *lthydro_s%Enthalpie(mph) &
                   + lthydro_s%DensiteMolaire(mph) &
                   *coats_s%Saturation(mph) &
                   *lthydro_s%SmEnthalpie(mph) &
                   )*MSWellDataNodebyMSWell(s)%vol

               Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) + Sm_s/Delta_t

#endif
            end do !m
         end do ! end of loop s

      end do! end of loop NbMSWellLocal

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Assemble by edges
      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         !Save the well_head_idx
         well_head_idx = NodebyMSWellLocal%Pt(k + 1)

         !For all nodes at the well except the head
         !s is the well node numbering
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            sp = NodeDatabyMSWellLocal%Val(s)%PtParent ! parent of s (well numbering)
            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the parent in the reservoir

            xk_s = XNodeLocal(1, num_s)    ! x-cordinate of node s
            xk_sp = XNodeLocal(1, num_sp)  ! x-cordinate of parent of s
            esize = (xk_s - xk_sp)**2

            xk_s = XNodeLocal(2, num_s)    ! y-cordinate of node s
            xk_sp = XNodeLocal(2, num_sp)   ! y-cordinate of parent of s
            esize = esize + (xk_s - xk_sp)**2

            xk_s = XNodeLocal(3, num_s)    ! z-cordinate of node s
            xk_sp = XNodeLocal(3, num_sp)   ! z-cordinate of parent of s

            esize = esize + (xk_s - xk_sp)**2
            esize = dsqrt(esize) !Edge-Size

            ecross_sec = MSWellDataNodebyMSWell(s)%edata%cross_section

            coats_s = IncMSWell(s)%coats
            coats_sp = IncMSWell(sp)%coats

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
            !Prepare indices to be used for the Jacobian
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown  idx for son node
            unk_idx_sp = IncIdxNodebyMSWellLocal%Num(sp)!unkown idx for parent node

            is_node_s_own = .false.
            is_node_sp_own = .false.

            if (IdNodeLocal(num_s)%Proc == "o") is_node_s_own = .true.
            if (IdNodeLocal(num_sp)%Proc == "o") is_node_sp_own = .true.

            if (is_node_s_own) then
               !!For all columns at the row of node s
               do m = JacA%Pt(unk_idx_s) + 1, JacA%Pt(unk_idx_s + 1)
                  csrS(JacA%Num(m)) = m
               end do

               Nz_ss = csrS(unk_idx_s)
               Nz_ssp = csrS(unk_idx_sp)
            endif

            if (is_node_sp_own) then
               !For all columns at the row of the parent node of s
               do m = JacA%Pt(unk_idx_sp) + 1, JacA%Pt(unk_idx_sp + 1)
                  csrSP(JacA%Num(m)) = m
               end do
               Nz_spsp = csrSP(unk_idx_sp)
               Nz_sps = csrSP(unk_idx_s)
            endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Compute Jacobian Fluxes
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

               !Upwind indices to be used for the JacA
               Nz_sup = csrS(IncIdxNodebyMSWellLocal%Num(uidx))
               Nz_spup = csrSP(IncIdxNodebyMSWellLocal%Num(uidx))

               do icp = 1, NbComp

                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !Derivatives with respect Vs taken at node_s
                     do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                        A_s = &
                           VSHydroMSWell(s)%edata%divSvel_ss(j, mph, 1) &
                           *ecross_sec &
                           *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                           *IncMSWell(uidx)%coats%Comp(icp, mph)

                        if (is_node_s_own) then
                           JacA%Val(j, icp, Nz_ss) = JacA%Val(j, icp, Nz_ss) + A_s
                        endif

                        if (is_node_sp_own) then
                           JacA%Val(j, icp, Nz_sps) = JacA%Val(j, icp, Nz_sps) - A_s
                        endif
                     end do !j

                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !Derivatives with respect Vs taken at node_sp
                     do j = 1, NbIncTotalPrim_ctx(coats_sp%ic)

                        A_sp = &
                           VSHydroMSWell(s)%edata%divSvel_ss(j, mph, 2) &
                           *ecross_sec &
                           *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                           *IncMSWell(uidx)%coats%Comp(icp, mph)

                        if (is_node_s_own) then
                           JacA%Val(j, icp, Nz_ssp) = JacA%Val(j, icp, Nz_ssp) + A_sp
                        endif

                        if (is_node_sp_own) then
                           JacA%Val(j, icp, Nz_spsp) = JacA%Val(j, icp, Nz_spsp) - A_sp
                        endif

                     end do !j

                     !Vector Sm: Derivative with respect Vs taken at both nodes
                     Sm_s_sp = VSHydroMSWell(s)%edata%SmSvel(mph) &
                               *ecross_sec &
                               *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                               *IncMSWell(uidx)%coats%Comp(icp, mph)

                     if (is_node_s_own) &
                        Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) + Sm_s_sp
                     if (is_node_sp_own) &
                        Sm(icp, unk_idx_sp) = Sm(icp, unk_idx_sp) - Sm_s_sp

                     !For BC at well head node
                     if (sp == well_head_idx) then
                        SmQw(k) = SmQw(k) + Sm_s_sp
                     endif

                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !Derivatives with respect density and Comp at node upw
                     do j = 1, NbIncTotalPrim_ctx(coats_up%ic)

                        A_up = &
                           LoisThermoHydroMSWell(uidx)%divDensiteMolaire(j, mph) &
                           *IncMSWell(uidx)%coats%Comp(icp, mph) &
                           + &
                           LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                           *LoisThermoHydroMSWell(uidx)%divComp(j, icp, mph)

                        A_up = VsbyPhase(mph) &
                               *ecross_sec &
                               *A_up

                        if (is_node_s_own) then
                           JacA%Val(j, icp, Nz_sup) = JacA%Val(j, icp, Nz_sup) + A_up
                        endif

                        if (is_node_sp_own) then
                           JacA%Val(j, icp, Nz_spup) = JacA%Val(j, icp, Nz_spup) - A_up
                        endif

                     end do !j

                     Sm_up = &
                        LoisThermoHydroMSWell(uidx)%SmDensiteMolaire(mph) &
                        *IncMSWell(uidx)%coats%Comp(icp, mph) &
                        + &
                        LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                        *LoisThermoHydroMSWell(uidx)%SmComp(icp, mph)

                     Sm_up = VsbyPhase(mph) &
                             *ecross_sec &
                             *Sm_up
                     if (is_node_s_own) &
                        Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) + Sm_up
                     if (is_node_sp_own) &
                        Sm(icp, unk_idx_sp) = Sm(icp, unk_idx_sp) - Sm_up

                     !For BC at well head node
                     if (sp == well_head_idx) then
                        SmQw(k) = SmQw(k) + Sm_up
                     endif

                  end if !MCP

               end do !icp

#ifdef _THERMIQUE_

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Derivatives with respect Vs at node s
               do j = 1, NbIncTotalPrim_ctx(coats_s%ic)
                  A_s = &
                     VSHydroMSWell(s)%edata%divSvel_ss(j, mph, 1) &
                     *ecross_sec &
                     *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                     *LoisThermoHydroMSWell(uidx)%Enthalpie(mph)

                  if (is_node_s_own) then
                     JacA%Val(j, NbComp + 1, Nz_ss) = JacA%Val(j, NbComp + 1, Nz_ss) + A_s
                  endif

                  if (is_node_sp_own) then
                     JacA%Val(j, NbComp + 1, Nz_sps) = JacA%Val(j, NbComp + 1, Nz_sps) - A_s
                  endif
               end do !j

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Derivatives with respect Vs at node sp
               do j = 1, NbIncTotalPrim_ctx(coats_sp%ic)

                  A_sp = &
                     VSHydroMSWell(s)%edata%divSvel_ss(j, mph, 2) &
                     *ecross_sec &
                     *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                     *LoisThermoHydroMSWell(uidx)%Enthalpie(mph)

                  if (is_node_s_own) then
                     JacA%Val(j, NbComp + 1, Nz_ssp) = JacA%Val(j, NbComp + 1, Nz_ssp) + A_sp
                  endif

                  if (is_node_sp_own) then
                     JacA%Val(j, NbComp + 1, Nz_spsp) = JacA%Val(j, NbComp + 1, Nz_spsp) - A_sp
                  endif

               end do !j

               !Vector Sm: Derivative with respect Vs taken at both nodes
               Sm_s_sp = VSHydroMSWell(s)%edata%SmSvel(mph) &
                         *ecross_sec &
                         *LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                         *LoisThermoHydroMSWell(uidx)%Enthalpie(mph)
               if (is_node_s_own) &
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) + Sm_s_sp
               if (is_node_sp_own) &
                  Sm(NbComp + 1, unk_idx_sp) = Sm(NbComp + 1, unk_idx_sp) - Sm_s_sp

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !Derivatives with respect density and enthalpy at node upw
               do j = 1, NbIncTotalPrim_ctx(coats_up%ic)
                  A_up = &
                     LoisThermoHydroMSWell(uidx)%divDensiteMolaire(j, mph) &
                     *LoisThermoHydroMSWell(uidx)%Enthalpie(mph) &
                     + &
                     LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                     *LoisThermoHydroMSWell(uidx)%divEnthalpie(j, mph)

                  A_up = VsbyPhase(mph) &
                         *ecross_sec &
                         *A_up

                  if (is_node_s_own) then
                     JacA%Val(j, NbComp + 1, Nz_sup) = JacA%Val(j, NbComp + 1, Nz_sup) + A_up
                  endif

                  if (is_node_sp_own) then
                     JacA%Val(j, NbComp + 1, Nz_spup) = JacA%Val(j, NbComp + 1, Nz_spup) - A_up
                  endif

               end do !j

               Sm_up = &
                  LoisThermoHydroMSWell(uidx)%SmDensiteMolaire(mph) &
                  *LoisThermoHydroMSWell(uidx)%Enthalpie(mph) &
                  + &
                  LoisThermoHydroMSWell(uidx)%DensiteMolaire(mph) &
                  *LoisThermoHydroMSWell(uidx)%SmEnthalpie(mph)

               Sm_up = VsbyPhase(mph) &
                       *ecross_sec &
                       *Sm_up
               if (is_node_s_own) &
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) + Sm_up
               if (is_node_sp_own) &
                  Sm(NbComp + 1, unk_idx_sp) = Sm(NbComp + 1, unk_idx_sp) - Sm_up

            end do !m

            !Derivatives of Fourier Fluxes at node s
            do j = 1, NbIncTotalPrim_ctx(coats_s%ic)
               A_s = LoisThermoHydroMSWell(s)%divTemperature(j) &
                     *e_thcond*ecross_sec/esize

               if (is_node_s_own) then
                  JacA%Val(j, NbComp + 1, Nz_ss) = JacA%Val(j, NbComp + 1, Nz_ss) + A_s
               endif

               if (is_node_sp_own) then
                  JacA%Val(j, NbComp + 1, Nz_sps) = JacA%Val(j, NbComp + 1, Nz_sps) - A_s
               endif
            end do !j

            Sm_s = LoisThermoHydroMSWell(s)%SmTemperature &
                   *e_thcond*ecross_sec/esize

            if (is_node_s_own) &
               Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) + Sm_s

            !Derivatives of Fourier Fluxes at node sp
            do j = 1, NbIncTotalPrim_ctx(coats_sp%ic)
               A_sp = -LoisThermoHydroMSWell(sp)%divTemperature(j) &
                      *e_thcond*ecross_sec/esize

               if (is_node_s_own) then
                  JacA%Val(j, NbComp + 1, Nz_ssp) = JacA%Val(j, NbComp + 1, Nz_ssp) + A_sp
               endif

               if (is_node_sp_own) then
                  JacA%Val(j, NbComp + 1, Nz_spsp) = JacA%Val(j, NbComp + 1, Nz_spsp) - A_sp
               endif
            end do    !j

            Sm_sp = -LoisThermoHydroMSWell(sp)%SmTemperature &
                    *e_thcond*ecross_sec/esize

            if (is_node_s_own) &
               Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) + Sm_sp

#endif

         end do ! end of loop s

      end do! end of loop  NbMSWellLocal

      call JacobianMSWells_SetBDCs_prod
#endif

   end subroutine JacobianMSWells_JacA_Sm_prod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     The following subroutine puts the BDCs as  follows:
!!!
!!!     CL au head node
!!!
!!!
!!!     flux sortant et pression imposee a Pout (outer pressure)
!!!
!!!     zero Fourier flux
!!!
!!!     Systeme a resoudre au head node (si diphasique)
!!!
!!!         Vsl = Um-Vsg
!!!
!!!         Vsg = Fg(sgk,sgk,Um)
!!!
!!!         Vsl*zetal + Vsg*zetag = qw  (molar conservation) -> Um  fct de qw et variables thermo
!!!
!!!         P_s = Pout
!!!
!!!         cons energy
!!!
!!!     On elimine Vsg, Vsl et Um
!!!
!!!
!!!     Restera les equations primaires
!!!
!!!          Pk = pwmin  (eq 1) +  molar cons eq for icp=2,...,NbComp et energy cons eq
!!!
!!!          qw = qwmax -> molar cons eq for icp=1,...,NbComp et energy cons eq
!!!
!!!    Note: This subroutine modifies ResiduMSWell
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
   subroutine JacobianMSWells_SetBDCs_prod

      !tmp
      integer :: k, s, num_s, unk_idx_s, col, icp, j, m, mph
      integer :: ncols_s, JacRow_end_idx_s, JacRow_begin_idx_s
      double precision ::  ecross_sec, qw, SmQw_s, Sm_s
      double precision, allocatable :: divQwRow(:, :), divVsHeadRow(:, :, :)
      type(TYPE_IncCVReservoir) :: coats_s
      type(TYPE_VSHydroMSWells)     :: vshydro_s
      type(TYPE_VSHydroMSWellsHead) :: vshydro_h
      type(TYPE_LoisThermoHydro)    :: lthydro_s
      double precision ::  VsHeadbyPhase(NbPhase), SmVsHeadbyPhase(NbPhase)
      double precision ::  bdc_T, bdc_Comp
      integer :: mswell_idx_s, leaf_data_idx_s, Nz
      logical :: is_node_s_own

      !For each producer

      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         s = NodebyMSWellLocal%Pt(k + 1)  !Well head node idx
         coats_s = IncMSWell(s)%coats
         num_s = NodebyMSWellLocal%Num(s) !numbering in mesh
         unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown idx
         is_node_s_own = .false.
         if (IdNodeLocal(num_s)%Proc == "o") is_node_s_own = .true.

         if (.NOT. is_node_s_own) then ! node not own
            ResiduMSWell(:, unk_idx_s) = 0.d0
            cycle
         end if

         !First of all update the status of the producer well

         if (DataMSWellLocal(k)%IndWell == 'p') then  !MSWell-Producer is in Fix-Pressure-Mode
            if (QwByMSWell(k) .gt. DataMSWellLocal(k)%ImposedFlowrate) then
               DataMSWellLocal(k)%IndWell = 'f'
            end if
         end if

         if (DataMSWellLocal(k)%IndWell == 'f') then  !MSWell-Producer is in Fix-Pressure-Mode
            if (coats_s%Pression .lt. DataMSWellLocal(k)%PressionMin) then
               DataMSWellLocal(k)%IndWell = 'p'
            end if
         end if

         ecross_sec = MSWellDataNodebyMSWell(s)%edata%cross_section
         lthydro_s = LoisThermoHydroMSWell(s)

         JacRow_begin_idx_s = JacA%Pt(unk_idx_s) + 1 !First column at the row of node s
         JacRow_end_idx_s = JacA%Pt(unk_idx_s + 1)   !Last column at the row of node s
         ncols_s = JacRow_end_idx_s - JacRow_begin_idx_s + 1

         allocate (divVsHeadRow(NbIncTotalPrimMax, NbPhase, ncols_s))

         if (DataMSWellLocal(k)%IndWell == 'p') then  !MSWell-Producer is in Fix-Pressure-Mode

            !qw = - \sum_icp ResWell(k,icp)/secw  -> molar flow rate at well node
            qw = QwByMSWell(k)/ecross_sec
            SmQw_s = SmQw(k)/ecross_sec

            allocate (divQwRow(NbIncTotalPrimMax, ncols_s))
            do col = JacRow_begin_idx_s, JacRow_end_idx_s

               divQwRow(:, col - JacRow_begin_idx_s + 1) = 0.d0

               do icp = 1, NbComp
                  do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives

                     divQwRow(j, col - JacRow_begin_idx_s + 1) = divQwRow(j, col - JacRow_begin_idx_s + 1) &
                                                                 - (1.d0/ecross_sec)*JacA%Val(j, icp, col)
                  end do !j
               end do !icp
            enddo !col

            !Compute Vs (superficial vel), div and Sm vectors
            !   for each phase at the head node
            !   using  qw and thermo variables at head node,
            !   Recall VSHydroMSWell(s)%edata%divSvel_ss(:,:,1)
            !      has the Vs data for the top
            call VSHydroMSWells_compute_bdcs_prod_well(k, qw)

            !Vs data
            vshydro_s = VSHydroMSWell(s)
            vshydro_h = VSHydroHeadMSWell(k)

            VsHeadbyPhase(:) = 0.d0
            divVsHeadRow(:, :, :) = 0.d0
            SmVsHeadbyPhase(:) = 0.d0

            ! Div of Vs + Sm
            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               VsHeadbyPhase(mph) = vshydro_s%edata%Svel(mph)
               !Watch it here: we have upwinded derivatives
               divVsHeadRow(:, mph, :) = divQwRow(:, :)*vshydro_h%divQwVsHead(mph)
               SmVsHeadbyPhase(mph) = vshydro_s%edata%SmSvel(mph) + vshydro_h%divQwVsHead(mph)*SmQw_s

               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                     do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives
                        divVsHeadRow(j, mph, 1) = divVsHeadRow(j, mph, 1) &
                                                  + vshydro_s%edata%divSvel_ss(j, mph, 1)
                     end do !j

                  end if !MCP
               end do !icp

               if (NbComp > 1) then
                  ! molar cons eq  icp = 1,...,NbComp
                  do icp = 1, NbComp
                     if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                        ! column which corresponds to node_s
                        ! (recall we built the diagonals terms first for each row)
                        !BDC
                        bdc_Comp = VsHeadbyPhase(mph) &
                                   *ecross_sec &
                                   *lthydro_s%DensiteMolaire(mph) &
                                   *coats_s%Comp(icp, mph)

                        ResiduMSWell(icp, unk_idx_s) = ResiduMSWell(icp, unk_idx_s) &
                                                       + bdc_comp

                        if (is_node_s_own) then

                           Sm_s = VsHeadbyPhase(mph) &
                                  *ecross_sec &
                                  *lthydro_s%SmDensiteMolaire(mph) &
                                  *coats_s%Comp(icp, mph)

                           Sm_s = Sm_s + (VsHeadbyPhase(mph) &
                                          *ecross_sec &
                                          *lthydro_s%DensiteMolaire(mph) &
                                          *lthydro_s%SmComp(icp, mph))

                           Sm_s = Sm_s + (vshydro_s%edata%SmSvel(mph) &
                                          *ecross_sec &
                                          *lthydro_s%DensiteMolaire(mph) &
                                          *coats_s%Comp(icp, mph))

                           !Fist modify Sm accordingly as ResiduMSWell
                           Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) &
                                                - bdc_comp

                           Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) + Sm_s

                           do col = JacRow_begin_idx_s, JacRow_end_idx_s

                              do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives
                                 JacA%Val(j, icp, col) = JacA%Val(j, icp, col) &
                                                         + divVsHeadRow(j, mph, col - JacRow_begin_idx_s + 1) &
                                                         *lthydro_s%DensiteMolaire(mph) &
                                                         *coats_s%Comp(icp, mph) &
                                                         *ecross_sec

                              end do!j
                           end do !col

                           do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                              !Recall the first column is the diagonal
                              col = JacRow_begin_idx_s
                              JacA%Val(j, icp, col) = JacA%Val(j, icp, col) &
                                                      + VsHeadbyPhase(mph)*ecross_sec*( &
                                                      lthydro_s%divDensiteMolaire(j, mph) &
                                                      *coats_s%Comp(icp, mph) &
                                                      + &
                                                      lthydro_s%DensiteMolaire(mph) &
                                                      *lthydro_s%divComp(j, icp, mph))

                           end do !j

                        end if !node s own

                     end if !MCP
                  end do !icp
               end if !NbComp > 1
            end do !m
#ifdef _THERMIQUE_
            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               bdc_T = VsHeadbyPhase(mph) &
                       *ecross_sec &
                       *lthydro_s%DensiteMolaire(mph) &
                       *lthydro_s%Enthalpie(mph)

               !BDC
               ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                     + bdc_T
               if (is_node_s_own) then

                  !First modify Sm accordingly as ResiduMSWell
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              - bdc_T

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + VsHeadbyPhase(mph) &
                                              *ecross_sec &
                                              *lthydro_s%SmDensiteMolaire(mph) &
                                              *lthydro_s%Enthalpie(mph)

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + vshydro_s%edata%SmSvel(mph) &
                                              *ecross_sec &
                                              *lthydro_s%DensiteMolaire(mph) &
                                              *lthydro_s%SmEnthalpie(mph)

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + SmVsHeadbyPhase(mph) &
                                              *ecross_sec &
                                              *lthydro_s%DensiteMolaire(mph) &
                                              *lthydro_s%Enthalpie(mph)

                  do col = JacRow_begin_idx_s, JacRow_end_idx_s
                     do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives
                        JacA%Val(j, NbComp + 1, col) = JacA%Val(j, NbComp + 1, col) &
                                                       + divVsHeadRow(j, mph, col - JacRow_begin_idx_s + 1) &
                                                       *ecross_sec &
                                                       *lthydro_s%DensiteMolaire(mph) &
                                                       *lthydro_s%Enthalpie(mph)
                     end do !j
                  end do !col

                  col = JacRow_begin_idx_s
                  do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                     JacA%Val(j, NbComp + 1, col) = JacA%Val(j, NbComp + 1, col) &
                                                    + VsHeadbyPhase(mph) &
                                                    *ecross_sec &
                                                    *( &
                                                    lthydro_s%divDensiteMolaire(j, mph) &
                                                    *lthydro_s%Enthalpie(mph) &
                                                    + &
                                                    lthydro_s%DensiteMolaire(mph) &
                                                    *lthydro_s%divEnthalpie(j, mph) &
                                                    )
                  end do !j
               end if ! node s own

            end do !m
#endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!eq 1: P_s-Pout = 0.d0   (total mole cons eq already eliminated -> cons of 1st comp replaced by P_s-Pout=0)
            !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! CL pression en sortie
            IncMSWell(s)%coats%Pression = DataMSWellLocal(k)%PressionMin
            ResiduMSWell(1, unk_idx_s) = 0.d0 ! P_s-Pout = 0
            if (is_node_s_own) then

               !Modify Sm accordingly
               Sm(1, unk_idx_s) = 0.d0

               !For the pressure eq., set zero  the whole row associated with node_s
               JacA%Val(:, 1, JacRow_begin_idx_s:JacRow_end_idx_s) = 0.d0
               JacA%Val(1, 1, JacRow_begin_idx_s) = 1.d0
            end if

            deallocate (divQwRow)
            !###########################################################
         else  !MSWell-Producer is in Fix-FlowRate-Mode

            !get well flowrate
            qw = DataMSWellLocal(k)%ImposedFlowrate/ecross_sec

            !Compute Vs (superficial vel), div and Sm vectors
            !   for each phase at the head node
            !   using  qw and thermo variables at head node,
            !   Recall VSHydroMSWell(s)%edata%divSvel_ss(:,:,1)
            !      has the Vs data for the top
            call VSHydroMSWells_compute_bdcs_prod_well(k, qw)

            !Vs data
            vshydro_s = VSHydroMSWell(s)
            vshydro_h = VSHydroHeadMSWell(k)

            VsHeadbyPhase(:) = 0.d0
            divVsHeadRow(:, :, :) = 0.d0
            SmVsHeadbyPhase(:) = 0.d0

            ! Div of Vs + Sm
            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               VsHeadbyPhase(mph) = vshydro_s%edata%Svel(mph)
               SmVsHeadbyPhase(mph) = vshydro_s%edata%SmSvel(mph)

               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! Q_k \cap P_i

                     do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives
                        divVsHeadRow(j, mph, 1) = divVsHeadRow(j, mph, 1) &
                                                  + vshydro_s%edata%divSvel_ss(j, mph, 1)
                     end do !j

                     ! column which corresponds to node_s
                     ! (recall we built the diagonals terms first for each row)
                     !BDC
                     bdc_Comp = VsHeadbyPhase(mph) &
                                *ecross_sec &
                                *lthydro_s%DensiteMolaire(mph) &
                                *coats_s%Comp(icp, mph)

                     ResiduMSWell(icp, unk_idx_s) = ResiduMSWell(icp, unk_idx_s) &
                                                    + bdc_comp

                     if (is_node_s_own) then

                        Sm_s = VsHeadbyPhase(mph) &
                               *ecross_sec &
                               *lthydro_s%SmDensiteMolaire(mph) &
                               *coats_s%Comp(icp, mph)

                        Sm_s = Sm_s + (VsHeadbyPhase(mph) &
                                       *ecross_sec &
                                       *lthydro_s%DensiteMolaire(mph) &
                                       *lthydro_s%SmComp(icp, mph))

                        Sm_s = Sm_s + (vshydro_s%edata%SmSvel(mph) &
                                       *ecross_sec &
                                       *lthydro_s%DensiteMolaire(mph) &
                                       *coats_s%Comp(icp, mph))

                        !Fist modify Sm accordingly as ResiduMSWell
                        Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) &
                                             - bdc_comp

                        Sm(icp, unk_idx_s) = Sm(icp, unk_idx_s) + Sm_s
                        !Recall the first column is the diagonal
                        col = JacRow_begin_idx_s

                        do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives
                           JacA%Val(j, icp, col) = JacA%Val(j, icp, col) &
                                                   + divVsHeadRow(j, mph, 1) &
                                                   *lthydro_s%DensiteMolaire(mph) &
                                                   *coats_s%Comp(icp, mph) &
                                                   *ecross_sec

                        end do !j
                        do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                           JacA%Val(j, icp, col) = JacA%Val(j, icp, col) &
                                                   + VsHeadbyPhase(mph)*ecross_sec*( &
                                                   lthydro_s%divDensiteMolaire(j, mph) &
                                                   *coats_s%Comp(icp, mph) &
                                                   + &
                                                   lthydro_s%DensiteMolaire(mph) &
                                                   *lthydro_s%divComp(j, icp, mph))

                        end do  !j
                     end if ! s_own
                  end if !MCP
               end do !icp
            end do !m

#ifdef _THERMIQUE_
            do m = 1, NbPhasePresente_ctx(coats_s%ic) ! Q_k, k is node
               mph = NumPhasePresente_ctx(m, coats_s%ic)

               bdc_T = VsHeadbyPhase(mph) &
                       *ecross_sec &
                       *lthydro_s%DensiteMolaire(mph) &
                       *lthydro_s%Enthalpie(mph)

               !BDC
               ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                     + bdc_T
               if (is_node_s_own) then

                  !First modify Sm accordingly as ResiduMSWell
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              - bdc_T

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + VsHeadbyPhase(mph) &
                                              *ecross_sec &
                                              *lthydro_s%SmDensiteMolaire(mph) &
                                              *lthydro_s%Enthalpie(mph)

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + vshydro_s%edata%SmSvel(mph) &
                                              *ecross_sec &
                                              *lthydro_s%DensiteMolaire(mph) &
                                              *lthydro_s%SmEnthalpie(mph)

                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + SmVsHeadbyPhase(mph) &
                                              *ecross_sec &
                                              *lthydro_s%DensiteMolaire(mph) &
                                              *lthydro_s%Enthalpie(mph)

                  col = JacRow_begin_idx_s
                  do j = 1, NbIncTotalPrimMax !Watch it here: we have upwinded derivatives

                     JacA%Val(j, NbComp + 1, col) = JacA%Val(j, NbComp + 1, col) &
                                                    + divVsHeadRow(j, mph, 1) &
                                                    *ecross_sec &
                                                    *lthydro_s%DensiteMolaire(mph) &
                                                    *lthydro_s%Enthalpie(mph)
                  end do !j

                  do j = 1, NbIncTotalPrim_ctx(coats_s%ic)

                     JacA%Val(j, NbComp + 1, col) = JacA%Val(j, NbComp + 1, col) &
                                                    + VsHeadbyPhase(mph) &
                                                    *ecross_sec &
                                                    *( &
                                                    lthydro_s%divDensiteMolaire(j, mph) &
                                                    *lthydro_s%Enthalpie(mph) &
                                                    + &
                                                    lthydro_s%DensiteMolaire(mph) &
                                                    *lthydro_s%divEnthalpie(j, mph) &
                                                    )
                  end do !j
               end if ! node s own

            end do !m
#endif
         end if !MSWell-Producer is in Fix-FlowRate-Mode

         deallocate (divVsHeadRow)

      end do ! k

      !If coupling with the reservoir then nothing else to be done
      if (IncCVMSWells_compute_coupling()) return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       CL de flux entrant aux leaf nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do k = 1, NbMSWellLocal

         !If the mswell is an injector or is closed, do nothing
         if ((DataMSWellLocal(k)%IndWell == 'c') .or. (DataMSWellLocal(k)%Well_type == 'i')) cycle

         !For all leaf nodes of the well
         do s = LeafNodebyMSWellLocal%Pt(k) + 1, LeafNodebyMSWellLocal%Pt(k + 1)
            mswell_idx_s = LeafNodebyMSWellLocal%Num(s) ! idx at mswell
            num_s = NodebyMSWellLocal%Num(mswell_idx_s) ! idx at the rewervoir mesh
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(mswell_idx_s)!unkown  idx for son node
            is_node_s_own = .false.
            if (IdNodeLocal(num_s)%Proc == "o") is_node_s_own = .true.

            if (.NOT. is_node_s_own) then ! node not own
               ResiduMSWell(:, unk_idx_s) = 0.d0
               cycle
            end if

            leaf_data_idx_s = LeafNodebyMSWellLocal%Val(s) ! idx at the LeafData
            ecross_sec = MSWellDataNodebyMSWell(mswell_idx_s)%edata%cross_section
            coats_s = IncMSWell(mswell_idx_s)%coats

            if (LeafDataMSWell(leaf_data_idx_s)%BDC == 'f') then ! Impose flowrate without resevoir data

               ResiduMSWell(1:NbComp, unk_idx_s) = ResiduMSWell(1:NbComp, unk_idx_s) &
                                                   - ecross_sec*LeafDataMSWell(leaf_data_idx_s)%Q(1:NbComp)
               !Modify Sm accordingly
               if (is_node_s_own) then
                  Sm(1:NbComp, unk_idx_s) = Sm(1:NbComp, unk_idx_s) &
                                            + ecross_sec*LeafDataMSWell(leaf_data_idx_s)%Q(1:NbComp)
               end if
#ifdef _THERMIQUE_
               ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                     - ecross_sec*LeafDataMSWell(leaf_data_idx_s)%Q(NbComp + 1)
               !Modify Sm accordingly as ResiduMSWell
               if (is_node_s_own) then
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + ecross_sec*LeafDataMSWell(leaf_data_idx_s)%Q(NbComp + 1)
               end if
#endif

            else if (LeafDataMSWell(leaf_data_idx_s)%BDC == 'r') then ! Impose flowrate using resevoir data

               ResiduMSWell(1:NbComp, unk_idx_s) = ResiduMSWell(1:NbComp, unk_idx_s) &
                                                   - LeafDataMSWell(leaf_data_idx_s)%TFluxRes(1:NbComp) &
                                                   *(LeafDataMSWell(leaf_data_idx_s)%PressionRes - coats_s%Pression)
               !Modify Sm accordingly as ResiduMSWell
               if (is_node_s_own) then
                  Sm(1:NbComp, unk_idx_s) = Sm(1:NbComp, unk_idx_s) &
                                            + LeafDataMSWell(leaf_data_idx_s)%TFluxRes(1:NbComp) &
                                            *(LeafDataMSWell(leaf_data_idx_s)%PressionRes - coats_s%Pression)

               endif

#ifdef _THERMIQUE_
               ResiduMSWell(NbComp + 1, unk_idx_s) = ResiduMSWell(NbComp + 1, unk_idx_s) &
                                                     - LeafDataMSWell(leaf_data_idx_s)%TFluxRes(NbComp + 1) &
                                                     *(LeafDataMSWell(leaf_data_idx_s)%PressionRes - coats_s%Pression)
               !Modify Sm accordingly as ResiduMSWell
               if (is_node_s_own) then
                  Sm(NbComp + 1, unk_idx_s) = Sm(NbComp + 1, unk_idx_s) &
                                              + LeafDataMSWell(leaf_data_idx_s)%TFluxRes(NbComp + 1) &
                                              *(LeafDataMSWell(leaf_data_idx_s)%PressionRes - coats_s%Pression)

               endif

#endif
               !Modify JacA accordingly
               if (is_node_s_own) then
                  !Nz: Diagonal idx in Jac (Recall we built the diagonals terms first for each row)
                  Nz = JacA%Pt(unk_idx_s) + 1

                  JacA%Val(1, 1:NbComp, Nz) = JacA%val(1, 1:NbComp, Nz) &
                                              + LeafDataMSWell(leaf_data_idx_s)%TFluxRes(1:NbComp)
#ifdef _THERMIQUE_
                  JacA%Val(1, NbComp + 1, Nz) = JacA%val(1, NbComp + 1, Nz) &
                                                + LeafDataMSWell(leaf_data_idx_s)%TFluxRes(NbComp + 1)
#endif
               end if

            else
               call CommonMPI_abort("Error at JacobianMSWells_SetBDCs_prod: Unkown BDC at leaf node!")

            end if !Leaf BDC

         end do

      end do

   end subroutine JacobianMSWells_SetBDCs_prod

   ! Alignment of Jacobian: manually method
   ! operation on JacA, Sm
   subroutine JacobianMSWells_Alignment_man

      integer :: k, row_s, s, unk_idx_s, num_s

      do k = 1, NbMSWellLocal
         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            num_s = NodebyMSWellLocal%Num(s)
            !Unkown indices of node s
            unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)

            if (IdNodeLocal(num_s)%Proc == "g") cycle ! node not own

            row_s = unk_idx_s! row of s
            call JacobianMSWells_Alignment_man_row(row_s, IncMSWell(s)%coats%ic)
         end do
      end do

   end subroutine JacobianMSWells_Alignment_man

   ! sub subroutine of JacobianMSWell_Alignment: manually method
   ! used for Alignment for mswell nodes
   subroutine JacobianMSWells_Alignment_man_row(rowk, ic)

      integer, intent(in) :: rowk, ic
      integer :: i

      double precision, dimension(NbCompThermique, NbCompThermique) :: &
         AA, BB
      double precision, dimension(NbCompThermique) :: &
         Smk

      ! the index order of JacA%Val(:,:,nz) is (col, row)

      BB(:, :) = aligmat(:, :, ic)

      do i = JacA%Pt(rowk) + 1, JacA%Pt(rowk + 1)

         AA(:, :) = JacA%Val(:, :, i)
         call dgemm('N', 'N', NbCompThermique, NbCompThermique, NbCompThermique, &
                    1.d0, AA, NbCompThermique, BB, NbCompThermique, 0.d0, &
                    JacA%Val(:, :, i), NbCompThermique)
      end do

      ! Sm(:,rowk) = BB * Sm(:,rowk), rowk
      Smk(:) = Sm(:, rowk)
      call dgemv('T', NbCompThermique, NbCompThermique, 1.d0, &
                 BB, NbCompThermique, Smk, 1, &
                 0.d0, Sm(:, rowk), 1)

   end subroutine JacobianMSWells_Alignment_man_row

   subroutine JacobianMSWells_free

      deallocate (SmQw)
      deallocate (Sm)

      deallocate (JacA%Pt)
      deallocate (JacA%Num)
      deallocate (JacA%Val)
#if  _DEBUG_JAC_MSWELLS_ == 1
      deallocate (JacANoAlig%Val)
#endif

      deallocate (csrSP)
      deallocate (csrS)

   end subroutine JacobianMSWells_free

   subroutine JacobianMSWells_print_IP_info_to_file(t) &
      bind(C, name="JacobianMSWells_print_IP_info_to_file")

      real(c_double), intent(in), value :: t
      integer :: s, k, nbwells
#if _DEBUG_JAC_IP_MSWELLS_ == 1
      if (commRank == 0) then !Master proc
         !Flowrate
         open (unit=110, file='qwwell.val', status='old', position='append', action='write')

         !Pressure at head
         open (unit=310, file='pwwell.val', status='old', position='append', action='write')

         !IP
         open (unit=410, file='QPwell.val', status='old', position='append', action='write')

         nbwells = NbMSWellLocal_Ncpus(commRank + 1)
         do k = 1, nbwells

            write (110, *) t, QwByMSWell(k)
            if (DataMSWellLocal(k)%IndWell == 'p') then
               write (410, *) t, 1
            else
               write (410, *) t, 2
            end if

            s = NodebyMSWellLocal%Pt(k + 1) !Head
            write (310, *) t, IncMSWell(s)%coats%Pression

         end do
         close (110)
         close (310)
         close (410)

      end if
#endif
   end subroutine JacobianMSWells_print_IP_info_to_file

   subroutine JacobianMSWells_print_LA_info_to_file(Delta_t, Iter) &
      bind(C, name="JacobianMSWells_print_LA_info_to_file")

      real(c_double), intent(in), value :: Delta_t
      integer(c_int), intent(in), value :: Iter
      integer :: s, k, num_s, unk_idx_s, m, kc
#if _DEBUG_JAC_MSWELLS_ == 1

      if (commRank == 0) then !Master proc
         !Residual
         open (unit=110, file='reswell_vec.val', status='old', position='append', action='write')

         !-Sm
         open (unit=210, file='sm_vec.val', status='old', position='append', action='write')

         write (110, *) '# dt,iternewtontotal:', Delta_t, Iter
         write (210, *) '# dt,iternewtontotal:', Delta_t, Iter
         do k = 1, NbMSWellLocal
            ! looping from  queue to  head
            do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
               unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown  idx for son node

               write (110, *) unk_idx_s, ResiduMSWell(1, unk_idx_s)
               write (110, *) unk_idx_s, ResiduMSWell(2, unk_idx_s)
               if (NbComp .eq. 2) then
                  write (110, *) unk_idx_s, ResiduMSWell(3, unk_idx_s)
               end if
               write (210, *) unk_idx_s, -Sm(1, unk_idx_s)
               write (210, *) unk_idx_s, -Sm(2, unk_idx_s)
               if (NbComp .eq. 2) then
                  write (210, *) unk_idx_s, -Sm(3, unk_idx_s)
               end if
            end do

         end do
         write (110, *) " "
         close (110)

         write (210, *) " "
         close (210)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Jacobian Not Aligned
         open (unit=310, file='jacobian.val', status='old', position='append', action='write')

         write (310, *) '# dt,iternewtontotal:', Delta_t, Iter
         do k = 1, NbMSWellLocal

            ! looping from  queue to  head
            do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
               unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown  idx for son node
               write (310, *) ' Row s ', unk_idx_s
               !!For all columns at the row of node s
               do m = JacA%Pt(unk_idx_s) + 1, JacA%Pt(unk_idx_s + 1)
                  kc = JacA%Num(m)

                  if (NbComp .eq. 1) then
                     !water2phase
                     write (310, *) ' JacMSWell 11 ', s, kc, JacANoAlig%Val(1, 1, m)
                     write (310, *) ' JacMSWell 12 ', s, kc, JacANoAlig%Val(1, 2, m)
                     write (310, *) ' JacMSWell 21 ', s, kc, JacANoAlig%Val(2, 1, m)
                     write (310, *) ' JacMSWell 22 ', s, kc, JacANoAlig%Val(2, 2, m)
                     write (310, *)
                  else if (NbComp .eq. 2) then
                     !immiscible
                     write (310, *) ' JacMSWell 11 ', s, kc, JacANoAlig%Val(1, 1, m)
                     write (310, *) ' JacMSWell 12 ', s, kc, JacANoAlig%Val(1, 2, m)
                     write (310, *) ' JacMSWell 13 ', s, kc, JacANoAlig%Val(1, 3, m)
                     write (310, *) ' JacMSWell 21 ', s, kc, JacANoAlig%Val(2, 1, m)
                     write (310, *) ' JacMSWell 22 ', s, kc, JacANoAlig%Val(2, 2, m)
                     write (310, *) ' JacMSWell 23 ', s, kc, JacANoAlig%Val(2, 3, m)
                     write (310, *) ' JacMSWell 31 ', s, kc, JacANoAlig%Val(3, 1, m)
                     write (310, *) ' JacMSWell 32 ', s, kc, JacANoAlig%Val(3, 2, m)
                     write (310, *) ' JacMSWell 33 ', s, kc, JacANoAlig%Val(3, 3, m)
                     write (310, *)
                  endif
               end do
            end do
         end do
         write (310, *) " "
         close (310)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Jacobian Aligned
         open (unit=310, file='jacobian_ali.val', status='old', position='append', action='write')

         write (310, *) '# dt,iternewtontotal:', Delta_t, Iter
         do k = 1, NbMSWellLocal

            ! looping from  queue to  head
            do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
               unk_idx_s = IncIdxNodebyMSWellLocal%Num(s)!unkown  idx for son node
               write (310, *) ' Row s ', unk_idx_s
               !!For all columns at the row of node s
               do m = JacA%Pt(unk_idx_s) + 1, JacA%Pt(unk_idx_s + 1)
                  kc = JacA%Num(m)
#if ComPASS_NUMBER_OF_COMPONENTS == 1
                  !water2phase
                  write (310, *) ' JacMSWell 11 ', s, kc, JacA%Val(1, 1, m)
                  write (310, *) ' JacMSWell 12 ', s, kc, JacA%Val(1, 2, m)
                  write (310, *) ' JacMSWell 21 ', s, kc, JacA%Val(2, 1, m)
                  write (310, *) ' JacMSWell 22 ', s, kc, JacA%Val(2, 2, m)
                  write (310, *)
#elif ComPASS_NUMBER_OF_COMPONENTS == 2
                  !immiscible
                  write (310, *) ' JacMSWell 11 ', s, kc, JacA%Val(1, 1, m)
                  write (310, *) ' JacMSWell 12 ', s, kc, JacA%Val(1, 2, m)
                  write (310, *) ' JacMSWell 13 ', s, kc, JacA%Val(1, 3, m)
                  write (310, *) ' JacMSWell 21 ', s, kc, JacA%Val(2, 1, m)
                  write (310, *) ' JacMSWell 22 ', s, kc, JacA%Val(2, 2, m)
                  write (310, *) ' JacMSWell 23 ', s, kc, JacA%Val(2, 3, m)
                  write (310, *) ' JacMSWell 31 ', s, kc, JacA%Val(3, 1, m)
                  write (310, *) ' JacMSWell 32 ', s, kc, JacA%Val(3, 2, m)
                  write (310, *) ' JacMSWell 33 ', s, kc, JacA%Val(3, 3, m)
                  write (310, *)
#endif
               end do
            end do
         end do
         write (310, *) " "
         close (310)

      end if
#endif
   end subroutine JacobianMSWells_print_LA_info_to_file

end module JacobianMSWells
