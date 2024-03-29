!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

#ifndef _THERMIQUE_
#error DefWellFlash: thermal transfer must be activated
#endif

!> \brief Define the flash to determine the phases
!! which are actualy present in each cell, and
!! the mode of the well (flowrate or pressure).
module DefFlashWells

#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use mpi
   use CommonMPI
   use DefModel
   use Physics
   use Thermodynamics
   use IncCVReservoir
   use MeshSchema
   use CommonType
   use CommonMPI
   use GlobalMesh
   use IncCVWells
   use IncCVReservoir
   use LoisThermoHydro
   use Thermodynamics
   use MeshSchema
   use Physics
   use DefModel
   use IncPrimSecd
   use WellState
#else
   use CommonType, only: CSRdble
   use CommonMPI, only: commRank, CommonMPI_abort
   use GlobalMesh, only: NodeRocktype
   use IncCVWells, only: &
      PerfoWellInj, &
      PerfoWellProd, &
      IncPressionWellProd, &
      IncPressionWellInj, &
      NodebyWellProdLocal, &
      NodebyWellInjLocal
   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, &
      IncNode
   use LoisThermoHydro, only: &
      DensitemolaireKrViscoCompNode, &
      DensitemolaireKrViscoEnthalpieNode, &
      DensitemolaireKrViscoCompWellInj, &
      LoisThermoHydro_divPrim_nodes, LoisThermoHydro_divP_wellinj
   use Thermodynamics, only: &
      ! f_MolarDensity,
      f_VolumetricMassDensity, f_MolarEnthalpy
   use MeshSchema, only: &
      DataWellInjLocal, &
      NodeDatabyWellInjLocal, &
      DataWellProdLocal, &
      NodeDatabyWellProdLocal, &
      XNodeLocal, &
      NbWellProdLocal_Ncpus, &
      NbWellInjLocal_Ncpus
   use Physics, only: gravity
   use DefModel, only: &
      IndThermique, NbPhase, NbComp, LIQUID_PHASE, MCP, &
      NumPhasePresente_ctx, NbPhasePresente_ctx
   use IncPrimSecd, only: IncPrimSecd_compPrim_nodes
   use WellState, only: &
      WellState_FlowrateWellProd, WellState_solve_for_temperature, &
      summolarFluxProd, sumnrjFluxProd
#endif

#ifndef ComPASS_SINGLE_PHASE
   use Thermodynamics, only: FluidThermodynamics_Tsat
   use DefModel, only: GAS_PHASE
#endif

   implicit none

   !integer, parameter, private:: WellsNslice = 100 !< Number of discretization for the approximation of the integral
   ! integer, private :: fdFl !< debug output

   type(CSRdble) :: ZSortedInj !< CSR vector storing the nodes of each injection well in Z coordinate order

   ! molar fluxes for injection well(s), head node
   double precision, allocatable, dimension(:) :: headmolarFluxInj !< Molar flux for injection well: headmolarFluxProd(num_well)

   public :: &
      DefFlashWells_allocate, & ! Allocation, initialization and deallocation (in NN.F90)
      DefFlashWells_NewtonFlashLinWells, & ! Flash after each Newton iteration
      DefFlashWells_TimeFlash_injectors, & ! Flash after each time step
      DefFlashWells_TimeFlash_producers_single_phase, &
      DefFlashWells_TimeFlash_producers_two_phases, &
      DefFlashWells_free

   private :: &
      DefFlashWells_NewtonFlashLinWellInj, & ! Flash after each Newton iteratio
      DefFlashWells_NewtonFlashLinWellProd, & ! Flash after each Newton iteration
      DefFlashWells_PressureToFlowrateWellProd, &
      DefFlashWells_PressureToFlowrateWellInj, &
      QuickSortCSR, &
      DefFlashWells_SortHeights_and_Init

contains

   subroutine DefFlashWells_NewtonFlashLinWells

      call DefFlashWells_NewtonFlashLinWellInj
      call DefFlashWells_NewtonFlashLinWellProd

   end subroutine DefFlashWells_NewtonFlashLinWells

   !> \brief Main subourtine, after each time iteration
   subroutine DefFlashWells_TimeFlash_injectors() &
      bind(C, name="DefFlashWells_TimeFlash_injectors")

      integer :: num_Well, s
      double precision :: Tw

      ! compute
      do num_Well = 1, NbWellInjLocal_Ncpus(commRank + 1)
         call DefFlashWells_PressureToFlowrateWellInj(num_Well, headmolarFluxInj(num_Well))
         ! print*, "head ", headmolarFluxInj(num_Well)
         ! update injection temperature
         ! taken to be the same as well head temperature for the whole well
         ! in the absence of thermal diffusion
         ! FIXME: the effect of pressure is not taken into account
         Tw = DataWellInjLocal(num_Well)%InjectionTemperature
         do s = NodebyWellInjLocal%Pt(num_Well) + 1, NodebyWellInjLocal%Pt(num_Well + 1)
            PerfoWellInj(s)%Temperature = Tw
         end do
      end do

   end subroutine DefFlashWells_TimeFlash_injectors

   !> Allocate global vectors used only in this file
   subroutine DefFlashWells_allocate

      ! character(len=300) :: fn, pid

!       ! put debugging magic here, open a file descriptor per processor
!       ! for more readability
! #define _DEBUG_LVL1_
! #if defined _DEBUG_ && defined _DEBUG_LVL1_
!       if (.false.) then
!          fdFl = 6 ! stdout: 6, scratch (temprary discarded file): 12
!       else
!          fdFl = 12
!          write (pid, *) getpid()
!          write (fn, '(a,a,a)') 'proc.', trim(adjustl(pid)), '.log'
!          open (unit=fdFl, file=trim(fn), action='WRITE')
!       end if
! #else
!       fdFl = 12
!       open (unit=fdFl, status='SCRATCH')
! #endif

      ! allocate flowrate
      allocate (headmolarFluxInj(NbWellInjLocal_Ncpus(commRank + 1)))

      ZSortedInj%Nb = NodebyWellInjLocal%Nb
      allocate (ZSortedInj%Pt(ZSortedInj%Nb + 1))
      ZSortedInj%Pt = NodebyWellInjLocal%Pt
      allocate (ZSortedInj%Num(ZSortedInj%Pt(ZSortedInj%Nb + 1)))
      allocate (ZSortedInj%Val(ZSortedInj%Pt(ZSortedInj%Nb + 1)))

      ! init sort height
      call DefFlashWells_SortHeights_and_Init

   end subroutine DefFlashWells_allocate

   subroutine DefFlashWells_free()

      ! call CommonType_deallocCSRdble(ZSortedInj)
      ! close(unit=12)
      ! call CommonType_deallocCSRdble(ZSortedInj)
      deallocate (headmolarFluxInj)

   end subroutine DefFlashWells_free

   !> \brief Performs some initialization, and sorts the injector well nodes
   subroutine DefFlashWells_SortHeights_and_Init
      integer :: k, s

      ! init these module global variables
      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)
         do s = ZSortedInj%Pt(k) + 1, ZSortedInj%Pt(k + 1)
            ZSortedInj%Num(s) = s
            ZSortedInj%Val(s) = XNodeLocal(3, NodebyWellInjLocal%Num(s))
         end do
      end do

      !! sorting the well nodes according to the heights
      do k = 1, NbWellInjLocal_Ncpus(commRank + 1)
         call QuickSortCSR(ZSortedInj, ZSortedInj%Pt(k) + 1, ZSortedInj%Pt(k + 1), 'd')
      end do

      ! write (fdFl, *) '{'
      ! write (fdFl, *) ZSortedInj%Num
      ! write (fdFl, *) ZSortedInj%Val
      ! write (fdFl, *) '}{'
      ! do s = 1, ZSortedInj%Pt(ZSortedInj%Nb + 1)
      !    write (fdFl, *) 'ptL', ZSortedInj%Num(s), 'numL', NodebyWellInjLocal%Num(ZSortedInj%Num(s)), &
      !       ZSortedInj%Val(s)
      ! end do
      ! write (fdFl, *) '}'

   end subroutine DefFlashWells_SortHeights_and_Init

   !> \brief Determine the mode of the injection well
   !! (flowrate or pressure).
   !!
   !! As long as the pressure is less or egal to
   !! the pressure max, the flowrate of the well
   !! is imposed. If the pressure is too high,
   !! the flowrate is no more fixed and the pressure
   !! is set as Pressure max.
   subroutine DefFlashWells_NewtonFlashLinWellInj

      double precision :: Flowrate_head
      integer :: num_Well

      ! print*, NbWellInjLocal_Ncpus

      do num_Well = 1, NbWellInjLocal_Ncpus(commRank + 1)

         if (DataWellInjLocal(num_Well)%IndWell == 'c') cycle ! well is closed

         if (DataWellInjLocal(num_Well)%IndWell == 'f') then ! flowrate mode

            if (IncPressionWellInj(num_Well) > DataWellInjLocal(num_Well)%PressionMax) then

               DataWellInjLocal(num_Well)%IndWell = 'p' ! change to pressure mode
               IncPressionWellInj(num_Well) = DataWellInjLocal(num_Well)%PressionMax ! Pw = PwMax

#ifdef COMPASS_LOG_WELL_INFO
               write (*, *) '[Well-Monitoring] Injector has changed to pressure mode, well number: ', num_Well
#endif

            endif

         else if (DataWellInjLocal(num_Well)%IndWell == 'p') then ! pressure mode

            if (IncPressionWellInj(num_Well) > DataWellInjLocal(num_Well)%PressionMax) then
               IncPressionWellInj(num_Well) = DataWellInjLocal(num_Well)%PressionMax ! With Newton inc, Pw>Pmax, then change it
            endif

            ! compute the new flowrate at the head node of well num_Well
            call LoisThermoHydro_divP_wellinj(num_Well) ! update thermo Laws of nodes in well num_Well
            call DefFlashWells_PressureToFlowrateWellInj(num_Well, Flowrate_head)

            if (abs(Flowrate_head) > abs(DataWellInjLocal(num_Well)%ImposedFlowrate)) then ! inj well then DataWellInjLocal(num_Well)%flowrate < 0
               DataWellInjLocal(num_Well)%IndWell = 'f' ! change to flowrate mode

#ifdef COMPASS_LOG_WELL_INFO
               write (*, *) '[Well-Monitoring] Injector has changed to flowrate mode, well number: ', num_Well
#endif

            endif
         else

            print *, "Error in Newton Flash Injection Well", num_Well, "  no such index well: ", DataWellInjLocal(num_Well)%IndWell
         end if

      end do ! well

   end subroutine DefFlashWells_NewtonFlashLinWellInj

   !> \brief Determine the mode of the production well
   !! (flowrate or pressure).
   !!
   !! As long as the pressure is greater or egal to
   !! the pressure min, the flowrate of the well
   !! is imposed. If the pressure is too low,
   !! the flowrate is no more fixed and the pressure
   !! is set as Pressure min.
   subroutine DefFlashWells_NewtonFlashLinWellProd

      double precision :: Flowrate_head
      integer :: num_Well

      do num_Well = 1, NbWellProdLocal_Ncpus(commRank + 1)

         if (DataWellProdLocal(num_Well)%IndWell == 'c') cycle ! well is closed

         if (DataWellProdLocal(num_Well)%IndWell == 'f') then ! flowrate mode

            if (IncPressionWellProd(num_Well) < DataWellProdLocal(num_Well)%PressionMin) then

               DataWellProdLocal(num_Well)%IndWell = 'p' ! change to pressure mode
               IncPressionWellProd(num_Well) = DataWellProdLocal(num_Well)%PressionMin ! Pw = PwMin

#ifdef COMPASS_LOG_WELL_INFO
               write (*, *) '[Well-Monitoring] Producer has changed to pressure mode, well number: ', num_Well
#endif

            endif

         else if (DataWellProdLocal(num_Well)%IndWell == 'p') then ! pressure mode

            if (IncPressionWellProd(num_Well) < DataWellProdLocal(num_Well)%PressionMin) then
               IncPressionWellProd(num_Well) = DataWellProdLocal(num_Well)%PressionMin ! With Newton inc, Pw<Pmin, then change it
            endif

            ! compute the new flowrate at the head node of the well num_Well
            call DefFlashWells_PressureToFlowrateWellProd(num_Well, Flowrate_head)

            if (abs(Flowrate_head) > abs(DataWellProdLocal(num_Well)%ImposedFlowrate)) then ! Prod well then DataWellProdLocal(num_Well)%flowrate > 0
               DataWellProdLocal(num_Well)%IndWell = 'f' ! change to flowrate mode

#ifdef COMPASS_LOG_WELL_INFO
               write (*, *) '[Well-Monitoring] Producer has changed to flowrate mode, well number: ', num_Well
#endif

            endif
         else

            print *, "Error in Newton Flash Production Well", num_Well, &
               "  no such index well: ", DataWellProdLocal(num_Well)%IndWell
         end if

      enddo ! well

   end subroutine DefFlashWells_NewtonFlashLinWellProd

   !> \brief Compute the flowrate at the head node of the production well num_Well
   !! with the new value of the unknows
   !!
   !! Use IncPressionWellProd, PerfoWellProd\%PressureDrop and IncNode\%Pression.
   !!
   !! \param[in]      num_Well             Numero of the production well
   !! \param[out]     Qw                   Flowrate at the head node
   subroutine DefFlashWells_PressureToFlowrateWellProd(num_Well, Qw)

      integer, intent(in) :: num_Well
      double precision, intent(out) :: Qw

      integer :: s, nums, icp, m, mph
      double precision :: Pws, Ps, WIDws
      double precision:: Flux_ks(NbComp)

      Qw = 0.d0

      ! update prim/secd arrays of all nodes
      call IncPrimSecd_compPrim_nodes
      ! update thermo Laws of all nodes
      call LoisThermoHydro_divPrim_nodes

      ! nodes of well num_Well
      do s = NodebyWellProdLocal%Pt(num_Well) + 1, NodebyWellProdLocal%Pt(num_Well + 1)
         nums = NodebyWellProdLocal%Num(s)

         Pws = IncPressionWellProd(num_Well) + PerfoWellProd(s)%PressureDrop ! P_{w,s}
         Ps = IncNode(nums)%Pression ! P_s
         WIDws = NodeDatabyWellProdLocal%Val(s)%WID

         Flux_ks(:) = 0.d0

         do m = 1, NbPhasePresente_ctx(IncNode(nums)%ic) ! Q_s
            mph = NumPhasePresente_ctx(m, IncNode(nums)%ic)

            !
            if ((Ps - Pws) > 0.d0) then
               do icp = 1, NbComp
                  if (MCP(icp, mph) == 1) then ! \cap P_i
                     Flux_ks(icp) = Flux_ks(icp) + DensiteMolaireKrViscoCompNode(icp, mph, nums)*WIDws*(Ps - Pws)
                  end if
               end do
            end if
         end do

         do icp = 1, NbComp
            Qw = Qw + Flux_ks(icp)
         end do

      end do

   end subroutine DefFlashWells_PressureToFlowrateWellProd

   !> \brief Compute the flowrate of the injection well num_Well.
   !! That is the total flowrate that goes through the injector well to the reservoir.
   !! Use PerfoWellInj, IncPressionWellInj, WID and IncNode\%Pression.
   !!
   !! \param[in]      num_Well             Numero of the injection well
   !! \param[out]     Qw    Flowrate at the head node of the injection well
   subroutine DefFlashWells_PressureToFlowrateWellInj(num_Well, Qw)

      integer, intent(in) :: num_Well
      double precision, intent(inout) :: Qw

      integer :: s, nums, icp
      double precision :: Pws, Ps, deltaPs, qws
      double precision:: Flux_ks(NbComp)

      Qw = 0.d0
      ! nodes of well num_Well
      do s = NodebyWellInjLocal%Pt(num_Well) + 1, NodebyWellInjLocal%Pt(num_Well + 1)
         nums = NodebyWellInjLocal%Num(s)
         Pws = IncPressionWellInj(num_Well) + PerfoWellInj(s)%PressureDrop
         Ps = IncNode(nums)%Pression
         deltaPs = Ps - Pws
         if (deltaPs < 0.d0) then
! #ifndef NDEBUG
!             if (NbPhasePresente_ctx(IncNode(nums)%ic) /= 1) &
!                call CommonMPI_abort("Injectors are supposed to be monophasic.")
! #endif
            qws = NodeDatabyWellInjLocal%Val(s)%WID*deltaPs
            do icp = 1, NbComp
               Flux_ks(icp) = DensiteMolaireKrViscoCompWellInj(icp, s)*qws
            end do
            Qw = Qw + sum(Flux_ks)
         end if
      end do

   end subroutine DefFlashWells_PressureToFlowrateWellInj

   !> \brief Sorting the heights contained in mycsr%Value, and update
   !! the corresponding node values (indexes) stored in mycsr%Num
   !! mycsr%Pt is constructed on NodebyWellInjLocal
   recursive subroutine QuickSortCSR(myCSR, left, right, mode)
      use commontype, only: CSRdble
      implicit none
      type(CSRdble), intent(inout) :: myCSR
      integer, intent(in) :: left, right
      character(len=1), intent(in) :: mode
      double precision :: x, tmp_val
      integer :: i, j, tmp_num

      ! write(fdFl, *) '--> quicksortCSR called with left: ', left, ' right: ', right
      ! we are sorting CSR%Num(f), f=CSR%Pt(k)+1,CSR%Pt(k+1), for well k

      x = myCSR%Val((left + right)/2) ! pivot point
      ! write(fdFl, *) 'pivot is ', x
      i = left; j = right ! left and right boundaries
      do
         if (mode == 'd') then
            do while (x < myCSR%Val(i))
               i = i + 1
            end do
            do while (x > myCSR%Val(j))
               j = j - 1
            end do
         elseif (mode == 'i') then
            do while (x > myCSR%Val(i))
               i = i + 1
            end do
            do while (x < myCSR%Val(j))
               j = j - 1
            end do
         else
            write (*, *) 'WARNING: sorting mode unknown !'
            return
         end if
         if (i >= j) then
            exit
         end if

         ! swapping the %Val field and %Num fields
         tmp_val = myCSR%Val(i); tmp_num = myCSR%Num(i)
         myCSR%Val(i) = myCSR%Val(j); myCSR%Num(i) = myCSR%Num(j)
         myCSR%Val(j) = tmp_val; myCSR%Num(j) = tmp_num
         i = i + 1; j = j - 1
      end do
      if (left < i - 1) then
         call QuickSortCSR(myCSR, left, i - 1, mode)
      end if
      if (j + 1 < right) then
         call QuickSortCSR(myCSR, j + 1, right, mode)
      end if
   end subroutine QuickSortCSR

   !> \brief Update PressureDrop and Pression of all the nodes of a given well
   ! we start by the head node and update the new value using rho*gravity*delta_z increments
   ! LIMITATIONS: the head node of a given well must be at the top height
   ! pressure update loop:  p^{n+1} = p^{n} + rho(p^{n}) * g * (z^{n+1} - z^{n})
   ! (z axis)
   ! ^
   ! | ----- (PerfoWellInj(pts1)%Pression)
   ! |    ^
   ! |    .
   ! |    .
   ! |    delta(P)
   ! |    .
   ! |    .
   ! |    .
   ! |    v
   ! | ----- (PerfoWellInj(pts2)%Pression)

   ! subroutine DefFlashWells_PressureDropInj

   !   double precision :: Sat(NbPhase)
   !   double precision :: Rhotmp, dz, Ptmp, ztmp
   !   integer :: i, s, k, nums1, nums2, pts1, pts2

   !   do k=1, NbWellInjLocal_Ncpus(commRank+1)

   !      ! initialize the head pressure drop to zero, and the head pressure
   !     ! with the Unknown value for the Well
   !     PerfoWellInj(ZSortedInj%Num(ZSortedInj%Pt(k)+1))%PressureDrop = 0.d0
   !     PerfoWellInj(ZSortedInj%Num(ZSortedInj%Pt(k)+1))%Pression = IncPressionWellInj(k)

   !     do s=ZSortedInj%Pt(k)+1, ZSortedInj%Pt(k+1) - 1
   !       write(fdFl, *) '.. .. .. .. .. .. .. .. .. .. .. .. .. .. ..'
   !       ! node pointer
   !       pts1 = ZsortedInj%Num(s)
   !       pts2 = ZsortedInj%Num(s+1)
   !       ! local node number
   !       nums1 = NodebyWellInjLocal%Num(pts1)
   !       nums2 = NodebyWellInjLocal%Num(pts2)

   !       Sat(GAS_PHASE) = 0.d0
   !       Sat(LIQUID_PHASE) = 1.d0

   !       Ptmp = PerfoWellInj(pts1)%Pression
   !       ztmp = XNodeLocal(3, nums1) ! just to check we are OK

   !       dz = (XNodeLocal(3, nums2) - XNodeLocal(3, nums1)) / dble(WellsNslice)
   !       ! pressure update loop:  p^{n+1} = p^{n} + rho(p^{n}) * g * (z^{n+1} - z^{n})
   !       do i=1, WellsNslice
   !         ztmp = ztmp + dz
   !         Rhotmp = f_MolarDensity(LIQUID_PHASE, Ptmp, PerfoWellInj(pts1)%Temperature, &
   !             DataWellInjLocal(k)%CompTotal, Sat)
   !         ! gravity points downwards, heights points upwards, hence the negative sign
   !         Ptmp = Ptmp - gravity * Rhotmp * dz
   !       end do

   !       PerfoWellInj(pts2)%Pression = Ptmp
   !       PerfoWellInj(pts2)%PressureDrop = PerfoWellInj(pts2)%Pression - IncPressionWellInj(k)

   !       ! debugging
   !       write(fdFl, '(4(a,i2),a)') 'nodes (', pts1, 'ptL,', nums1, &
   !           'numL) -> (', pts2, 'ptL,', nums2, 'numL)'

   !       write(fdFl, '(a,f6.2,a,f6.2,f6.2,a,e14.6,a,e14.6)'), &
   !           'integrate z1: ', XNodeLocal(3, nums1), ' -> z2: ', XNodeLocal(3, nums2), &
   !           ztmp, ' P1: ', PerfoWellInj(pts1)%Pression, ' P2: ', Ptmp
   !       write(fdFl, '(a,e14.6)') 'pressure drop', PerfoWellInj(pts2)%PressureDrop
   !     end do
   !   end do
   ! end subroutine DefFlashWells_PressureDropInj

   !! \brief Execute the flash to determine T and the molar fractions
   !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
   !! This Flash is performed for a monophasic multicomponent fluid.
   !!
   !! Loop over the nodes s from head to tail to
   !! to determine which phases are present, then computes the temperature
   !! and the mean density.
   !! The well pressure and the pressure drop are taken from
   !! the previous Newton iteration and the previous time step respectively.
   subroutine DefFlashWells_TimeFlash_producers_single_phase() &
      bind(C, name="DefFlashWells_TimeFlash_producers_single_phase")

      double precision :: T, ResT, Pws, Ci(NbComp), sumni, E
      integer :: k, s
      logical :: converged

      call WellState_FlowrateWellProd

      do k = 1, NodebyWellProdLocal%Nb

         ! looping from head to queue
         do s = NodebyWellProdLocal%Pt(k + 1), NodebyWellProdLocal%Pt(k) + 1, -1

            Pws = PerfoWellProd(s)%Pression

            ! Newton method to compute T: R = E-Enthalpie*n = 0
            E = sumnrjFluxProd(s) ! energy

            sumni = sum(summolarFluxProd(:, s))

            ! initialize newton with reservoir temperature
            T = IncNode(NodebyWellProdLocal%Num(s))%Temperature
            call WellState_solve_for_temperature(LIQUID_PHASE, E, Pws, T, sumni, converged, ResT)
            if (.not. converged) then
               print *, "Warning: Newton in DefFlashWells_TimeFlash_producers_single_phase has not converged"
               print *, "Residue is", abs(ResT), "Well_idx= ", k, "node_idx= ", s
            end if

            ! update PhysPerfoWell
            PerfoWellProd(s)%Temperature = T
            PerfoWellProd(s)%Density = f_VolumetricMassDensity(LIQUID_PHASE, Pws, T, Ci)

         end do ! end of well k
      end do

   end subroutine DefFlashWells_TimeFlash_producers_single_phase

#ifndef ComPASS_SINGLE_PHASE
   subroutine DefFlashWells_copy_reservoir_states(wk)
      integer, intent(in) :: wk

      double precision :: Tws, Pws, Sg, Sl
      double precision :: rhog, rhol, rho
      integer :: s, sr

      do s = NodebyWellProdLocal%Pt(wk) + 1, NodebyWellProdLocal%Pt(wk)
         sr = NodebyWellProdLocal%Num(s) ! reservoir node index
         Pws = IncNode(sr)%Pression
         Tws = IncNode(sr)%Temperature
         Sg = IncNode(sr)%Saturation(GAS_PHASE)
         Sl = IncNode(sr)%Saturation(LIQUID_PHASE)
         PerfoWellProd(s)%Pression = Pws
         PerfoWellProd(s)%Temperature = Tws
         PerfoWellProd(s)%Saturation = IncNode(sr)%Saturation
         rho = 0.d0
         if (Sg > 0) then
            rhog = f_VolumetricMassDensity(GAS_PHASE, Pws, Tws, IncNode(sr)%Comp(:, GAS_PHASE))
            rho = rho + Sg*rhog
         end if
         if (Sl > 0) then
            rhol = f_VolumetricMassDensity(LIQUID_PHASE, Pws, Tws, IncNode(sr)%Comp(:, LIQUID_PHASE))
            rho = rho + Sl*rhol
         end if
         PerfoWellProd(s)%Density = rho
      enddo ! loop on perforations s of well wk

   end subroutine DefFlashWells_copy_reservoir_states
#endif

   !! \brief Execute the flash to determine T and the molar fractions
   !! to update PerfoWellProd(s)%Temperature and PerfoWellProd(s)%Density.
   !! This Flash is performed for a diphasique monocomponent fluid (H20).
   !!
   !! Loop over the nodes s from head to tail to
   !! to determine which phases are present, then computes the temperature
   !! and the mean density.
   !! The well pressure and the pressure drop are taken from
   !! the previous Newton iteration and the previous time step respectively.
   subroutine DefFlashWells_TimeFlash_producers_two_phases() &
      bind(C, name="DefFlashWells_TimeFlash_producers_two_phases")

#ifdef ComPASS_SINGLE_PHASE

      call CommonMPI_abort("You cannot have a two phase well model with a single phase physics!")

#else

      double precision :: Tws, Pws, xg, Sg, Sl, ResT
      double precision :: Hgas, Hliq, rhog, rhol
      double precision :: sumni, E
      double precision :: dP_Tsat, Cw(NbComp)
      integer :: wk, s, ID_PHASE ! ID_PHASE=(-1 if diphasique, GAS_PHASE if gas, LIQUID_PHASE if liq)
      logical :: converged

      ! compute flowrate of well wk (fill summolarFluxProd and sumnrjFluxProd)
      call WellState_FlowrateWellProd

      ! loop over production well
      do wk = 1, NodebyWellProdLocal%Nb

         if (DataWellProdLocal(wk)%IndWell == 'c') then ! well is closed

            call DefFlashWells_copy_reservoir_states(wk)

         else ! well is open

            ! looping from head to queue
            do s = NodebyWellProdLocal%Pt(wk + 1), NodebyWellProdLocal%Pt(wk) + 1, -1

               Pws = PerfoWellProd(s)%Pression
               E = sumnrjFluxProd(s) ! energy
               sumni = sum(summolarFluxProd(:, s))
               if (sumni < 1.0D-12) then
                  cycle ! FIXME: keep everything as the previous timestep
               end if
               Cw = summolarFluxProd(:, s)/sumni
               ! suppose that the two phases are present at the node
               ID_PHASE = -1
               ! set temperature to saturation temperature at Pws
               call FluidThermodynamics_Tsat(Pws, Tws, dP_Tsat)
               ! thus compute liq_molarfrac thanks to the energy, and the enthalpies
               ! molarFrac is not used in the computation of the enthalpies
               Hgas = f_MolarEnthalpy(GAS_PHASE, Pws, Tws, Cw)
               Hliq = f_MolarEnthalpy(LIQUID_PHASE, Pws, Tws, Cw)

               xg = (E/sumni - Hliq)/(Hgas - Hliq)

               if (xg <= 0.d0) then ! the hypothesis that the two phases are present is wrong: only liquid
                  xg = 0.d0
                  ID_PHASE = LIQUID_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 0.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 1.d0
               else if (xg >= 1.d0) then ! the hypothesis that the two phases are present is wrong: only gas
                  xg = 1.d0
                  ID_PHASE = GAS_PHASE
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = 1.d0
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = 0.d0
               endif

               if (ID_PHASE > 0) then ! perforation is monophasic
                  ! compute temperature: initialize newton with reservoir temperature
                  Tws = IncNode(NodebyWellProdLocal%Num(s))%Temperature
                  call WellState_solve_for_temperature(ID_PHASE, E, Pws, Tws, sumni, converged, ResT)
                  if (.not. converged) then
                     print *, "Warning: Newton in DefFlashWells_TimeFlash_producers_two_phases has not converged"
                     print *, "Residue is", abs(ResT), "Well_idx= ", wk, "node_idx= ", s
                  end if
                  PerfoWellProd(s)%Temperature = Tws
                  PerfoWellProd(s)%Density = f_VolumetricMassDensity(ID_PHASE, Pws, Tws, Cw)
               else ! perforation is diphasic
                  PerfoWellProd(s)%Temperature = Tws
                  ! molarFrac is not used in the computation of the massique densities
                  rhog = f_VolumetricMassDensity(GAS_PHASE, Pws, Tws, Cw)
                  rhol = f_VolumetricMassDensity(LIQUID_PHASE, Pws, Tws, Cw)
                  Sg = (xg/rhog)/((xg/rhog) + ((1.d0 - xg)/rhol))
                  Sl = 1.d0 - Sg
                  PerfoWellProd(s)%Density = Sl*rhol + Sg*rhog
                  PerfoWellProd(s)%Saturation(GAS_PHASE) = Sg
                  PerfoWellProd(s)%Saturation(LIQUID_PHASE) = Sl
               end if ! perforation is monophasic / diphasic

            enddo ! loop on perforations s of well wk

         end if ! well is closed / open

      enddo ! loop on wells wk

#endif

   end subroutine DefFlashWells_TimeFlash_producers_two_phases

end module DefFlashWells
