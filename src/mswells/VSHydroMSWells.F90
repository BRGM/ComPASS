!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!INFO:
!     Model for mswells: Phase superficial velocities.
!     See
!        i) "Multi-segmented non-isothermal compositional liquid gas wellmodel" of Castanon Quiroz & Masson
!
!*This module has been adapted from R. Masson's mswells code.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module VSHydroMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding, only: c_double
   use CommonMPI
   use DefModel
   use NumbyContext
   use IncCVReservoir
   use IncCVMSWells
   use LoisThermoHydroMSWells
   use MeshSchema
   use LoisThermoHydro
   use Physics
   use DFMHydroMSWells
   use MSWellsData
#else

   use iso_c_binding, only: c_double
   use CommonMPI, only: commRank, CommonMPI_abort

   use DefModel, only: &
      NbPhase, NbComp, IndThermique, LIQUID_PHASE, MCP, &
      NbEqEquilibreMax, NbIncPTCMax, NbIncTotalPrimMax, &
      NbIncTotalMax, NbEqFermetureMax, &
      NbPhasePresente_ctx, NumPhasePresente_ctx

#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: &
      GAS_CONTEXT, LIQUID_CONTEXT, DIPHASIC_CONTEXT, GAS_PHASE
#endif

   use NumbyContext, only: &
      NumCompEqEquilibre_ctx, Num2PhasesEqEquilibre_ctx, NumIncComp2NumIncPTC_ctx, &
      NbEqEquilibre_ctx, NbEqFermeture_ctx, NumCompCtilde_ctx, &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
      NbCompCtilde_ctx, NbIncPTC_ctx, NbIncTotalPrim_ctx

   use IncCVReservoir, only: &
      TYPE_IncCVReservoir, NumPhasePresente_ctx, NbPhasePresente_ctx

   use IncCVMSWells, only: &
      IncMSWell

   use MeshSchema, only: &
      XNodeLocal, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus

   use LoisThermoHydro, only: &
      ContextInfo, &
      LoisThermoHydro_init_cv

   use LoisThermoHydroMSWells, only: &
      LoisThermoHydroMSWell, &
      TYPE_LoisThermoHydro

   use Physics, only: gravity
   use DFMHydroMSWells, only: &
      DFMHydroMSWells_MixtureVelocity, DFMHydroMSWells_FluxVsg, &
      DFMHydroMSWells_f_C0, DFMHydroMSWells_df_C0, &
      DFMHydroMSWells_df_sgC0, DFMHydroMSWells_f_sgC0K, DFMHydroMSWells_df_sgC0K, &
      DFMHydroMSWells_f_G, DFMHydroMSWells_df_G, DFMHydroMSWells_DPhi

   use MSWellsData, only: MSWellDataNodebyMSWell

#endif

   implicit none

   double precision, parameter ::  sigmagl = 71.97d-3 !Gas-Liquid interfacial tension
   !Special data to compute the superficial velocities stuff
   !For all input/output arguments withe postfix _ss
   ! indicates that the 1-index of the last entry is w.r s_prime and the 2-index is w.r the node s (the parent of s_prime)
   !Moreover, take account that these derivatives are taking with respect the primary unkowns of each node (which can be different)
   type :: TYPE_VSHydroMSWEdge

      real(c_double) :: MixVel ! mixture vel
      real(c_double) :: Svel(NbPhase), & !   superficial velocity
                        divSvel_ss(NbIncTotalPrimMax, NbPhase, 2), &
                        SmSvel(NbPhase)

   end type TYPE_VSHydroMSWEdge

   type :: TYPE_VSHydroMSWells

      type(TYPE_VSHydroMSWEdge) :: edata ! node son has the edge data
      ! the well head has the bdc data
   end type TYPE_VSHydroMSWells

   type :: TYPE_VSHydroMSWellsHead ! Extra well head data

      real(c_double) :: divQwVsHead(NbPhase)

   end type TYPE_VSHydroMSWellsHead

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   type(TYPE_VSHydroMSWells), allocatable, dimension(:), target, public :: &
      VSHydroMSWell   !<  Injector-Producer-MSWells-VSHydro

   ! Info at the head
   type(TYPE_VSHydroMSWellsHead), allocatable, dimension(:), target, protected :: &
      VSHydroHeadMSWell   !<  Injector/Producer-MSWells-VSHydro-HeadData

   public :: &
      VSHydroMSWells_allocate, &
      VSHydroMSWells_free, &
      VSHydroMSWells_compute, &
      VSHydroMSWells_init, &
      VSHydroMSWells_compute_bdcs_prod_well
   private :: &
      VSHydroMSWells_compute_bdc_prod_cv, &
      VSHydroMSWells_compute_cv, &
      VSHydroMSWells_init_cv, &
      VSHydroMSWells_ammix_densitemassique_cv, &
      VSHydroMSWells_ammix_viscosite_cv, &
      VSHydroMSWells_amean_densitemassique_cv

contains

   subroutine VSHydroMSWells_allocate()

      integer :: Nb, Nnz

      Nb = NodebyMSWellLocal%Nb
      Nnz = NodebyMSWellLocal%Pt(Nb + 1)
      allocate (VSHydroMSWell(Nnz))

      allocate (VSHydroHeadMSWell(NbMSWellLocal_Ncpus(commRank + 1)))

   end subroutine VSHydroMSWells_allocate

   subroutine VSHydroMSWells_free()

      deallocate (VSHydroMSWell)
      deallocate (VSHydroHeadMSWell)

   end subroutine VSHydroMSWells_free

   !Computes initial pressure
   subroutine VSHydroMSWells_init() &
      bind(C, name="VSHydroMSWells_init")

      integer :: s, sparent, k, nbwells, scoords_idx_s, scoords_idx_sp

      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do k = 1, nbwells

         !Only for producers
         if (DataMSWellLocal(k)%Well_type == 'i') cycle

         ! looping from  head to queue
         do s = NodebyMSWellLocal%Pt(k + 1), NodebyMSWellLocal%Pt(k) + 1, -1

            !set default vals to zero execpt MixVel
            VSHydroMSWell(s)%edata%Svel(:) = 0.d0
            VSHydroMSWell(s)%edata%divSvel_ss(:, :, :) = 0.d0
            VSHydroMSWell(s)%edata%SmSvel(:) = 0.d0

            if (s == NodebyMSWellLocal%Pt(k + 1)) cycle !head node

            sparent = NodeDatabyMSWellLocal%Val(s)%PtParent ! parent idx
            scoords_idx_s = NodebyMSWellLocal%Num(s)    ! index of the space coords of node s
            scoords_idx_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the space coords of parent node

            call VSHydroMSWells_init_cv( &
               IncMSWell(s)%coats, &
               IncMSWell(sparent)%coats, &
               s, &
               sparent, &
               scoords_idx_s, &
               scoords_idx_sp, &
               MSWellDataNodebyMSWell(s)%edata%cross_section, &
               MSWellDataNodebyMSWell(s)%edata%orientation, &
               VSHydroMSWell(s))  !recall: the son node has the data for each edge.

         end do
      end do

   end subroutine VSHydroMSWells_init

   subroutine VSHydroMSWells_compute() &
      bind(C, name="VSHydroMSWells_compute")

      integer :: s, sparent, k, nbwells, scoords_idx_s, scoords_idx_sp

      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do k = 1, nbwells

         !Only for producers
         if (DataMSWellLocal(k)%Well_type == 'i') cycle

         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            if (s == NodebyMSWellLocal%Pt(k + 1)) cycle !head node

            sparent = NodeDatabyMSWellLocal%Val(s)%PtParent ! parent idx
            scoords_idx_s = NodebyMSWellLocal%Num(s)    ! index of the space coords of node s
            scoords_idx_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the space coords of parent node

            call VSHydroMSWells_compute_cv( &
               IncMSWell(s)%coats, &
               IncMSWell(sparent)%coats, &
               s, &
               sparent, &
               scoords_idx_s, &
               scoords_idx_sp, &
               MSWellDataNodebyMSWell(s)%edata%cross_section, &
               MSWellDataNodebyMSWell(s)%edata%orientation, &
               VSHydroMSWell(s) &  !recall: the son node has the data for each edge.
               )

         end do
      end do

   end subroutine VSHydroMSWells_compute

   !Compute BDCs for a producer well which has mswell_idx
   subroutine VSHydroMSWells_compute_bdcs_prod_well(mswell_idx, qw)
      integer, intent(in)                     :: mswell_idx
      double precision, intent(in)           :: qw
      integer :: s, scoords_idx_s

      ! For only the  head
      s = NodebyMSWellLocal%Pt(mswell_idx + 1)

      scoords_idx_s = NodebyMSWellLocal%Num(s)    ! index of the space coords of node s

      call VSHydroMSWells_compute_bdc_prod_cv( &
         qw, &
         IncMSWell(s)%coats, &
         s, &
         scoords_idx_s, &
         MSWellDataNodebyMSWell(s)%edata%cross_section, &
         VSHydroMSWell(s), &  !recall: the head node has the bdc data
         VSHydroHeadMSWell(mswell_idx) &
         )

   end subroutine VSHydroMSWells_compute_bdcs_prod_well

   subroutine VSHydroMSWells_compute_bdc_prod_cv( &
      qw, &
      inc_s, &
      well_idx_s, &
      scoords_idx_s, &
      cross_section, &
      vshydro_s, &
      vshydro_h &
      )
      double precision, intent(in)           :: qw
      type(TYPE_IncCVReservoir), intent(in)  :: inc_s
      integer, intent(in)                    :: well_idx_s, scoords_idx_s
      double precision, intent(in)           :: cross_section
      type(TYPE_VSHydroMSWells), intent(inout) :: vshydro_s
      type(TYPE_VSHydroMSWellsHead), intent(inout) :: vshydro_h

#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort("Multi-segmented wells are only implemented for diphasic physics!")

#else

      !tmp
      double precision ::ratiod, zeta_s, Um, SmUm, divUm(NbIncTotalPrimMax), ss, Uc
      double precision :: dRhoUc(NbPhase), divUc(NbIncTotalPrimMax), SmUc
      double precision :: divratiod(NbIncTotalPrimMax), Smratiod, sC0, divsC0(NbIncTotalPrimMax), sC0K, divsC0K(NbIncTotalPrimMax)
      double precision ::fG, dsG, dratioG, sUd, divG(NbIncTotalPrimMax), SmG, divsUd(NbIncTotalPrimMax), SmsUd, Vsg, Vsl
      double precision Zeta(NbPhase), gas_sat, divZeta(NbIncTotalPrimMax, NbPhase), SmZeta(NbPhase)
      double precision Rho(NbPhase), divRho(NbIncTotalPrimMax, NbPhase), SmRho(NbPhase)
      type(TYPE_LoisThermoHydro) lthydro

      !Init outpus to zero
      vshydro_s%edata%Svel = 0.d0
      vshydro_s%edata%divSvel_ss = 0.d0
      vshydro_s%edata%SmSvel = 0.d0
      vshydro_h%divQwVsHead = 0.d0

      lthydro = LoisThermoHydroMSWell(well_idx_s)

      if (inc_s%ic == GAS_CONTEXT) then

         zeta_s = lthydro%DensiteMolaire(GAS_PHASE)
         Um = qw/zeta_s

         vshydro_s%edata%MixVel = Um
         vshydro_s%edata%Svel(GAS_PHASE) = Um
         vshydro_s%edata%SmSvel(GAS_PHASE) = -lthydro%SmDensiteMolaire(GAS_PHASE)*qw/zeta_s**2
         vshydro_s%edata%divSvel_ss(:, GAS_PHASE, 1) = -lthydro%divDensiteMolaire(:, GAS_PHASE)*qw/zeta_s**2
         vshydro_h%divQwVsHead(GAS_PHASE) = 1.d0/zeta_s

      else if (inc_s%ic == LIQUID_CONTEXT) then

         zeta_s = lthydro%DensiteMolaire(LIQUID_PHASE)
         Um = qw/zeta_s

         vshydro_s%edata%MixVel = Um
         vshydro_s%edata%Svel(LIQUID_PHASE) = Um
         vshydro_s%edata%SmSvel(LIQUID_PHASE) = -lthydro%SmDensiteMolaire(LIQUID_PHASE)*qw/zeta_s**2
         vshydro_s%edata%divSvel_ss(:, LIQUID_PHASE, 1) = -lthydro%divDensiteMolaire(:, LIQUID_PHASE)*qw/zeta_s**2
         vshydro_h%divQwVsHead(LIQUID_PHASE) = 1.d0/zeta_s

      else if (inc_s%ic == DIPHASIC_CONTEXT) then

         gas_sat = inc_s%Saturation(GAS_PHASE)
         Zeta(GAS_PHASE) = lthydro%DensiteMolaire(GAS_PHASE)
         Zeta(LIQUID_PHASE) = lthydro%DensiteMolaire(LIQUID_PHASE)
         divZeta(:, GAS_PHASE) = lthydro%divDensiteMolaire(:, GAS_PHASE)
         divZeta(:, LIQUID_PHASE) = lthydro%divDensiteMolaire(:, LIQUID_PHASE)
         SmZeta(GAS_PHASE) = lthydro%SmDensiteMolaire(GAS_PHASE)
         SmZeta(LIQUID_PHASE) = lthydro%SmDensiteMolaire(LIQUID_PHASE)

         Rho(GAS_PHASE) = lthydro%DensiteMassique(GAS_PHASE)
         Rho(LIQUID_PHASE) = lthydro%DensiteMassique(LIQUID_PHASE)
         divRho(:, GAS_PHASE) = lthydro%divDensiteMassique(:, GAS_PHASE)
         divRho(:, LIQUID_PHASE) = lthydro%divDensiteMassique(:, LIQUID_PHASE)
         SmRho(GAS_PHASE) = lthydro%SmDensiteMassique(GAS_PHASE)
         SmRho(LIQUID_PHASE) = lthydro%SmDensiteMassique(LIQUID_PHASE)

         Uc = (sigmagl*gravity*(Rho(LIQUID_PHASE) - Rho(GAS_PHASE))/Rho(LIQUID_PHASE)**2)**0.25d0

         ss = 0.25d0*((sigmagl*gravity)**0.25d0) &
              *(1.d0/Rho(LIQUID_PHASE) - Rho(GAS_PHASE)/Rho(LIQUID_PHASE)**2)**(-0.75d0)

         dRhoUc(LIQUID_PHASE) = ss*(2.d0*Rho(GAS_PHASE) - Rho(LIQUID_PHASE))/Rho(LIQUID_PHASE)**3
         dRhoUc(GAS_PHASE) = -ss/Rho(LIQUID_PHASE)**2

         divUc(:) = dRhoUc(LIQUID_PHASE)*divRho(:, LIQUID_PHASE)
         divUc(:) = divUc(:) + dRhoUc(GAS_PHASE)*divRho(:, GAS_PHASE)

         SmUc = dRhoUc(LIQUID_PHASE)*SmRho(LIQUID_PHASE) + dRhoUc(GAS_PHASE)*SmRho(GAS_PHASE)

         ratiod = Rho(GAS_PHASE)/Rho(LIQUID_PHASE)
         divratiod(:) = divRho(:, GAS_PHASE)/Rho(LIQUID_PHASE)
         divratiod(:) = divratiod(:) - divRho(:, LIQUID_PHASE)*Rho(GAS_PHASE)/Rho(LIQUID_PHASE)**2
         Smratiod = SmRho(GAS_PHASE)/Rho(LIQUID_PHASE) - SmRho(LIQUID_PHASE)*Rho(GAS_PHASE)/Rho(LIQUID_PHASE)**2

         sC0 = gas_sat*DFMHydroMSWells_f_C0(gas_sat)
         divsC0(:) = DFMHydroMSWells_df_sgC0(gas_sat)*lthydro%divSaturation(:, GAS_PHASE)

         sC0K = DFMHydroMSWells_f_sgC0K(gas_sat)
         divsC0K(:) = DFMHydroMSWells_df_sgC0K(gas_sat)*lthydro%divSaturation(:, GAS_PHASE)

         fG = DFMHydroMSWells_f_G(gas_sat, ratiod)

         call DFMHydroMSWells_df_G(gas_sat, ratiod, dsG, dratioG)

         sUd = sC0K*fG*Uc

         divG(:) = dsG*lthydro%divSaturation(:, GAS_PHASE) &
                   + dratioG*divratiod(:)
         SmG = dratioG*Smratiod

         divsUd(:) = divsC0K(:)*fG*Uc &
                     + sC0K*Uc*divG(:) + sC0K*fG*divUc(:)

         SmsUd = sC0K*Uc*SmG + sC0K*fG*SmUc

         ss = Zeta(LIQUID_PHASE) - (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*sC0

         Um = (qw + (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*sUd)/ss

         divUm(:) = ((divZeta(:, LIQUID_PHASE) - divZeta(:, GAS_PHASE))*sUd &
                     + (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*divsUd(:))/ss

         divUm(:) = divUm(:) - (qw + (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*sUd) &
                    *(divZeta(:, LIQUID_PHASE) - (divZeta(:, LIQUID_PHASE) - divZeta(:, GAS_PHASE))*sC0 &
                      - (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*divsC0(:))/ss**2

         SmUm = ((SmZeta(LIQUID_PHASE) - SmZeta(GAS_PHASE))*sUd &
                 + (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*SmsUd)/ss

         SmUm = SmUm - (qw + (Zeta(LIQUID_PHASE) - Zeta(GAS_PHASE))*sUd) &
                *(SmZeta(LIQUID_PHASE) - (SmZeta(LIQUID_PHASE) - SmZeta(GAS_PHASE))*sC0)/ss**2

         Vsg = sC0*Um + sUd
         vshydro_h%divQwVsHead(GAS_PHASE) = sC0/ss
         vshydro_s%edata%divSvel_ss(:, GAS_PHASE, 1) = divsC0(:)*Um + sC0*divUm(:) + divsUd(:)
         vshydro_s%edata%SmSvel(GAS_PHASE) = sC0*SmUm + SmsUd

         Vsl = Um - Vsg
         vshydro_h%divQwVsHead(LIQUID_PHASE) = (1.d0 - sC0)/ss
         vshydro_s%edata%divSvel_ss(:, LIQUID_PHASE, 1) = -divsC0(:)*Um + (1 - sC0)*divUm(:) - divsUd(:)
         vshydro_s%edata%SmSvel(LIQUID_PHASE) = (1.d0 - sC0)*SmUm - SmsUd

         vshydro_s%edata%MixVel = Um
         vshydro_s%edata%Svel(GAS_PHASE) = Vsg
         vshydro_s%edata%Svel(LIQUID_PHASE) = Vsl

      end if

#endif

   end subroutine VSHydroMSWells_compute_bdc_prod_cv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Init Pressure from mean velocity using friction law
!     and mean density/viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine VSHydroMSWells_init_cv( &
      inc_s, &
      inc_sp, &
      well_idx_s, &
      well_idx_sp, &
      scoords_idx_s, &
      scoords_idx_sp, &
      cross_section, &
      eorientation, &
      vshydro_s)

      type(TYPE_IncCVReservoir), intent(inout)  :: inc_s
      type(TYPE_IncCVReservoir), intent(in)  ::  inc_sp
      double precision, intent(in)           :: cross_section, eorientation
      integer, intent(in)                    :: well_idx_s, well_idx_sp, scoords_idx_s, scoords_idx_sp
      type(TYPE_VSHydroMSWells), intent(in)  :: vshydro_s
      !tmp,
      !TODO: we are copying things here (not so efficient), probably constructing an array of pointers is better
      type(TYPE_IncCVReservoir) :: inc_ss(2)
      type(ContextInfo)         :: ctxinfo_ss(2)
      integer                   :: well_idx_ss(2)
      double precision :: &
         Rhom, & ! arithmetic mean of the mixture
         divRhom_ss(NbIncTotalPrimMax, 2), &
         SmRhom, &
         Viscom, &  ! arithmetic mean of the mixture
         divViscom_ss(NbIncTotalPrimMax, 2), &
         SmViscom
      !tmp
      double precision :: DPhi, SizeEdge, Redge, s, Um, zk_s, zk_sp
      double precision, parameter :: pi = datan(1.0d0)*4.0d0
      type(TYPE_LoisThermoHydro) lthydro

      ! init tmp values for each cv

      inc_ss(1) = inc_s  !node s
      inc_ss(2) = inc_sp !parent node of s
      well_idx_ss(1) = well_idx_s
      well_idx_ss(2) = well_idx_sp

      !         and the second entry for the node s which is the parent of s_prime
      call LoisThermoHydro_init_cv(inc_ss(1), ctxinfo_ss(1))
      call LoisThermoHydro_init_cv(inc_ss(2), ctxinfo_ss(2))

      ! densite massique
      call VSHydroMSWells_ammix_densitemassique_cv( &
         inc_ss, &
         ctxinfo_ss, &
         well_idx_ss, &
         Rhom, &
         divRhom_ss, &
         SmRhom)

      ! viscosite
      call VSHydroMSWells_ammix_viscosite_cv( &
         inc_ss, &
         ctxinfo_ss, &
         well_idx_ss, &
         Viscom, &
         divViscom_ss, &
         SmViscom)

      Redge = dsqrt(cross_section/pi)
      s = (XNodeLocal(1, scoords_idx_s) - XNodeLocal(1, scoords_idx_sp))**2
      s = s + (XNodeLocal(2, scoords_idx_s) - XNodeLocal(2, scoords_idx_sp))**2
      s = s + (XNodeLocal(3, scoords_idx_s) - XNodeLocal(3, scoords_idx_sp))**2
      SizeEdge = dsqrt(s)

      Um = vshydro_s%edata%MixVel
      call DFMHydroMSWells_DPhi(Um, Rhom, Viscom, Redge, SizeEdge, DPhi)

      zk_s = XNodeLocal(3, scoords_idx_s)    ! z-cordinate of node s
      zk_sp = XNodeLocal(3, scoords_idx_sp) ! z-cordinate of parent of s

      inc_s%Pression = inc_sp%Pression - DPhi + gravity*Rhom*(zk_sp - zk_s)
      !Set phase pressure for all phases
      inc_s%phase_pressure(:) = inc_s%Pression

   end subroutine VSHydroMSWells_init_cv

   subroutine VSHydroMSWells_compute_cv( &
      inc_s, &
      inc_sp, &
      well_idx_s, &
      well_idx_sp, &
      scoords_idx_s, &
      scoords_idx_sp, &
      cross_section, &
      eorientation, &
      vshydro_s &
      )

      type(TYPE_IncCVReservoir), intent(in)  :: inc_s, inc_sp
      double precision, intent(in)           :: cross_section, eorientation
      integer, intent(in)                     :: well_idx_s, well_idx_sp, scoords_idx_s, scoords_idx_sp
      type(TYPE_VSHydroMSWells), intent(inout) :: vshydro_s
      !tmp,
      !TODO: we are copying things here (not so efficient), probably constructing an array of pointers is better
      type(TYPE_IncCVReservoir) :: inc_ss(2)
      type(ContextInfo)         :: ctxinfo_ss(2)
      integer                   :: well_idx_ss(2)
      double precision :: &
         Rhom, & ! arithmetic mean of the mixture
         divRhom_ss(NbIncTotalPrimMax, 2), &
         SmRhom, &
         Viscom, &  ! arithmetic mean of the mixture
         divViscom_ss(NbIncTotalPrimMax, 2), &
         SmViscom

#ifndef ComPASS_DIPHASIC_CONTEXT

      call CommonMPI_abort("Multi-segmented wells are only implemented for diphasic physics!")

#else

      double precision :: DPhi, Pk_s, Pk_sp, SmDphi, zk_s, zk_sp, s, SizeEdge, Redge
      double precision :: divDPhi_ss(NbIncTotalPrimMax, 2), divUm_ss(NbIncTotalPrimMax, 2)
      double precision :: DDPhiUm, DRhoUm, DViscoUm, Um, SmUm, gas_sat_s, gas_sat_sp
      double precision :: Rho(NbPhase), divRho_ss(NbIncTotalPrimMax, NbPhase, 2), SmRho(NbPhase)
      double precision, parameter :: pi = datan(1.0d0)*4.0d0
      double precision :: flux_Vsg, dUm_flux, drhol_flux, drhog_flux, dsigma_flux, flux_Vsl, dsg_flux_ss(2)
      type(TYPE_LoisThermoHydro) lthydro
      integer          ::node_idx

      ! init tmp values for each cv

      inc_ss(1) = inc_s  !node s
      inc_ss(2) = inc_sp !parent node of s
      well_idx_ss(1) = well_idx_s
      well_idx_ss(2) = well_idx_sp

      !         and the second entry for the node s which is the parent of s_prime
      call LoisThermoHydro_init_cv(inc_ss(1), ctxinfo_ss(1))
      call LoisThermoHydro_init_cv(inc_ss(2), ctxinfo_ss(2))

      ! densite massique
      call VSHydroMSWells_ammix_densitemassique_cv( &
         inc_ss, &
         ctxinfo_ss, &
         well_idx_ss, &
         Rhom, &
         divRhom_ss, &
         SmRhom)

      ! viscosite
      call VSHydroMSWells_ammix_viscosite_cv( &
         inc_ss, &
         ctxinfo_ss, &
         well_idx_ss, &
         Viscom, &
         divViscom_ss, &
         SmViscom)

      Pk_s = inc_ss(1)%Pression
      Pk_sp = inc_ss(2)%Pression

      zk_s = XNodeLocal(3, scoords_idx_s)    ! z-cordinate of node s
      zk_sp = XNodeLocal(3, scoords_idx_sp) ! z-cordinate of parent of s

      DPhi = Pk_sp - Pk_s + gravity*(zk_sp - zk_s)*Rhom

      divDPhi_ss(:, 1) = gravity*(zk_sp - zk_s)*divRhom_ss(:, 1)
      divDPhi_ss(1, 1) = divDPhi_ss(1, 1) - 1.d0
      divDPhi_ss(:, 2) = gravity*(zk_sp - zk_s)*divRhom_ss(:, 2)
      divDPhi_ss(1, 2) = divDPhi_ss(1, 2) + 1.d0
      SmDPhi = gravity*SmRhom*(zk_sp - zk_s)

      Redge = dsqrt(cross_section/pi)
      s = (XNodeLocal(1, scoords_idx_s) - XNodeLocal(1, scoords_idx_sp))**2
      s = s + (XNodeLocal(2, scoords_idx_s) - XNodeLocal(2, scoords_idx_sp))**2
      s = s + (XNodeLocal(3, scoords_idx_s) - XNodeLocal(3, scoords_idx_sp))**2
      SizeEdge = dsqrt(s)

      call DFMHydroMSWells_MixtureVelocity(DPhi, Rhom, Viscom, Redge, SizeEdge, &
                                           vshydro_s%edata%MixVel, &
                                           Um, &
                                           DDPhiUm, DRhoUm, DViscoUm)

      vshydro_s%edata%MixVel = Um

      divUm_ss(:, 1) = DRhoUm*divRhom_ss(:, 1) &
                       + DViscoUm*divViscom_ss(:, 1) + DDPhiUm*divDPhi_ss(:, 1)

      divUm_ss(:, 2) = DRhoUm*divRhom_ss(:, 2) &
                       + DViscoUm*divViscom_ss(:, 2) + DDPhiUm*divDPhi_ss(:, 2)

      SmUm = DRhoUm*SmRhom + DViscoUm*SmViscom + DDPhiUm*SmDphi

      !Init outputs to Zero
      vshydro_s%edata%Svel = 0.d0
      vshydro_s%edata%divSvel_ss = 0.d0
      vshydro_s%edata%SmSvel = 0.d0

      if ((inc_s%ic == GAS_CONTEXT) .and. (inc_sp%ic == GAS_CONTEXT)) then

         vshydro_s%edata%Svel(GAS_PHASE) = Um
         vshydro_s%edata%SmSvel(GAS_PHASE) = SmUm
         vshydro_s%edata%divSvel_ss(:, GAS_PHASE, :) = divUm_ss(:, :)

      else if ((inc_s%ic == LIQUID_CONTEXT) .and. (inc_sp%ic == LIQUID_CONTEXT)) then

         vshydro_s%edata%Svel(LIQUID_PHASE) = Um
         vshydro_s%edata%SmSvel(LIQUID_PHASE) = SmUm
         vshydro_s%edata%divSvel_ss(:, LIQUID_PHASE, :) = divUm_ss(:, :)

      else
         call VSHydroMSWells_amean_densitemassique_cv( &
            inc_ss, &
            ctxinfo_ss, &
            well_idx_ss, &
            LIQUID_PHASE, &
            Rho(LIQUID_PHASE), &
            divRho_ss(:, LIQUID_PHASE, :), &
            SmRho(LIQUID_PHASE))

         call VSHydroMSWells_amean_densitemassique_cv( &
            inc_ss, &
            ctxinfo_ss, &
            well_idx_ss, &
            GAS_PHASE, &
            Rho(GAS_PHASE), &
            divRho_ss(:, GAS_PHASE, :), &
            SmRho(GAS_PHASE))

         gas_sat_s = inc_s%Saturation(GAS_PHASE)
         gas_sat_sp = inc_sp%Saturation(GAS_PHASE)

         call DFMHydroMSWells_FluxVsg(gas_sat_s, gas_sat_sp, Um, &
                                      Rho(LIQUID_PHASE), Rho(GAS_PHASE), sigmagl, gravity, eorientation, &
                                      flux_Vsg, &
                                      dsg_flux_ss(1), dsg_flux_ss(2), dUm_flux, drhol_flux, drhog_flux, dsigma_flux)
         flux_Vsl = Um - flux_Vsg

         vshydro_s%edata%Svel(GAS_PHASE) = flux_Vsg
         vshydro_s%edata%Svel(LIQUID_PHASE) = flux_Vsl

         do node_idx = 1, 2 !node s_prime and  its parent

            lthydro = LoisThermoHydroMSWell(well_idx_ss(node_idx))

            vshydro_s%edata%divSvel_ss(:, GAS_PHASE, node_idx) = &
               +dsg_flux_ss(node_idx)*lthydro%divSaturation(:, GAS_PHASE) &
               + dUm_flux*divUm_ss(:, node_idx) &
               + drhol_flux*divRho_ss(:, LIQUID_PHASE, node_idx) &
               + drhog_flux*divRho_ss(:, GAS_PHASE, node_idx)
         end do

         vshydro_s%edata%SmSvel(GAS_PHASE) = dUm_flux*SmUm &
                                             + drhog_flux*SmRho(GAS_PHASE) &
                                             + drhol_flux*SmRho(LIQUID_PHASE)

         vshydro_s%edata%divSvel_ss(:, LIQUID_PHASE, :) = &
            divUm_ss(:, :) - vshydro_s%edata%divSvel_ss(:, GAS_PHASE, :)

         vshydro_s%edata%SmSvel(LIQUID_PHASE) = SmUm &
                                                - vshydro_s%edata%SmSvel(GAS_PHASE)

      end if

#endif

   end subroutine VSHydroMSWells_compute_cv

   !Computes the value of the   arithmetic mean of the mixture
   !For all input/output arguments withe postfix _ss
   !indicates that the 1-index of the last entry is w.r s_prime and the 2-index is w.r the node s (the parent of s_prime)
   subroutine VSHydroMSWells_ammix_densitemassique_cv(inc_ss, ctxinfo_ss, well_idx_ss, &
                                                      val, dval_ss, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc_ss(2)
      type(ContextInfo), intent(in) :: ctxinfo_ss(2)
      integer, intent(in) :: well_idx_ss(2)

      ! output
      double precision, intent(out) :: val
      double precision, intent(out) :: dval_ss(NbIncTotalPrimMax, 2)
      double precision, intent(out) :: Smval

      ! tmp
      double precision :: f, sat_iph
      integer :: iph, i, node_idx
      type(TYPE_LoisThermoHydro) lthydro

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val = 0.d0
      dval_ss(:, :) = 0.d0
      Smval = 0.d0

      do node_idx = 1, 2 !node s_prime and  its parent

         do i = 1, ctxinfo_ss(node_idx)%NbPhasePresente
            iph = ctxinfo_ss(node_idx)%NumPhasePresente(i)

            sat_iph = inc_ss(node_idx)%Saturation(iph)

            lthydro = LoisThermoHydroMSWell(well_idx_ss(node_idx))

            f = lthydro%DensiteMassique(iph)
            val = val + 0.5*sat_iph*f ! val

            dval_ss(:, node_idx) = dval_ss(:, node_idx) &
                                   + 0.5*sat_iph*lthydro%divDensiteMassique(:, iph) &
                                   + 0.5*f*lthydro%divSaturation(:, iph)

            Smval = Smval &
                    + 0.5*sat_iph*lthydro%SmDensiteMassique(iph)

         end do

      end do

   end subroutine VSHydroMSWells_ammix_densitemassique_cv

   !Computes the value of the   arithmetic mean of the mixture
   !For all input/output arguments withe postfix _ss
   !indicates that the 1-index of the last entry is w.r s_prime and the 2-index is w.r the node s (the parent of s_prime)
   subroutine VSHydroMSWells_ammix_viscosite_cv(inc_ss, ctxinfo_ss, well_idx_ss, &
                                                val, dval_ss, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc_ss(2)
      type(ContextInfo), intent(in) :: ctxinfo_ss(2)
      integer, intent(in) :: well_idx_ss(2)

      ! output
      double precision, intent(out) :: val
      double precision, intent(out) :: dval_ss(NbIncTotalPrimMax, 2)
      double precision, intent(out) :: Smval

      ! tmp
      double precision :: f, sat_iph
      integer :: iph, i, node_idx
      type(TYPE_LoisThermoHydro) lthydro

      ! 1. val
      ! 2. dval
      ! 3. Smval

      val = 0.d0
      dval_ss(:, :) = 0.d0
      Smval = 0.d0

      do node_idx = 1, 2 !node s_prime and  its parent

         do i = 1, ctxinfo_ss(node_idx)%NbPhasePresente
            iph = ctxinfo_ss(node_idx)%NumPhasePresente(i)

            sat_iph = inc_ss(node_idx)%Saturation(iph)

            lthydro = LoisThermoHydroMSWell(well_idx_ss(node_idx))

            f = lthydro%Viscosite(iph)
            val = val + 0.5*sat_iph*f ! val

            dval_ss(:, node_idx) = dval_ss(:, node_idx) &
                                   + 0.5*sat_iph*lthydro%divViscosite(:, iph) &
                                   + 0.5*f*lthydro%divSaturation(:, iph)

            Smval = Smval &
                    + 0.5*sat_iph*lthydro%SmViscosite(iph)

         end do

      end do

   end subroutine VSHydroMSWells_ammix_viscosite_cv

   !Computes the value of the mean  arithmetic
   !It assumes we have only the phase with PHASE_ID in only one node (fahter xor son) of the edge
   subroutine VSHydroMSWells_amean_densitemassique_cv(inc_ss, ctxinfo_ss, &
                                                      well_idx_ss, PHASE_ID, &
                                                      val, dval_ss, Smval)

      ! input
      type(TYPE_IncCVReservoir), intent(in)  :: inc_ss(2)
      type(ContextInfo), intent(in) :: ctxinfo_ss(2)
      integer, intent(in) :: well_idx_ss(2)
      integer, intent(in) :: PHASE_ID

      ! output
      double precision, intent(out) :: val
      double precision, intent(out) :: dval_ss(NbIncTotalPrimMax, 2)
      double precision, intent(out) :: Smval

      ! tmp
      logical :: phases_both_nodes
      double precision ::  f, sats_sum
      integer :: i, j, iph, ns_idx, nsprime_idx, node_idx, iph_prime
      type(TYPE_LoisThermoHydro) lthydro

      val = 0.d0
      Smval = 0.d0
      dval_ss = 0.d0

      phases_both_nodes = .false.
      nsprime_idx = 1!idx of node s_prime
      ns_idx = 2    !idx of parent node of s_sprime

      !Let us first find out when a phase is in both nodes
      !do the outer loop with respect to node s_prime
      do i = 1, ctxinfo_ss(nsprime_idx)%NbPhasePresente
         iph_prime = ctxinfo_ss(nsprime_idx)%NumPhasePresente(i)

         if (iph_prime == PHASE_ID) then

            sats_sum = inc_ss(nsprime_idx)%Saturation(iph_prime)

            do j = 1, ctxinfo_ss(ns_idx)%NbPhasePresente
               iph = ctxinfo_ss(ns_idx)%NumPhasePresente(j)

               if (iph_prime .eq. iph) then !phase is in both nodes
                  phases_both_nodes = .true.
                  exit
               end if
            end do
         end if
      end do

      if (phases_both_nodes) then !Perform computation when PHASE_ID is in both nodes

         iph = PHASE_ID

         do node_idx = 1, 2 !node s_prime and  its parent

            lthydro = LoisThermoHydroMSWell(well_idx_ss(node_idx))

            f = lthydro%DensiteMassique(iph)
            val = val + 0.5*f ! val

            dval_ss(:, node_idx) = dval_ss(:, node_idx) &
                                   + 0.5*lthydro%divDensiteMassique(:, iph)

            Smval = Smval &
                    + 0.5*lthydro%SmDensiteMassique(iph)

         end do !node_idx

      else !phases_both_nodes, PHASE_ID is only in one node

         do node_idx = 1, 2 !node s_prime and  its parent

            do i = 1, ctxinfo_ss(node_idx)%NbPhasePresente
               iph = ctxinfo_ss(node_idx)%NumPhasePresente(i)

               if (iph .eq. PHASE_ID) then

                  lthydro = LoisThermoHydroMSWell(well_idx_ss(node_idx))

                  f = lthydro%DensiteMassique(iph)
                  val = val + 0.5*f ! val

                  dval_ss(:, node_idx) = dval_ss(:, node_idx) &
                                         + 0.5*lthydro%divDensiteMassique(:, iph)

                  Smval = Smval &
                          + 0.5*lthydro%SmDensiteMassique(iph)

                  return !Nothing to be done, we have found the node who has PHASE_ID

               end if!PHASE_ID
            end do!i
         end do!node_idx

      end if

   end subroutine VSHydroMSWells_amean_densitemassique_cv

end module VSHydroMSWells
