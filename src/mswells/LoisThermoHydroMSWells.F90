!
! This file is part of ComPASS.
!
! ComPAsSS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module LoisThermoHydroMSWells
   use, intrinsic :: iso_c_binding, only: c_double
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use CommonMPI
   use Thermodynamics
   use DefModel
   use NumbyContext
   use IncCVReservoir
   use IncCVMSWells
   use IncPrimSecdMSWells
   use MeshSchema
   use LoisThermoHydro

#else

   use CommonMPI, only: commRank
   use DefModel, only: &
      NbPhase, NbComp, IndThermique, &
      NbIncTotalPrimMax
   use NumbyContext, only: &
      NumIncPTC2NumIncComp_comp_ctx, NumIncPTC2NumIncComp_phase_ctx, &
      NbIncPTC_ctx
   use IncCVReservoir, only: &
      TYPE_IncCVReservoir

   use IncCVMSWells, only: &
      IncMSWell

   use IncPrimSecdMSWells, only: &
      TYPE_IncPrimSecd, &
      IncPrimSecdMSWell

   use MeshSchema, only: &
      NodebyMSWellLocal, &
      NbMSWellLocal_Ncpus

   use LoisThermoHydro, only: &
      ContextInfo, &
      LoisThermoHydro_init_cv, &
      LoisThermoHydro_viscosite_cv, &
      LoisThermoHydro_densitemassique_cv, &
      LoisThermoHydro_densitemolaire_cv, &
      LoisThermoHydro_Inc_cv, &
      LoisThermoHydro_Saturation_cv, &
      LoisThermoHydro_EnergieInterne_cv, &
      LoisThermoHydro_Enthalpie_cv
#endif

   implicit none

   type :: TYPE_LoisThermoHydro
      ! Rq important: it contains values for all phases, not only present phases
      ! densite massique
      double precision :: &
         DensiteMassique(NbPhase), &
         divDensiteMassique(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMassique(NbPhase)

      ! pression
      double precision :: &
         divPression(NbIncTotalPrimMax), &
         SmPression

      ! Saturation
      double precision ::  divSaturation(NbIncTotalPrimMax, NbPhase)

      ! the following  must be allocated even if there is no energy transfer
      double precision ::  SmTemperature

#ifdef _THERMIQUE_
      ! temperature
      double precision ::  divTemperature(NbIncTotalPrimMax)
#endif

      ! densite molaire
      double precision  :: &
         DensiteMolaire(NbPhase), &
         divDensiteMolaire(NbIncTotalPrimMax, NbPhase), &
         SmDensiteMolaire(NbPhase)

      !  viscosity
      double precision :: &
         Viscosite(NbPhase), &
         divViscosite(NbIncTotalPrimMax, NbPhase), &
         SmViscosite(NbPhase)

      ! Comp
      double precision  :: &
         divComp(NbIncTotalPrimMax, NbComp, NbPhase), &
         SmComp(NbComp, NbPhase)

      ! Energie & enthalpie
      double precision  :: &
         EnergieInterne(NbPhase), &
         divEnergieInterne(NbIncTotalPrimMax, NbPhase), &
         SmEnergieInterne(NbPhase), &
         Enthalpie(NbPhase), &
         divEnthalpie(NbIncTotalPrimMax, NbPhase), &
         SmEnthalpie(NbPhase)

   end type TYPE_LoisThermoHydro

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   type(TYPE_LoisThermoHydro), allocatable, dimension(:), target, protected :: &
      LoisThermoHydroMSWell !<  Injector/Prodicer-MSWells-LoisThermoHydro

   public :: &
      LoisThermoHydroMSWells_allocate, &
      LoisThermoHydroMSWells_free, &
      LoisThermoHydroMSWells_compute, &
      LoisThermoHydroMSWells_compute_cv

contains

   subroutine LoisThermoHydroMSWells_allocate()

      integer :: Nb, Nnz

      Nb = NodebyMSWellLocal%Nb
      Nnz = NodebyMSWellLocal%Pt(Nb + 1)
      allocate (LoisThermoHydroMSWell(Nnz))

   end subroutine LoisThermoHydroMSWells_allocate

   subroutine LoisThermoHydroMSWells_free()

      deallocate (LoisThermoHydroMSWell)

   end subroutine LoisThermoHydroMSWells_free

   subroutine LoisThermoHydroMSWells_compute() &
      bind(C, name="LoisThermoHydroMSWells_compute")

      integer :: s, k, nbwells

      !For each mswell
      nbwells = NbMSWellLocal_Ncpus(commRank + 1)
      do k = 1, nbwells
         ! looping from  queue to  head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)

            call LoisThermoHydroMSWells_compute_cv( &
               IncMSWell(s)%coats, &
               IncPrimSecdMSWell(s), &
               LoisThermoHydroMSWell(s))

         end do
      end do

   end subroutine LoisThermoHydroMSWells_compute

   subroutine LoisThermoHydroMSWells_compute_cv(inc, inc_primsecd, lois_thydro)

      type(TYPE_IncCVReservoir), intent(in) :: inc
      type(TYPE_IncPrimSecd)     :: inc_primsecd
      type(TYPE_LoisThermoHydro), intent(out)::lois_thydro

      integer ::  i, icp, iph, k, context
      real(c_double) :: dpadS(NbPhase)
      type(ContextInfo) :: ctxinfo

      dpadS = 0.d0 !Derivative of phase pressure with respect saturation
      ! init tmp values for each cv
      call LoisThermoHydro_init_cv(inc, ctxinfo)

      ! viscosite
      call LoisThermoHydro_viscosite_cv(inc, dpadS, &
                                        inc_primsecd%dXssurdXp, &
                                        inc_primsecd%SmdXs, &
                                        inc_primsecd%NumIncTotalPrimCV, &
                                        inc_primsecd%NumIncTotalSecondCV, &
                                        lois_thydro%Viscosite, &
                                        lois_thydro%divViscosite, &
                                        lois_thydro%SmViscosite, &
                                        .true.) !Use absolute ordering

      ! densite massique
      call LoisThermoHydro_densitemassique_cv(inc, dpadS, &
                                              inc_primsecd%dXssurdXp, &
                                              inc_primsecd%SmdXs, &
                                              inc_primsecd%NumIncTotalPrimCV, &
                                              inc_primsecd%NumIncTotalSecondCV, &
                                              lois_thydro%DensiteMassique, &
                                              lois_thydro%divDensiteMassique, &
                                              lois_thydro%SmDensiteMassique)!Here absolute ordering is by default

      ! deniste molaire
      call LoisThermoHydro_densitemolaire_cv(inc, dpadS, &
                                             inc_primsecd%dXssurdXp, &
                                             inc_primsecd%SmdXs, &
                                             inc_primsecd%NumIncTotalPrimCV, &
                                             inc_primsecd%NumIncTotalSecondCV, &
                                             lois_thydro%DensiteMolaire, &
                                             lois_thydro%divDensiteMolaire, &
                                             lois_thydro%SmDensiteMolaire, &
                                             .true.) !Use absolute ordering

      ! Reference Pressure (unknown index is 1)
      call LoisThermoHydro_Inc_cv(1, inc, &
                                  inc_primsecd%NumIncTotalPrimCV, &
                                  inc_primsecd%NumIncTotalSecondCV, &
                                  inc_primsecd%dXssurdXp, &
                                  inc_primsecd%SmdXs, &
                                  lois_thydro%divPression, &
                                  lois_thydro%SmPression)

      ! Comp
      context = inc%ic
      do i = 2 + IndThermique, NbIncPTC_ctx(context) ! loop over index of Components
         iph = NumIncPTC2NumIncComp_phase_ctx(i, context) ! phase corresponding to unknown i
         icp = NumIncPTC2NumIncComp_comp_ctx(i, context) ! component corresponding to unknown i
         call LoisThermoHydro_Inc_cv(i, inc, &
                                     inc_primsecd%NumIncTotalPrimCV, inc_primsecd%NumIncTotalSecondCV, &
                                     inc_primsecd%dXssurdXp, inc_primsecd%SmdXs, &
                                     lois_thydro%divComp(:, icp, iph), lois_thydro%SmComp(icp, iph))
      enddo

      ! Saturation div FIXME: not done with LoisThermoHydro_Inc_cv because last saturation is eliminated
      call LoisThermoHydro_Saturation_cv(inc, ctxinfo, &
                                         inc_primsecd%NumIncTotalPrimCV, &
                                         inc_primsecd%NumIncTotalSecondCV, &
                                         inc_primsecd%dXssurdXp, &
                                         lois_thydro%divSaturation, &
                                         .true.) !Use absolute ordering
#ifdef _THERMIQUE_
      ! FIXME: Temperature (unknown index is 2)
      call LoisThermoHydro_Inc_cv(2, inc, &
                                  inc_primsecd%NumIncTotalPrimCV, &
                                  inc_primsecd%NumIncTotalSecondCV, &
                                  inc_primsecd%dXssurdXp, &
                                  inc_primsecd%SmdXs, &
                                  lois_thydro%divTemperature, &
                                  lois_thydro%SmTemperature)

      ! energie interne
      call LoisThermoHydro_EnergieInterne_cv(inc, dpadS, ctxinfo, &
                                             inc_primsecd%dXssurdXp, &
                                             inc_primsecd%SmdXs, &
                                             inc_primsecd%NumIncTotalPrimCV, &
                                             inc_primsecd%NumIncTotalSecondCV, &
                                             lois_thydro%EnergieInterne, &
                                             lois_thydro%divEnergieInterne, &
                                             lois_thydro%SmEnergieInterne, &
                                             .true.) !Use absolute ordering
      ! Enthalpie
      call LoisThermoHydro_Enthalpie_cv(inc, dpadS, ctxinfo, &
                                        inc_primsecd%dXssurdXp, &
                                        inc_primsecd%SmdXs, &
                                        inc_primsecd%NumIncTotalPrimCV, &
                                        inc_primsecd%NumIncTotalSecondCV, &
                                        lois_thydro%Enthalpie, &
                                        lois_thydro%divEnthalpie, &
                                        lois_thydro%SmEnthalpie, &
                                        .true.) !Use absolute ordering
#endif

   end subroutine LoisThermoHydroMSWells_compute_cv

end module LoisThermoHydroMSWells
