!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module NumbyContext

   use iso_c_binding, only: c_bool
   use CommonType, only: ModelConfiguration

   implicit none

   ! C^tilde_Q: ensembles des composants ds les phases absentes fonction du contexte
   integer, dimension(:), allocatable, protected :: &
      NbCompCtilde_ctx
   integer, dimension(:, :), allocatable, protected :: &
      NumCompCtilde_ctx

   ! ***** Equation ***** !

   ! info of problem size
   integer, dimension(:), allocatable, protected :: &
      NbEqFermeture_ctx

   ! f_i^alpha * C_i^alpha = f_i^beta * C_i^beta, for Q
   integer, dimension(:), allocatable, protected :: &
      NbEqEquilibre_ctx ! Number of equillibrium eq fonction of the context

   integer, dimension(:, :), allocatable, protected :: &
      NumCompEqEquilibre_ctx     ! Index of comp involved in equillibrium eq

   integer, dimension(:, :, :), allocatable, protected :: &
      Num2PhasesEqEquilibre_ctx  ! Index of phases involved in equillibrium eq

   ! ***** Inc ***** !

   ! nb of IncPTC
   integer, dimension(:), allocatable, protected :: &
      NbIncPTC_ctx

   ! nb of total unknowns P (T) C S ... (whithout Ctilde)
   integer, dimension(:), allocatable, protected :: &
      NbIncTotal_ctx

   ! nb of primary unknowns (whithout Ctilde) (= NbComp + IndThermique - NbCompCtilde_ctx)
   integer, dimension(:), allocatable, protected :: &
      NbIncTotalPrim_ctx

   ! from IncComp to IncPTC (num)
   integer, dimension(:, :, :), allocatable, target :: &
      NumIncComp2NumIncPTC_ctx

   ! from IncPTC to IncComp (i and alpha)
   integer, dimension(:, :), allocatable, protected :: &
      NumIncPTC2NumIncComp_comp_ctx, &
      NumIncPTC2NumIncComp_phase_ctx

   public :: &
      NumbyContext_make, &
      NumbyContext_free

   private :: &
      NumbyContext_is_phase_present, &
      NumbyContext_PhaseComp, &
      NumbyContext_Inc, &
      NumbyContext_Eq

contains

   ! main subroutine of this module
   subroutine NumbyContext_make(configuration)
      type(ModelConfiguration), intent(in) :: configuration

      ! info phase and comp
      call NumbyContext_PhaseComp(configuration)

      ! info equation
      call NumbyContext_Eq(configuration)

      ! Inc
      call NumbyContext_Inc(configuration)

   end subroutine NumbyContext_make

   ! free everything
   subroutine NumbyContext_Free

      deallocate (NbCompCtilde_ctx)

      deallocate (NbIncPTC_ctx)
      deallocate (NbIncTotal_ctx)
      deallocate (NbIncTotalPrim_ctx)

      deallocate (NumIncComp2NumIncPTC_ctx)
      deallocate (NumIncPTC2NumIncComp_comp_ctx)
      deallocate (NumIncPTC2NumIncComp_phase_ctx)

      deallocate (Num2PhasesEqEquilibre_ctx)

      deallocate (NbEqEquilibre_ctx)
      deallocate (NbEqFermeture_ctx)

   end subroutine NumbyContext_Free

   logical(c_bool) function NumbyContext_is_phase_present(configuration, iph, ic) result(here)
      type(ModelConfiguration), intent(in) :: configuration
      integer, intent(in)  :: iph, ic
      integer :: k

      here = .false.
      do k = 1, configuration%NbPhasePresente_ctx(ic)
         if (configuration%NumPhasePresente_ctx(k, ic) == iph) then
            here = .true.
            exit
         endif
      end do

   end function NumbyContext_is_phase_present

   subroutine NumbyContext_PhaseComp(configuration)
      type(ModelConfiguration), intent(in) :: configuration
      integer :: NbPhase, NbComp, NbContexte
      integer :: ic, iph, icp, n
      logical :: IsCtidle

      NbPhase = configuration%nb_phases
      NbComp = configuration%nb_components
      NbContexte = configuration%nb_contexts

      ! 2. Ensembles Ctilde fct du contexte

      allocate (NbCompCtilde_ctx(NbContexte))
      allocate (NumCompCtilde_ctx(NbComp, NbContexte))
      NumCompCtilde_ctx(:, :) = 0

      do ic = 1, NbContexte

         n = 0
         do icp = 1, NbComp ! loop of components

            ! check if icp in C_tidle
            IsCtidle = .true.
            do iph = 1, NbPhase

               if (configuration%MCP(icp, iph) == 1 .and. &
                   NumbyContext_is_phase_present(configuration, iph, ic)) then
                  IsCtidle = .false.
                  exit
               end if
            end do

            if (IsCtidle .eqv. .true.) then
               n = n + 1
               NumCompCtilde_ctx(n, ic) = icp
            end if
         end do

         NbCompCtilde_ctx(ic) = n
      end do

   end subroutine NumbyContext_PhaseComp

   subroutine NumbyContext_Inc(configuration)
      type(ModelConfiguration), intent(in) :: configuration
      integer :: NbPhase, NbComp, NbContexte
      integer :: i, ic, iph, icp, n
      integer :: NbIncPTCMax

      NbPhase = configuration%nb_phases
      NbComp = configuration%nb_components
      NbContexte = configuration%nb_contexts
      NbIncPTCMax = configuration%NbIncPTCMax

      ! 1. Nb of unknowns
      allocate (NbIncPTC_ctx(NbContexte))    ! P (T) C depending on the context
      allocate (NbIncTotal_ctx(NbContexte))  ! total nb: P (T) C S_Principal ...

      do ic = 1, NbContexte

         n = 1 + configuration%IndThermique ! P (and T)

         do i = 1, configuration%NbPhasePresente_ctx(ic)
            iph = configuration%NumPhasePresente_ctx(i, ic) ! phase

            ! nb of components in iph, using MCP
            do icp = 1, NbComp

               if (configuration%MCP(icp, iph) == 1) then ! if component in phase iph
                  n = n + 1
               end if
            end do
         enddo

         NbIncPTC_ctx(ic) = n
         NbIncTotal_ctx(ic) = NbIncPTC_ctx(ic) + configuration%NbPhasePresente_ctx(ic) ! Warning: NbIncTotal_ctx=NbIncTotalPrim_ctx+NbEqFermeture + 1 because 1 saturation eliminated in hard
#ifdef ComPASS_WITH_diphasic_PHYSICS
         ! FIXME: Laurence triche pour avoir les molar flowrates comme inconnues supplémentaires
         if (ic > 3) then
            NbIncTotal_ctx(ic) = NbIncTotal_ctx(ic) + configuration%NbPhasePresente_ctx(ic)
         endif
#endif
      end do

      ! 2. Nb of primary unknowns
      allocate (NbIncTotalPrim_ctx(NbContexte))

      do ic = 1, NbContexte
         NbIncTotalPrim_ctx(ic) = NbComp + configuration%IndThermique - NbCompCtilde_ctx(ic)
      end do

      ! 3.1  from IncComp to IncPTC (num)
      allocate (NumIncComp2NumIncPTC_ctx(NbComp, NbPhase, NbContexte))

      ! 3.2  from IncPTC to IncComp (i and alpha)
      allocate (NumIncPTC2NumIncComp_comp_ctx(NbIncPTCMax, NbContexte))
      allocate (NumIncPTC2NumIncComp_phase_ctx(NbIncPTCMax, NbContexte))

      NumIncComp2NumIncPTC_ctx(:, :, :) = 0
      NumIncPTC2NumIncComp_comp_ctx(:, :) = 0
      NumIncPTC2NumIncComp_phase_ctx(:, :) = 0

      do ic = 1, NbContexte

         n = 1 + configuration%IndThermique ! Pref (and T)

         ! loop of phase
         do i = 1, configuration%NbPhasePresente_ctx(ic)
            iph = configuration%NumPhasePresente_ctx(i, ic)

            ! loop of component in phase iph
            do icp = 1, NbComp

               if (configuration%MCP(icp, iph) == 1) then ! if component in phase iph
                  n = n + 1
                  NumIncComp2NumIncPTC_ctx(icp, iph, ic) = n
                  NumIncPTC2NumIncComp_comp_ctx(n, ic) = icp
                  NumIncPTC2NumIncComp_phase_ctx(n, ic) = iph
               end if
            enddo

         end do
      enddo

   end subroutine NumbyContext_Inc

   subroutine NumbyContext_Eq(configuration)
      type(ModelConfiguration), intent(in) :: configuration
      integer :: NbPhase, NbComp, NbContexte
      integer :: ic, iph, icp, nphi, n, i
      integer :: PhPrComp(configuration%nb_phases)

      NbPhase = configuration%nb_phases
      NbComp = configuration%nb_components
      NbContexte = configuration%nb_contexts

      allocate (NumCompEqEquilibre_ctx(configuration%NbEqEquilibreMax, NbContexte))
      allocate (Num2PhasesEqEquilibre_ctx(2, configuration%NbEqEquilibreMax, NbContexte))

      allocate (NbEqEquilibre_ctx(NbContexte))
      allocate (NbEqFermeture_ctx(NbContexte))

      do ic = 1, NbContexte

         n = 0 ! used for NbEquilibre

         do icp = 1, NbComp

            nphi = 0
            PhPrComp(:) = 0

            ! PhPrComp: phase present and contains icp
            do iph = 1, NbPhase

               if (configuration%MCP(icp, iph) == 1 .and. &
                   NumbyContext_is_phase_present(configuration, iph, ic)) then
                  nphi = nphi + 1
                  PhPrComp(nphi) = iph
               end if
            end do

            if (nphi >= 2) then

               ! loop of phases in f_i^alpha * C_i^alpha = f_i^beta * C_i^beta, for Q
               do i = 1, nphi - 1
                  NumCompEqEquilibre_ctx(i + n, ic) = icp ! i
                  Num2PhasesEqEquilibre_ctx(1, i + n, ic) = PhPrComp(1)   ! alpha
                  Num2PhasesEqEquilibre_ctx(2, i + n, ic) = PhPrComp(i + 1) ! beta
               end do

               n = n + nphi - 1
            end if

         end do ! end of loop icp

         NbEqEquilibre_ctx(ic) = n
#ifdef ComPASS_WITH_diphasic_PHYSICS
         ! count the number of secd unknowns (=nb of closure laws), pssecd is filled in DefModel
         n = 0
         do i = 1, size(configuration%pssecd, 1)
            if (configuration%pssecd(i, ic) > 0) n = n + 1
         end do
         NbEqFermeture_ctx(ic) = n
#else
         NbEqFermeture_ctx(ic) = configuration%NbPhasePresente_ctx(ic) + n
#endif
      end do

   end subroutine NumbyContext_Eq

end module NumbyContext
