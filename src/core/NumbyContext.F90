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
  integer, dimension(:,:), allocatable, protected :: &
      NumCompCtilde_ctx


  ! ***** Equation ***** !

  ! info of problem size
  integer, dimension(:), allocatable, protected :: &
      NbEqFermeture_ctx

  ! f_i^alpha * C_i^alpha = f_i^beta * C_i^beta, for Q
  integer, dimension(:), allocatable, protected :: &
      NbEqEquilibre_ctx ! Nombre d'Equation d'Equilibre thermodynamique fct du contexte

  integer, dimension(:,:), allocatable, protected :: &
      NumCompEqEquilibre_ctx     ! Numero comp

  integer, dimension(:,:,:), allocatable, protected :: &
      Num2PhasesEqEquilibre_ctx  ! Numero du couple de phases


  ! ***** Inc ***** !

  ! nb of IncPTC
  integer, dimension(:), allocatable, protected :: &
      NbIncPTC_ctx

  ! nb of IncPTCSPrim
  integer, dimension(:), allocatable, protected :: &
      NbIncPTCSPrim_ctx

  ! from IncComp to IncPTC (num)
  integer, dimension(:,:,:), allocatable, protected :: &
      NumIncComp2NumIncPTC_ctx

  ! from IncPTC to IncComp (i and alpha)
  integer, dimension(:,:), allocatable, protected :: &
      NumIncPTC2NumIncComp_comp_ctx, &
      NumIncPTC2NumIncComp_phase_ctx

  public :: &
      NumbyContext_make,  &
      NumbyContext_free

  private :: &
      NumbyContext_is_phase_present, &
      NumbyContext_PhaseComp, &
      NumbyContext_Inc,   &
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

    deallocate( NbCompCtilde_ctx)

    deallocate( NbIncPTC_ctx)
    deallocate( NbIncPTCSPrim_ctx)

    deallocate( NumIncComp2NumIncPTC_ctx)
    deallocate( NumIncPTC2NumIncComp_comp_ctx)
    deallocate( NumIncPTC2NumIncComp_phase_ctx)

    deallocate( NumCompCtilde_ctx)
    deallocate( Num2PhasesEqEquilibre_ctx)

    deallocate( NbEqEquilibre_ctx)
    deallocate( NbEqFermeture_ctx)

  end subroutine NumbyContext_Free

  logical(c_bool) function NumbyContext_is_phase_present(configuration, iph, ic) result(here)
    type(ModelConfiguration), intent(in) :: configuration
    integer, intent(in)  :: iph, ic
    integer :: k

    here = .false.
    do k=1, configuration%NbPhasePresente_ctx(ic)
        if(configuration%NumPhasePresente_ctx(k, ic)==iph) then
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

    allocate( NbCompCtilde_ctx(NbContexte))
    allocate( NumCompCtilde_ctx(NbComp, NbContexte))
    NumCompCtilde_ctx(:,:) = 0

    do ic = 1, NbContexte

      n=0
      do icp=1, NbComp ! loop of components

        ! check if icp in C_tidle
        IsCtidle = .true.
        do iph=1, NbPhase

          if ( configuration%MCP(icp,iph)==1 .and. &
          NumbyContext_is_phase_present(configuration,iph,ic)) then
            IsCtidle = .false.
            exit
          end if
        end do

        if(IsCtidle .eqv. .true.) then
          n = n + 1
          NumCompCtilde_ctx(n,ic) = icp
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

    ! 1. Nb d'inconnues P (T) et C fct du contexte
    allocate( NbIncPTC_ctx(NbContexte))

    do ic = 1, NbContexte

      n = 1 ! P

      do i=1, configuration%NbPhasePresente_ctx(ic)
        iph = configuration%NumPhasePresente_ctx(i,ic) ! phase

        ! nb of components in iph, using MCP
        do icp=1, NbComp

          if(configuration%MCP(icp,iph)==1) then ! if component in phase iph
            n = n + 1
          end if
        end do
      enddo

      NbIncPTC_ctx(ic) = n + configuration%IndThermique ! if thermique
    end do

    ! 2. Nb d'inconnues P (T) et C S Prim fct du contexte
    allocate( NbIncPTCSPrim_ctx(NbContexte))

    do ic=1, NbContexte
      NbIncPTCSPrim_ctx(ic) = NbIncPTC_ctx(ic) &
          + configuration%NbPhasePresente_ctx(ic) - NbEqFermeture_ctx(ic) - 1
    end do

    ! 3.1  from IncComp to IncPTC (num)
    allocate( NumIncComp2NumIncPTC_ctx(NbComp, NbPhase, NbContexte))

    ! 3.2  from IncPTC to IncComp (i and alpha)
    allocate( NumIncPTC2NumIncComp_comp_ctx(NbIncPTCMax, NbContexte))
    allocate( NumIncPTC2NumIncComp_phase_ctx(NbIncPTCMax, NbContexte))

    NumIncComp2NumIncPTC_ctx(:,:,:) = 0
    NumIncPTC2NumIncComp_comp_ctx(:,:) = 0
    NumIncPTC2NumIncComp_phase_ctx(:,:) = 0

    do ic=1, NbContexte

      n = 1 + configuration%IndThermique ! P and T

      ! loop of phase
      do i=1,configuration%NbPhasePresente_ctx(ic)
        iph = configuration%NumPhasePresente_ctx(i,ic)

        ! loop of component in phase iph
        do icp=1, NbComp

          if(configuration%MCP(icp,iph)==1) then ! if component in phase iph

            n = n + 1

            NumIncComp2NumIncPTC_ctx(icp, iph, ic) = n
#ifdef DEBUG_LOISTHEMOHYDRO
            ! FIXME: Remove comment
            write(*,*) ic, 'NumIncPTC2NumIncComp_comp_ctx cp=', icp, 'ph=', iph, 'n=', n
#endif
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

    allocate( NumCompEqEquilibre_ctx(configuration%NbEqEquilibreMax, NbContexte))
    allocate( Num2PhasesEqEquilibre_ctx(2, configuration%NbEqEquilibreMax, NbContexte))

    allocate( NbEqEquilibre_ctx(NbContexte))
    allocate( NbEqFermeture_ctx(NbContexte))

    do ic = 1,NbContexte

      n = 0 ! used for NbEquilibre

      do icp=1, NbComp

        nphi = 0
        PhPrComp(:) = 0

        ! PhPrComp: phase present and contains icp
        do iph=1, NbPhase

          if( configuration%MCP(icp,iph)==1 .and. &
            NumbyContext_is_phase_present(configuration,iph,ic)) then
            nphi = nphi + 1
            PhPrComp(nphi) = iph
          end if
        end do

        if(nphi>=2) then

          ! loop of phases in f_i^alpha * C_i^alpha = f_i^beta * C_i^beta, for Q
          do i=1, nphi-1
            NumCompEqEquilibre_ctx(i+n, ic) = icp ! i
            Num2PhasesEqEquilibre_ctx(1, i+n, ic) = PhPrComp(1)   ! alpha
            Num2PhasesEqEquilibre_ctx(2, i+n, ic) = PhPrComp(i+1) ! beta
          end do

          n = n + nphi - 1
        end if

      end do ! end of loop icp

      NbEqEquilibre_ctx(ic) = n
      NbEqFermeture_ctx(ic) = configuration%NbPhasePresente_ctx(ic) + n
      ! NbIncPTCPrim(ic) = NbIncPTC(ic) - NbEqFermeture_ctx(ic)
    end do

  end subroutine NumbyContext_Eq

end module NumbyContext
