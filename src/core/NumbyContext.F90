module NumbyContext

  use DefModel

  implicit none

  ! P_Q: ensemble des phases presentes fct du contexte 
  integer, dimension(:), allocatable, protected :: &
      NbPhasePresente_ctx
  integer, dimension(:,:), allocatable, protected :: &
      NumPhasePresente_ctx

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

  ! Inc (P, T, C_i^alpha): IncPTC
  ! Inc C_i^alpha:         IncComp

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

  ! phase is present or not
  logical, allocatable, dimension(:,:), protected :: &
      IndPhase_ctx

  public :: &
      NumbyContext_make,  &
      NumbyContext_free

  private :: &
      NumbyContext_IndPhaseDeNumContexte, &
      NumbyContext_PhaseComp, &
      NumbyContext_Inc,   &
      NumbyContext_Eq

contains

  ! main subroutine of this module
  subroutine NumbyContext_make

    integer :: i, j, k

    ! IndPhase
    call NumbyContext_IndPhaseDeNumContexte

    ! info phase and comp
    call NumbyContext_PhaseComp

    ! info equation
    call NumbyContext_Eq

    ! Inc
    call NumbyContext_Inc

  end subroutine NumbyContext_make


  ! free everything
  subroutine NumbyContext_Free

    deallocate( NbPhasePresente_ctx)
    deallocate( NumPhasePresente_ctx)

    deallocate( NbCompCtilde_ctx)
    deallocate( NumCompCtilde_ctx)

    deallocate( NbIncPTC_ctx)
    deallocate( NbIncPTCSPrim_ctx)

    deallocate( NumIncComp2NumIncPTC_ctx)
    deallocate( NumIncPTC2NumIncComp_comp_ctx)
    deallocate( NumIncPTC2NumIncComp_phase_ctx)

    deallocate( NumCompCtilde_ctx)
    deallocate( Num2PhasesEqEquilibre_ctx)

    deallocate( NbEqEquilibre_ctx)
    deallocate( NbEqFermeture_ctx)

    deallocate( IndPhase_ctx)

  end subroutine NumbyContext_Free


  subroutine NumbyContext_PhaseComp

    ! 1. Nb/Num de phases presentes fct du contexte 
    ! 2. Ensembles Ctilde fct du contexte 

    integer :: ic, iph, icp, n
    logical :: IsCtidle

    ! 1. Nb/Num de phases presentes fct du contexte

    allocate( NbPhasePresente_ctx(NbContexte))
    allocate( NumPhasePresente_ctx(NbPhase, NbContexte))
    NumPhasePresente_ctx(:,:) = 0

    do ic=1, NbContexte

      ! from IndPh to Nb/NumPhasePresente_ctx
      n = 0
      do iph=1, NbPhase

        if ( IndPhase_ctx(iph,ic) .eqv. .true.) then 
          n = n + 1
          NumPhasePresente_ctx(n,ic) = iph
        endif
      enddo

      NbPhasePresente_ctx(ic) = n
    enddo

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

          if ((MCP(icp,iph)==1) .and. (IndPhase_ctx(iph,ic) .eqv. .true.)) then 
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


  subroutine NumbyContext_Inc

    integer :: i, ic, iph, icp, n

    ! 1. Nb d'inconnues P (T) et C fct du contexte
    allocate( NbIncPTC_ctx(NbContexte))

    do ic = 1, NbContexte 

      n = 1 ! P

      do i=1, NbPhasePresente_ctx(ic)
        iph = NumPhasePresente_ctx(i,ic) ! phase

        ! nb of components in iph, using MCP
        do icp=1, NbComp

          if(MCP(icp,iph)==1) then ! if component in phase iph
            n = n + 1
          end if
        end do
      enddo

      NbIncPTC_ctx(ic) = n + IndThermique ! if thermique
    end do

    ! 2. Nb d'inconnues P (T) et C S Prim fct du contexte
    allocate( NbIncPTCSPrim_ctx(NbContexte))

    do ic=1, NbContexte
      NbIncPTCSPrim_ctx(ic) = NbIncPTC_ctx(ic) &
          + NbPhasePresente_ctx(ic) - NbEqFermeture_ctx(ic) - 1
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

      n = 1 + IndThermique ! P and T

      ! loop of phase
      do i=1,NbPhasePresente_ctx(ic)
        iph = NumPhasePresente_ctx(i,ic)

        ! loop of component in phase iph
        do icp=1, NbComp

          if(MCP(icp,iph)==1) then ! if component in phase iph

            n = n + 1

            NumIncComp2NumIncPTC_ctx(icp, iph, ic) = n
#ifndef NDEBUG              
            ! FIXME: Remove comment
            write(*,*) ic, 'NumIncPTC2NumIncComp_comp_ctx cp=', icp, 'ph=', iph, 'n=', n
#endif !NDEBUG              
            NumIncPTC2NumIncComp_comp_ctx(n, ic) = icp
            NumIncPTC2NumIncComp_phase_ctx(n, ic) = iph
          end if
        enddo

      end do
    enddo

  end subroutine NumbyContext_Inc


  subroutine NumbyContext_Eq

    integer :: ic, iph, icp, nphi, n, i
    integer :: PhPrComp(NbPhase)

    allocate( NumCompEqEquilibre_ctx(NbEqEquilibreMax, NbContexte))
    allocate( Num2PhasesEqEquilibre_ctx(2, NbEqEquilibreMax, NbContexte))

    allocate( NbEqEquilibre_ctx(NbContexte))
    allocate( NbEqFermeture_ctx(NbContexte))

    do ic = 1,NbContexte 

      n = 0 ! used for NbEquilibre

      do icp=1, NbComp

        nphi = 0
        PhPrComp(:) = 0

        ! PhPrComp: phase present and contains icp 
        do iph=1, NbPhase

          if((MCP(icp,iph)==1) .and. (IndPhase_ctx(iph,ic) .eqv. .true.)) then
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
      NbEqFermeture_ctx(ic) = NbPhasePresente_ctx(ic) + n
      ! NbIncPTCPrim(ic) = NbIncPTC(ic) - NbEqFermeture_ctx(ic)
    end do

  end subroutine NumbyContext_Eq


  ! from ic to IndPh(NbPhase)
  ! IndPh(i)=1 -> phase i is presente in this contexte ic
  subroutine NumbyContext_IndPhaseDeNumContexte

    integer :: ic, mic, i

    allocate(IndPhase_ctx(NbPhase, NbContexte))

    do ic = 1, NbContexte

      mic = ic ! 

      ! from decimal to binary 
      do i=NbPhase, 1, -1

        if(mic >= 2**(i-1)) then
          IndPhase_ctx(i,ic) = .true.
          mic = mic - 2**(i-1)
        else
          IndPhase_ctx(i,ic) = .false.
        end if
      end do
    end do

  end subroutine NumbyContext_IndPhaseDeNumContexte

end module NumbyContext

  
    
    ! print*, "MCP ", MCP

    ! do i=1, NbContexte
    !    do j=1, NbPhase
    !       print*, IndPh(j,i)
    !    end do
    !    print*, ""
    ! end do

    ! print*, " "
    ! print*, NbCompCtilde_ctx

    ! do i=1,NbContexte
    !    do j=1,NbComp
    !       print*, NumCompCtilde_ctx(j,i)
    !    end do
    ! end do

    ! print*, NbPhasePresente_ctx
    ! do i=1,NbContexte
    !    do j=1,NbPhase
    !       print*, NumPhasePresente_ctx(j,i)
    !    end do
    ! end do

    ! print*, NbEqEquilibre_ctx
    ! print*, NbEqFermeture_ctx

    ! do i=1, NbContexte
    !    do j=1, NbEqEquilibreMax
    !       print*, NumCompEqEquilibre_ctx(j,i)
    !       print*, Num2PhasesEqEquilibre_ctx(1,j,i)
    !       print*, Num2PhasesEqEquilibre_ctx(2,j,i)
    !    end do
    ! end do

    ! print*, NbIncPTC_ctx
    ! do i=1, NbContexte
    !    do j=1, NbPhase
    !       do k=1, NbComp
    !          print*, NumIncComp2NumIncPTC_ctx(k,j,i)
    !       end do
    !    end do
    ! end do

    ! do i=1, NbContexte
    !    print*, NumIncPTC2NumIncComp_comp_ctx(:,i)
    ! end do

    ! do i=1, NbContexte
    !    print*, NumIncPTC2NumIncComp_phase_ctx(:,i)
    ! end do

    
    ! print*, " "
