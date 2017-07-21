module IncCV

  use MeshSchema
  use DefModel

  use NumbyContext  
  use CommonMPI

#ifdef _HDF5_
  use hdf5
  use h5lt
#endif

  use iso_c_binding

  implicit none

  !> Unknown for Degree Of Freedom (including thermal). DOF can be Cell, Fracture Face or Node.
  TYPE TYPE_IncCV

     integer(c_int) :: ic !< context (???)

     real(c_double) ::         & ! values of Inc
          Pression,              & !< Pressure of the element
          Temperature,           & !< Temperature of the element
          Comp(NbComp, NbPhase), & !< Molar composition of the element
          Saturation(NbPhase),   & !< Saturation of the element
          AccVol(NbCompThermique)  !< ??? of the element

  end TYPE TYPE_IncCV
  
  ! Data of well prod/inj, Ref=max for inj and min for prod
  ! TYPE TYPE_DataWell
  !    character :: IndWell ! indwell: p for pressure mode, f for flowrate mode
  !    double precision :: &
  !         PressionRef, &
  !         Flowrate ! < 0, reference flowrate
  ! end type TYPE_DataWell
  
  !> Type for the perforations, stores informations which are not constant in the well
  TYPE TYPE_PhysPerfoWell
     double precision :: &
          Pression,      & !< Pressure at the perforation
          Temperature,   & !< Temperature at the perforation
          Density,       & !< Density at the perforation: constant per edge, stored at node parent
          PressureDrop     !< Pressure drop at the perforation, used to construct Pressure from the head pressure
          ! FluxMolar(NbComp), & !< Molar flux at the perforation, q_{w,s,i}
          ! FluxEnergy           !< Energy flux at the perforation, q_{w,s,e}
  end TYPE TYPE_PhysPerfoWell

  
  !> to allow = between two TYPE_IncCV
  interface assignment(=)
     module procedure assign_type_inccv
  end interface assignment(=)

  ! Inc for current time step: cell, fracture faces and nodes
  TYPE(TYPE_IncCV), allocatable, dimension(:), target, public :: &
       IncCell, & !< Cell unknowns for current time step
       IncFrac, & !< Fracture Face unknowns for current time step
       IncNode    !< Node unknowns for current time step

  ! well pressure for current time step
  real(c_double), allocatable, dimension(:), target, public :: &
       IncPressionWellInj, & !< Injection Well unknown: head pressure for current time step
       IncPressionWellProd   !< Production Well unknown: head pressure for current time step

  ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
  TYPE(TYPE_PhysPerfoWell), allocatable, dimension(:), target, public :: &  
       PerfoWellInj, & !< Injection Well informations at each perforation for current time step
       PerfoWellProd   !< Production Well informations at each perforation for current time step 

  ! Dir BC
  TYPE(TYPE_IncCV), allocatable, dimension(:), target, public :: &
       IncNodeDirBC !< Dirichlet boundary unknowns for current time step (size NbNodeLocal)

  ! Inc for previous time step: current time step - 1
  TYPE(TYPE_IncCV), allocatable, dimension(:), public :: &
       IncCellPreviousTimeStep, & !< Cell unknowns for previous time step
       IncFracPreviousTimeStep, & !< Fracture Face unknowns for previous time step
       IncNodePreviousTimeStep    !< Node unknowns for previous time step
  ! well pressure for previous time step
  double precision, allocatable, dimension(:), target, public :: &
       IncPressionWellInjPreviousTimeStep, & !< Injection Well unknown: head pressure for previous time step
       IncPressionWellProdPreviousTimeStep   !< Production Well unknown: head pressure for previous time step

  ! Sorted NodebyWellInjLocal%Num and %Val according to z-coordinate for each well
  integer, allocatable, dimension(:), protected :: ZSortedInj_Znum
  double precision, allocatable, dimension(:), protected :: ZSortedInj_Zval

  public :: &
       IncCV_allocate, &
       IncCV_UpdateDirBCValue, &
       IncCV_NewtonRelax,      &
       IncCV_NewtonIncrement,  &
       IncCV_LoadIncPreviousTimeStep, &
       IncCV_SaveIncPreviousTimeStep, &
       IncCV_free

  ! The following subroutines are defined in:
  ! DefInitBCvalues.F90
  public :: &
       IncCV_SetInitialValue,   &
       IncCV_SetDirBCValue

  private :: &
       IncCV_NewtonIncrement_reservoir

contains

  ! IncCV_SetInitialvalue and
  ! IncCV_SetDirBCvalue are defined in: 
#include "DefInitBCvalues.F90"

  !> \brief Define operator = between two TYPE_IncCV:  inc2=inc1
  subroutine assign_type_inccv(inc2, inc1)

    type(TYPE_IncCV), intent(in) :: inc1
    type(TYPE_IncCV), intent(out) :: inc2

    inc2%ic = inc1%ic

    inc2%Pression = inc1%Pression
#ifdef _THERMIQUE_
    inc2%Temperature = inc1%Temperature
#endif
    inc2%Comp(:,:) = inc1%Comp(:,:)
    inc2%Saturation(:) = inc1%Saturation(:)
    inc2%AccVol(:) = inc1%AccVol(:)

  end subroutine assign_type_inccv

  !> \brief Allocate unknown vectors
  subroutine IncCV_allocate

    integer :: Nb, Nnz

    allocate(IncCell(NbCellLocal_Ncpus(commRank+1)))
    allocate(IncFrac(NbFracLocal_Ncpus(commRank+1)))
    allocate(IncNode(NbNodeLocal_Ncpus(commRank+1)))

    allocate(IncPressionWellInj(NbWellInjLocal_Ncpus(commRank+1)))
    allocate(IncPressionWellProd(NbWellProdLocal_Ncpus(commRank+1)))

    Nb = NodebyWellInjLocal%Nb
    Nnz = NodebyWellInjLocal%Pt(Nb+1)
    allocate(PerfoWellInj(Nnz))

    Nb = NodebyWellProdLocal%Nb
    Nnz = NodebyWellProdLocal%Pt(Nb+1)
    allocate(PerfoWellProd(Nnz))

    allocate(IncNodeDirBC(NbNodeLocal_Ncpus(commRank+1)))

    allocate(IncCellPreviousTimeStep(NbCellLocal_Ncpus(commRank+1)))
    allocate(IncFracPreviousTimeStep(NbFracLocal_Ncpus(commRank+1)))
    allocate(IncNodePreviousTimeStep(NbNodeLocal_Ncpus(commRank+1)))
    allocate(IncPressionWellInjPreviousTimeStep(NbWellInjLocal_Ncpus(commRank+1)))
    allocate(IncPressionWellProdPreviousTimeStep(NbWellProdLocal_Ncpus(commRank+1)))

    Nnz = NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)
    allocate( ZSortedInj_Znum(Nnz) )
    allocate( ZSortedInj_Zval(Nnz) )

    ! Sort injection well
    ! and save the results in vector ZsortedInj_Znum and ZsortedInj_Zval
    call IncCV_SortHeightWellInj 

  end subroutine IncCV_allocate


  !> \brief Update the values of the node IncNode(k) if k is Dirichlet
  !! with the Dirichlet boundary values IncNodeDirBC(k)
  subroutine IncCV_UpdateDirBCValue

    integer :: k

    ! Pressure
    do k=1, NbNodeLocal_Ncpus(commRank+1)

       if(IdNodeLocal(k)%P == "d") then
          IncNode(k)%ic = IncNodeDirBC(k)%ic
          IncNode(k)%Pression = IncNodeDirBC(k)%Pression
          IncNode(k)%Saturation(:) = IncNodeDirBC(k)%Saturation(:)
          IncNode(k)%Comp(:,:) = IncNodeDirBC(k)%Comp(:,:)
       end if

#ifdef _THERMIQUE_

       ! Temperature
       if(IdNodeLocal(k)%T == "d") then
          IncNode(k)%ic = IncNodeDirBC(k)%ic
          IncNode(k)%Temperature = IncNodeDirBC(k)%Temperature
       end if
#endif
    end do

  end subroutine IncCV_UpdateDirBCValue


  !> \brief Newton increment of Nodes, Fracture Faces, Cells and Wells unknows.
  subroutine IncCV_NewtonIncrement( &
       NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
       NewtonIncreWellInj, NewtonIncreWellProd, relax)

    double precision, dimension(:,:), intent(in) :: &
         NewtonIncreNode, &
         NewtonIncreFrac, &
         NewtonIncreCell
    double precision, dimension(:), intent(in) :: &
         NewtonIncreWellInj, &
         NewtonIncreWellProd

    double precision, intent(in) :: relax

    integer :: k

    ! nodes
    do k=1, NbNodeLocal_Ncpus(commRank+1)
       call IncCV_NewtonIncrement_reservoir(IncNode(k),NewtonIncreNode(:,k), relax)
    end do

    ! fracture faces
    do k=1, NbFracLocal_Ncpus(commRank+1)
       call IncCV_NewtonIncrement_reservoir(IncFrac(k),NewtonIncreFrac(:,k), relax)
    end do

    ! cells
    do k=1, NbCellLocal_Ncpus(commRank+1)
       call IncCV_NewtonIncrement_reservoir(IncCell(k),NewtonIncreCell(:,k), relax)
    end do

    ! injection wells (head Pressure)
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       IncPressionWellInj(k) = IncPressionWellInj(k) + relax * NewtonIncreWellInj(k)
    end do

    ! production wells
    do k=1, NbWellProdLocal_Ncpus(commRank+1)
       IncPressionWellProd(k) = IncPressionWellProd(k) + relax * NewtonIncreWellProd(k)
    end do

  end subroutine IncCV_NewtonIncrement


  !> \brief Compute relaxation in Newton.
  !!
  !! relax = min(1, IncreObj/NewtonIncreObjMax)                   <br>
  !! where IncreObj is set by the user in DefModel.F90            <br>
  !! and NewtonIncreObjMax is the maximum of the Nemton increment 
  !! in current iteration
  subroutine IncCV_NewtonRelax( &
       NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax)

    double precision, dimension(:,:), intent(in) :: &
         NewtonIncreNode, &
         NewtonIncreFrac, &
         NewtonIncreCell

    double precision, intent(out) :: relax

    double precision ::   &
         incremaxlocal_P, &
         incremaxlocal_T, &
         incremaxlocal_C(NbComp, NbPhase), &
         incremaxlocal_S(NbPhase), &
         relaxlocal

    integer :: k, i, ic, iph, icp, j, Ierr
    integer :: NbIncPTC

    incremaxlocal_P = 0.d0

#ifdef _THERMIQUE_
    incremaxlocal_T = 0.d0
#endif
    incremaxlocal_C(:,:) = 0.d0
    incremaxlocal_S(:) = 0.d0

    ! max Newton increment node
    do k=1, NbNodeOwn_Ncpus(commRank+1)

       if(IdNodeLocal(k)%P /= "d") then

          ic = IncNode(k)%ic
          NbIncPTC = NbIncPTC_ctx(ic)

          incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreNode(1,k)))

          do i=2+IndThermique, NbIncPTC
             icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
             iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

             incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), abs(NewtonIncreNode(i,k)))
          enddo

          do i=1, NbPhasePresente_ctx(ic)
             iph = NumPhasePresente_ctx(i,ic)

             incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreNode(iph+NbIncPTC,k)))
          end do
       end if

#ifdef _THERMIQUE_
       if(IdNodeLocal(k)%T /= "d") then
          incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreNode(2,k)))
       end if
#endif       

    end do

    ! max Newton increment fracture face
    do k=1, NbFracOwn_Ncpus(commRank+1)

       ic = IncFrac(k)%ic
       NbIncPTC = NbIncPTC_ctx(ic)

       incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreFrac(1,k)))

#ifdef _THERMIQUE_
       incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreFrac(2,k)))
#endif       

       do i=2+IndThermique, NbIncPTC
          icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
          iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

          incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), abs(NewtonIncreFrac(i,k)))
       enddo

       do i=1, NbPhasePresente_ctx(ic)
          iph = NumPhasePresente_ctx(i,ic)

          incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreFrac(iph+NbIncPTC,k)))
       end do
    end do

    ! max Newton increment cell
    do k=1, NbCellOwn_Ncpus(commRank+1)

       ic = IncCell(k)%ic
       NbIncPTC = NbIncPTC_ctx(ic)

       incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreCell(1,k)))

#ifdef _THERMIQUE_
       incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreCell(2,k)))
#endif       

       do i=2+IndThermique, NbIncPTC
          icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
          iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

          incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), abs(NewtonIncreCell(i,k)))
       enddo

       do i=1, NbPhasePresente_ctx(ic)
          iph = NumPhasePresente_ctx(i,ic)

          incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreCell(iph+NbIncPTC,k)))
       end do
    end do

    ! relax local = min(1, increobj/incremax)
    relaxlocal = min(1.d0, NewtonIncreObj_P/incremaxlocal_P) ! P

#ifdef _THERMIQUE_
    relaxlocal = min(relaxlocal, NewtonIncreObj_T/incremaxlocal_T) ! T
#endif

    do i=1, NbPhase ! C_i^alpha
       do j=1, NbComp

          if( MCP(j,i)==1 .and. abs(incremaxlocal_C(j,i))>eps) then
             relaxlocal = min(relaxlocal, NewtonIncreObj_C/incremaxlocal_C(j,i))
          end if
       end do
    end do

    do i=1, NbPhase ! S^alpha
       if(abs(incremaxlocal_S(i))>eps) then
          relaxlocal = min(relaxlocal, NewtonIncreObj_S/incremaxlocal_S(i) )
       end if
    end do

    ! relax global
    call MPI_Allreduce(relaxlocal, relax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)


    ! pmaxlocal = increobj_P/incremaxlocal_P
    ! !tmaxlocal = increobj_T/incremaxlocal_T
    ! do i=1, NbPhase
    !    smaxlocal(i) = increobj_S/incremaxlocal_S(i)
    ! end do

    ! call MPI_Allreduce(pmaxlocal, pmax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)
    ! call MPI_Allreduce(tmaxlocal, tmax, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)
    ! call MPI_Allreduce(smaxlocal, smax, NbPhase, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

    ! if(commRank==0) then
    !    write(*,'(F16.5,F16.5,F16.5)',advance="no"), pmax, smax(:)
    ! end if

  end subroutine IncCV_NewtonRelax


  !> \brief Realize Newton increment of each control volume
  subroutine IncCV_NewtonIncrement_reservoir(inc, incre, relax)

    type(TYPE_IncCV), intent(inout) :: inc
    double precision, intent(in) :: incre(NbIncPTCSMax), relax

    integer :: i, icp, iph
    integer :: NbIncPTC, NbIncPTCS
    integer :: ic

    ic = inc%ic
    NbIncPTC = NbIncPTC_ctx(ic)
    NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente_ctx(ic)

    ! increment Pressure
    inc%Pression = inc%Pression + relax * incre(1)

#ifdef _THERMIQUE_

    ! increment Temperature
    inc%Temperature = inc%Temperature + relax * incre(2)
#endif

    ! increment comp
    do i=2+IndThermique, NbIncPTC

       icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
       iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)
       inc%Comp(icp,iph) = inc%Comp(icp,iph) + relax * incre(i)
    enddo

    ! increment saturation
    do i=1, NbPhasePresente_ctx(ic)
       iph = NumPhasePresente_ctx(i,ic)

       inc%Saturation(iph) = inc%Saturation(iph) + relax * incre(iph+NbIncPTC)
    end do

    ! AccVol
    do i=1, NbCompCtilde_ctx(ic) 
       icp = NumCompCtilde_ctx(i,ic)
       inc%AccVol(icp) = incre(NbIncPTCS+i)
    end do

  end subroutine IncCV_NewtonIncrement_reservoir


  !> \brief Compute time step (Delta_t) for the next time iteration
  subroutine IncCV_ComputeTimeStep(Delta_t, TimeCurrent)

    double precision, intent(inout) :: Delta_t
    double precision, intent(in) :: TimeCurrent

    double precision :: Delta_tloc, alpha

    double precision ::   &
         incremaxlocal_P, &
         incremaxlocal_T, &
         incremaxlocal_C(NbComp, NbPhase), &
         incremaxlocal_S(NbPhase)

    integer :: Ierr

    incremaxlocal_P = 0.d0

#ifdef _THERMIQUE_
    incremaxlocal_T = 0.d0
#endif
    incremaxlocal_C(:,:) = 0.d0
    incremaxlocal_S(:) = 0.d0

    ! compute increments of inc of current time step and previsous time step

    !     ! loop of node
    !     do k=1, NbNodeLocal_Ncpus(commRank+1)

    !        if(IdNodeLocal(k)%P /= "d") then

    !           ic = IncNode(k)%ic
    !           NbIncPTC = NbIncPTC_ctx(ic)

    !           incremaxlocal_P = max(incremaxlocal_P, &
    !                abs(IncNode(k)%Pression-IncNodePreviousTimeStep(k)%Pression))

    !           do i=2+IndThermique, NbIncPTC
    !              icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
    !              iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

    !              incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), &
    !                   abs(IncNode(k)%Comp(icp,iph)-IncNodePreviousTimeStep(k)%Comp(icp,iph)))
    !           enddo

    !           do i=1, NbPhasePresente_ctx(ic)
    !              iph = NumPhasePresente_ctx(i,ic)

    !              incremaxlocal_S(iph) = max(incremaxlocal_S(iph), &
    !                   abs(IncNode(k)%Saturation(iph)-IncNodePreviousTimeStep(k)%Saturation(iph)))
    !           end do

    ! #ifdef _THERMIQUE_
    !           if(IdNodeLocal(k)%T /= "d") then
    !              incremaxlocal_T = max(incremaxlocal_T, &
    !                   abs(IncNode(k)%Temperature-IncNodePreviousTimeStep(k)%Temperature))
    !           end if
    ! #endif                 
    !        end if
    !     end do

    !     ! loop of frac
    !     do k=1, NbFracLocal_Ncpus(commRank+1)

    !        ic = IncFrac(k)%ic
    !        NbIncPTC = NbIncPTC_ctx(ic)

    !        incremaxlocal_P = max(incremaxlocal_P, &
    !             abs(IncFrac(k)%Pression-IncFracPreviousTimeStep(k)%Pression))

    !        do i=2+IndThermique, NbIncPTC
    !           icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
    !           iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

    !           incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), &
    !                abs(IncFrac(k)%Comp(icp,iph)-IncFracPreviousTimeStep(k)%Comp(icp,iph)))
    !        enddo

    !        do i=1, NbPhasePresente_ctx(ic)
    !           iph = NumPhasePresente_ctx(i,ic)

    !           incremaxlocal_S(iph) = max(incremaxlocal_S(iph), &
    !                abs(IncFrac(k)%Saturation(iph)-IncFracPreviousTimeStep(k)%Saturation(iph)))
    !        end do

    ! #ifdef _THERMIQUE_
    !        incremaxlocal_T = max(incremaxlocal_T, &
    !             abs(IncFrac(k)%Temperature-IncFracPreviousTimeStep(k)%Temperature))
    ! #endif       
    !     end do

    !     ! loop of cell
    !     do k=1, NbCellLocal_Ncpus(commRank+1)

    !        ic = IncCell(k)%ic
    !        NbIncPTC = NbIncPTC_ctx(ic)

    !        incremaxlocal_P = max(incremaxlocal_P, &
    !             abs(IncCell(k)%Pression-IncCellPreviousTimeStep(k)%Pression))

    !        do i=2+IndThermique, NbIncPTC
    !           icp = NumIncPTC2NumIncComp_comp_ctx(i,ic)
    !           iph = NumIncPTC2NumIncComp_phase_ctx(i,ic)

    !           incremaxlocal_C(icp,iph) = max(incremaxlocal_C(icp,iph), &
    !                abs(IncCell(k)%Comp(icp,iph)-IncCellPreviousTimeStep(k)%Comp(icp,iph)))
    !        enddo

    !        do i=1, NbPhasePresente_ctx(ic)
    !           iph = NumPhasePresente_ctx(i,ic)

    !           incremaxlocal_S(iph) = max(incremaxlocal_S(iph), &
    !                abs(IncCell(k)%Saturation(iph)-IncCellPreviousTimeStep(k)%Saturation(iph)))
    !        end do

    ! #ifdef _THERMIQUE_
    !        incremaxlocal_T = max(incremaxlocal_T, &
    !             abs(IncCell(k)%Temperature-IncCellPreviousTimeStep(k)%Temperature))
    ! #endif       
    !     end do

    ! ! print*, "increment"
    ! ! print*, incremaxlocal_P, incremaxlocal_T, incremaxlocal_S

    ! compute coffient alpha: min of obj/incrementmax
    alpha = 1.2d0
    !     alpha = min(alpha, TimeStepObj_P/incremaxlocal_P) ! P
    ! #ifdef _THERMIQUE_
    !     alpha = min(alpha, TimeStepObj_T/incremaxlocal_T) ! T
    ! #endif

    !     do i=1, NbPhase ! C_i^alpha
    !        do j=1, NbComp

    !           if( MCP(j,i)==1 .and. incremaxlocal_C(j,i)>eps) then
    !              alpha = min(alpha, TimeStepObj_C/incremaxlocal_C(j,i))
    !           end if
    !        end do
    !     end do

    !     do i=1, NbPhase ! S^alpha
    !        if(abs(incremaxlocal_S(i))>eps) then
    !           alpha = min(alpha, TimeStepObj_S/incremaxlocal_S(i) )
    !        end if
    !     end do

    ! time step: Delta_t
    Delta_tloc = Delta_t * alpha

    call MPI_Allreduce(Delta_tloc, Delta_t, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

    ! ! 1cp2ph-ggf, computed with obj
    ! Delta_t = min(Delta_t, TimeStepMax)

    ! ! 2cp2ph-tet6177
    ! if(TimeCurrent < 2000.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! end if

    ! ! 2cp2ph-tet6177 with Pc
    ! if(TimeCurrent < 200.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else if(TimeCurrent < 1800.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax3)
    ! end if

    ! ! BO-tet6177
    ! if(TimeCurrent < 180.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else if(TimeCurrent < 2000.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax3)
    ! end if


    ! ! 2cp2ph-cpgfrac-uniform
    ! if(TimeCurrent < 100.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! end if

    ! ! BO-cpgfrac-uniform
    ! if(TimeCurrent < 300.d0 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! end if

    ! ! 1cp2ph-instable-car100
    ! if(TimeCurrent < 8.2d3 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else if(TimeCurrent < 9.1d3 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! else if(TimeCurrent < 9.65d3 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax3)
    ! else if(TimeCurrent < 1.05d4 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax4)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax5)
    ! end if


    ! ! 1cp2ph-instable-car50
    ! if(TimeCurrent < 1.37d4 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else if(TimeCurrent < 1.42d4 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! else if(TimeCurrent < 1.44d4 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax3)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax4)
    ! end if


    ! ! 1cp2ph-instable-frac-car120
    ! if(TimeCurrent < 2.5d6 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax1)
    ! else if(TimeCurrent < 17.d6 * OneDay) then
    !    Delta_t = min(Delta_t, TimeStepMax2)
    ! else
    !    Delta_t = min(Delta_t, TimeStepMax3)
    ! end if


    ! mesh Simon 5M, one well inj and one well prod
    Delta_t = min(Delta_t, TimeStepMax)
    
    ! if(commRank==0) then
    !    print*, ""
    !    print*, incremaxlocal_P, incremaxlocal_T, incremaxlocal_S, alpha
    ! end if    

  end subroutine IncCV_ComputeTimeStep


  !> \brief Save current status if it is necessary to start again current time iteration.
  !!
  !! Copy IncObj to IncObjPreviousTimeStep
  subroutine IncCV_SaveIncPreviousTimeStep

    integer :: k

    ! save current status
    do k=1, NbNodeLocal_Ncpus(commRank+1)
       IncNodePreviousTimeStep(k) = IncNode(k)
    end do
    do k=1, NbFracLocal_Ncpus(commRank+1)
       IncFracPreviousTimeStep(k) = IncFrac(k)
    end do
    do k=1, NbCellLocal_Ncpus(commRank+1)
       IncCellPreviousTimeStep(k) = IncCell(k)
    end do
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       IncPressionWellInjPreviousTimeStep(k) = IncPressionWellInj(k)
    end do
    do k=1, NbWellProdLocal_Ncpus(commRank+1)
       IncPressionWellProdPreviousTimeStep(k) = IncPressionWellProd(k)
    end do

  end subroutine IncCV_SaveIncPreviousTimeStep


  !> \brief Load previous status to start again current time iteration.
  !!
  !! Copy IncObjPreviousTimeStep to IncObj
  subroutine IncCV_LoadIncPreviousTimeStep

    integer :: k

    do k=1, NbNodeLocal_Ncpus(commRank+1)
       IncNode(k) = IncNodePreviousTimeStep(k)
    end do
    do k=1, NbFracLocal_Ncpus(commRank+1)
       IncFrac(k) = IncFracPreviousTimeStep(k)
    end do
    do k=1, NbCellLocal_Ncpus(commRank+1)
       IncCell(k) = IncCellPreviousTimeStep(k)
    end do
    do k=1, NbWellInjLocal_Ncpus(commRank+1)
       IncPressionWellInj(k) = IncPressionWellInjPreviousTimeStep(k)
    end do
    do k=1, NbWellProdLocal_Ncpus(commRank+1)
       IncPressionWellProd(k) = IncPressionWellProdPreviousTimeStep(k)
    end do

  end subroutine IncCV_LoadIncPreviousTimeStep

  !> \brief Transform IncCV format to a vector, 
  !! necessary for visualization.
  !!
  !! Transform IncCell into output vector datacell 
  !! and IncFrac into output vector datafrac.                                       <br>
  !! The structure of the vector is                                                 <br>
  !!   (Pressure of all cells, Temperature of all cells,                            <br>
  !!        Comp(1) of all cells, ... , Comp(n) of all cells,                       <br>
  !!            Saturation(1) of all cells, ..., Saturation(n) of all cells)
  subroutine IncCV_ToVec( &
       datavisucell, datavisufrac, &
       datavisuwellinj, datavisuwellprod)

    double precision, dimension(:), intent(inout) :: &
         datavisucell, datavisufrac, datavisuwellinj, datavisuwellprod

    integer :: k, i, j, start
    integer :: NbCellOwn, NbFracOwn

    NbCellOwn = NbCellOwn_Ncpus(commRank+1)
    NbFracOwn = NbFracOwn_Ncpus(commRank+1)

    ! cell Pressure
    do k=1, NbCellOwn
       datavisucell(k) = IncCell(k)%Pression
    end do

    ! cell Temperature
#ifdef _THERMIQUE_

    start = NbCellOwn
    do k=1, NbCellOwn
       datavisucell(k+start) = IncCell(k)%Temperature
    end do
#endif

    ! cell Comp
    start = (1 + IndThermique) * NbCellOwn
    do j=1, NbPhase
       do i=1, NbComp

          if(MCP(i,j)==1) then
             do k=1, NbCellOwn
                datavisucell(k+start) = IncCell(k)%Comp(i,j)
             end do
             start = start + NbCellOwn
          end if

       end do
    end do

    ! cell Saturation
    do i=1, NbPhase
       do k=1, NbCellOwn
          datavisucell(k+start) = IncCell(k)%Saturation(i)
       end do
       start = start + NbCellOwn
    end do

    ! frac Pressure
    do k=1, NbFracOwn
       datavisufrac(k) = IncFrac(k)%Pression
    end do

    ! frac Temperature
#ifdef _THERMIQUE_

    start = NbFracOwn
    do k=1, NbFracOwn
       datavisufrac(k+start) = IncFrac(k)%Temperature
    end do
#endif

    ! frac Comp
    start = (1 + IndThermique) * NbFracOwn
    do j=1, NbPhase
       do i=1, NbComp

          if(MCP(i,j)==1) then
             do k=1, NbFracOwn
                datavisufrac(k+start) = IncFrac(k)%Comp(i,j)
             end do
             start = start + NbFracOwn
          end if

       end do
    end do

    ! frac Saturation
    do i=1, NbPhase
       do k=1, NbFracOwn
          datavisufrac(k+start) = IncFrac(k)%Saturation(i)
       end do
       start = start + NbFracOwn
    end do

    ! pressure at well edges, inj
    ! it is equal to the average of its two nodes
    do i=1, sum(NbEdgebyWellInjLocal(1:NbWellInjOwn_Ncpus(commRank+1)))
       datavisuwellinj(i) = 1.d0 ! not implemented
    end do

    ! pressure at well edges, prod
    ! it is equal to the average of its two nodes
    do i=1, sum(NbEdgebyWellProdLocal(1:NbWellProdOwn_Ncpus(commRank+1)))
       datavisuwellprod(i) = 1.d0 ! not implemented
    end do
    
  end subroutine IncCV_ToVec

  ! sort the nodes of wells by z-cordinate from the smallest to the largest
  ! the results are stored in ZSortedInj_Znum (num) and in ZSortedinj_Zval (z-cordinate)
  subroutine IncCV_SortHeightWellInj

    integer :: s, k, j, Nnz, nums
    integer :: tmp_num
    double precision :: tmp_val

    Nnz = NodebyWellInjLocal%Pt(NodebyWellInjLocal%Nb+1)
    do s=1, Nnz
       nums = NodebyWellInjLocal%Num(s)
       ZSortedInj_Znum(s) = s
       ZSortedInj_Zval(s) = XNodeLocal(3, nums)
    end do

    do k=1, NodebyWellInjLocal%Nb ! = NbWellInjLocal(commRank+1)

       do s=1, NodebyWellInjLocal%Pt(k+1)-NodebyWellInjLocal%Pt(k)
          do j=NodebyWellInjLocal%Pt(k)+1, NodebyWellInjLocal%Pt(k+1) - s

             if(ZSortedInj_Zval(j) > ZSortedInj_Zval(j+1) ) then

                tmp_num = ZSortedInj_Znum(j+1)
                tmp_val = ZSortedInj_Zval(j+1)

                ZSortedInj_Znum(j+1) = ZSortedInj_Znum(j)
                ZSortedInj_Zval(j+1) = ZSortedInj_Zval(j)

                ZSortedInj_Znum(j) = tmp_num
                ZSortedInj_Zval(j) = tmp_val
             end if
          end do
       end do
    end do

  end subroutine IncCV_SortHeightWellInj

  ! compute P_{w,s} using Pw (pressure head) and density
  subroutine IncCV_PressureDropWellProd

    integer :: k, s, m, mph, nums, sparent
    double precision :: Pws, zp, zs, Pdrop, Rhotmp
    double precision :: dPf, dTf, dCf(NbComp), dSf(NbPhase)

    do k=1, NbWellProdLocal_Ncpus(commRank+1)

       ! ! Init PressureDrop as zero
       ! do s=NodebyWellProdLocal%Pt(k)+1, NodebyWellProdLocal%Pt(k+1)
       !    PerfoWellProd(s)%PressureDrop = 0.d0
       ! end do

       ! looping from head to queue
       do s=NodebyWellProdLocal%Pt(k+1), NodebyWellProdLocal%Pt(k)+1, -1
          nums = NodebyWellProdLocal%Num(s)

          ! average density
          PerfoWellProd(s)%Density = 0.d0
          do m=1, NbPhasePresente_ctx(IncNode(nums)%ic)
             mph = NumPhasePresente_ctx(m,IncNode(nums)%ic)

             call f_DensiteMolaire(NbPhase, IncNode(nums)%Pression, IncNode(nums)%Temperature, &
                  IncNode(nums)%Comp(:, mph), IncNode(nums)%Saturation, Rhotmp, dPf, dTf, dCf, dSf)
             PerfoWellProd(s)%Density = PerfoWellProd(s)%Density + Rhotmp * IncNode(nums)%Saturation(m)
          end do

          if(s==NodebyWellProdLocal%Pt(k+1)) then ! head node, P = Pw

             Pws = IncPressionWellProd(k) ! P_{w,s} = Pw

             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = 0.d0 

          else ! Pws = P_{w,parent} + \Delta P_{w,parent}

             zs = XNodeLocal(3,nums) ! z-cordinate of node s          
             zp = XNodeLocal(3,NodeDatabyWellProdLocal%Val(s)%Parent) ! z-cordinate of parent of s

             sparent = NodeDatabyWellProdLocal%Val(s)%PtParent ! parent pointer

             Pdrop = PerfoWellProd(sparent)%Density * Gravite * (zp - zs)
             Pws = PerfoWellProd(sparent)%Pression + Pdrop ! Pws

             PerfoWellProd(s)%Pression = Pws
             PerfoWellProd(s)%PressureDrop = PerfoWellProd(sparent)%PressureDrop + Pdrop
          end if

       end do
    end do

  end subroutine IncCV_PressureDropWellProd


  subroutine IncCV_PressureDropWellInj_integrate(sfirst, slast, direction, Pfirst, T, C)

  integer, intent(in) :: sfirst, slast
  ! direction = 1 upwards / -1 downwards
  integer, intent(in) :: direction
  ! T and C are constant !!!
  double precision, intent(in) :: Pfirst, T, C(NbComp)

  integer :: n, s
  ! nb pieces for discrete integration
  ! FIXME: call quad or something similar
  integer, parameter :: Npiece = 100
  double precision :: Ptmp, Stmp(NbPhase), &
      z1, z2, dz, Pdrop, Rhotmp, dPf, dTf, dCf(NbComp), dSf(NbPhase)
#ifndef NDEBUG
    integer :: Ierr, errcode ! used for MPI_Abort
#endif

  Ptmp = Pfirst
  ! Saturation is fixed to liquid  
  Stmp(PHASE_GAS) = 0.d0
  Stmp(PHASE_WATER) = 1.d0
  PerfoWellInj(sfirst)%Pression  = Pfirst
  do s=sfirst, slast-direction, direction
      z1 = ZSortedInj_Zval(s)
      z2 = ZSortedInj_Zval(s+direction)
#ifndef NDEBUG
      if(direction*(z2-z1)<0) then
          write(*,*) 'Nodes are badly sorted.'
          call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
#endif
      dz = (z2-z1) / Npiece
      do n=1, Npiece
          call f_DensiteMolaire(PHASE_WATER, Ptmp, T, C, Stmp, &
              Rhotmp, dPf, dTf, dCf, dSf)
          Ptmp = Ptmp - direction * Gravite * Rhotmp * dz
      end do
      PerfoWellInj(s+direction)%Pression  = Ptmp
  end do

  end subroutine IncCV_PressureDropWellInj_integrate

  ! compute pressure drop of injection well
  ! integration from node head (w) to node (s)
  subroutine IncCV_PressureDropWellInj

    integer :: s, wk, nbwells, bottom, head, headsortedpos
    double precision :: Phead, T, C(NbComp)

    nbwells = NbWellInjLocal_Ncpus(commRank+1)
    do s=NodebyWellInjLocal%Pt(1)+1, NodebyWellInjLocal%Pt(nbwells+1)
        PerfoWellInj(s)%PressureDrop = 0.d0
    end do
    do wk=1, nbwells
        ! Locate head
        bottom = NodebyWellInjLocal%Pt(wk)+1
        head = NodebyWellInjLocal%Pt(wk+1)
        do headsortedpos=bottom, head
            if(ZsortedInj_Znum(headsortedpos)==head) exit
        end do
        Phead = IncPressionWellInj(wk)
        T = DataWellInjLocal(wk)%Temperature
        C = DataWellInjLocal(wk)%CompTotal
        call IncCV_PressureDropWellInj_integrate(headsortedpos,   head,  1, Phead, T, C)
        call IncCV_PressureDropWellInj_integrate(headsortedpos, bottom, -1, Phead, T, C)
        do s=bottom, head
            PerfoWellInj(s)%PressureDrop = Phead - PerfoWellInj(s)%Pression
        end do
    end do

  end subroutine IncCV_PressureDropWellInj


  !> \brief Deallocate unknowns vectors
  subroutine IncCV_free

    deallocate(IncCell)
    deallocate(IncFrac)
    deallocate(IncNode)

    deallocate(IncPressionWellInj)
    deallocate(IncPressionWellProd)

    deallocate(PerfoWellInj)
    deallocate(PerfoWellProd)

    deallocate(IncNodeDirBC)

    deallocate(IncCellPreviousTimeStep)
    deallocate(IncFracPreviousTimeStep)
    deallocate(IncNodePreviousTimeStep)
    deallocate(IncPressionWellInjPreviousTimeStep)
    deallocate(IncPressionWellProdPreviousTimeStep)

    deallocate(ZSortedInj_Znum)
    deallocate(ZSortedInj_Zval)

  end subroutine IncCV_free



#ifdef _HDF5_
  ! hyperslab visual example
  ! +------------+
  ! | xx xx xx xx|
  ! | xx xx xx xx|
  ! | xx xx xx xx|
  ! |            |
  ! | xx xx xx xx|
  ! | xx xx xx xx|
  ! | xx xx xx xx|
  ! |            |
  ! +------------+
  ! start: a starting location for the hyperslab. In the example start is (0,1).
  ! stride: the number of elements to separate each element or block to be selected. In the example stride is (4,3). If the stride parameter is set to NULL, the stride size defaults to 1 in each dimension.
  ! count: the number of elements or blocks to select along each dimension. In the example, count is (2,4).
  ! block: the size of the block selected from the dataspace. In the example, block is (3,2). If the block parameter is set to NULL, the block size defaults to a single element in each dimension, as if the block array was set to all 1s.

  subroutine h5_init_datatype(inctype_id, perfotype_id)
    integer(HID_T), intent(inout) :: inctype_id, perfotype_id

    ! tmp variables used to define new type
    integer(HID_T)  :: arrayC_type_id, arrayS_type_id, arrayAcc_type_id
    integer(SIZE_T) :: sizetypeint, sizetypedouble, sizetypeinc, sizetypeperfo
    integer(SIZE_T) :: sizetypearray(1)
    integer(SIZE_T) :: offsettype

    integer :: Ierr

    ! C_LOC(X) determines the C address of the argument
    sizetypeint = H5OFFSETOF(C_loc(IncCell(1)%ic), C_loc(IncCell(1)%Pression)) ! int
    sizetypedouble = H5OFFSETOF(C_loc(IncCell(1)%Pression), C_loc(IncCell(1)%Temperature)) ! double

    ! size of new data type
    ! pressure: 1 double, temperature: 1 double, + all fields defined in Type_IncCV
    sizetypeinc = sizetypeint + sizetypedouble * (2 + NbComp*NbPhase + NbPhase + NbCompThermique)

    sizetypeperfo = 3 * sizetypedouble

    ! create and insert
    call H5Tcreate_f(H5T_COMPOUND_F, sizetypeinc, inctype_id, Ierr)

    offsettype = 0
    call H5Tinsert_f(inctype_id, "ic", offsettype, H5T_NATIVE_INTEGER, Ierr)

    offsettype = offsettype + sizetypeint 
    call H5Tinsert_f(inctype_id, "Pression", offsettype, H5T_NATIVE_DOUBLE, Ierr)

    offsettype = offsettype + sizetypedouble
    call H5Tinsert_f(inctype_id, "Temperature", offsettype, H5T_NATIVE_DOUBLE, Ierr)

    offsettype = offsettype + sizetypedouble
    sizetypearray = NbComp * NbPhase
    call h5Tarray_create_f(H5T_NATIVE_DOUBLE, 1, sizetypearray, arrayC_type_id, Ierr)
    call H5Tinsert_f(inctype_id, "Comp", offsettype, arrayC_type_id, Ierr)

    offsettype = offsettype + sizetypedouble * NbPhase*NbComp
    sizetypearray = NbPhase
    call h5Tarray_create_f(H5T_NATIVE_DOUBLE, 1, sizetypearray, arrayS_type_id, Ierr)
    call H5Tinsert_f(inctype_id, "Saturation", offsettype, arrayS_type_id, Ierr)

    offsettype = offsettype + sizetypedouble * NbPhase
    sizetypearray = NbCompThermique 
    call h5Tarray_create_f(H5T_NATIVE_DOUBLE, 1, sizetypearray, arrayAcc_type_id, Ierr)
    call H5Tinsert_f(inctype_id, "AccVol", offsettype, arrayAcc_type_id, Ierr)


    ! TYPE_Perfo
    call H5Tcreate_f(H5T_COMPOUND_F, sizetypeperfo, perfotype_id, Ierr)

    offsettype = 0
    call H5Tinsert_f(perfotype_id, "Pression", offsettype, H5T_NATIVE_DOUBLE, Ierr)

    offsettype = offsettype + sizetypedouble 
    call H5Tinsert_f(perfotype_id, "Temperature", offsettype, H5T_NATIVE_DOUBLE, Ierr)

    offsettype = offsettype + sizetypedouble
    call H5Tinsert_f(perfotype_id, "Density", offsettype, H5T_NATIVE_DOUBLE, Ierr)

  end subroutine h5_init_datatype

  subroutine h5_create_and_write_dset(filename, fieldname, local_Ncpus, type_id, targetptr, file_id)
    type(C_PTR), intent(in) :: targetptr
    character(*), intent(in) :: fieldname, filename
    integer, allocatable, dimension(:), intent(in) :: local_Ncpus
    integer(HID_T), intent(in) :: type_id       ! Compound datatype identifier
    integer(HID_T), intent(inout) :: file_id

    integer(HID_T) :: filespace     ! Filespace identifier in file
    integer(HID_T) :: memspace      ! Dataspace identifier in memory
    integer(HID_T) :: plist_id
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HSIZE_T),  dimension(1)  :: sizeinc    ! Dataset dimensions global.
    integer(HSIZE_T),  dimension(1)  :: sizeincloc ! Dataset dimensions local.
    integer(HSSIZE_T), dimension(1)  :: offset 

    integer :: i, Ierr

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr) ! Stores MPI IO communicator information to the file access property list.

    ! Create the file collectively.
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, Ierr, access_prp = plist_id)
    call h5pclose_f(plist_id, Ierr)

    ! sizeinc, sizeincloc
    sizeinc = 0
    do i=1, commSize
       sizeinc = sizeinc + local_Ncpus(i)
    end do
    sizeincloc = local_Ncpus(commRank+1)

    ! offset
    offset = 0
    do i=1, commRank
       offset = offset + local_Ncpus(i)
    end do

    ! Create the data space for the dataset: filespace 
    call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! Create the dataset with default properties.
    call h5dcreate_f(file_id, trim(fieldname), type_id, filespace, dset_id, Ierr)
    call h5sclose_f(filespace, Ierr)

    ! Select hyperslab in the file: filespace
    call h5dget_space_f(dset_id, filespace, Ierr)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr) ! here, we only specify start and count, not stride nor block (which then defaults to 1)

    ! memspace, Creates a new simple dataspace and opens it for access
    call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! Create property list for collective dataset read: plist_id
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! Write the dataset collectively.
    call h5dwrite_f(dset_id, type_id, targetptr, Ierr, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    call h5dclose_f(dset_id, Ierr)

  end subroutine h5_create_and_write_dset

  subroutine h5_open_and_read_dset(filename, fieldname, local_Ncpus, type_id, targetptr, TimeIter, file_id)
    type(C_PTR), intent(inout) :: targetptr
    character(*), intent(in) :: fieldname, filename
    integer, intent(in) :: TimeIter
    integer, allocatable, dimension(:), intent(in) :: local_Ncpus
    integer(HID_T), intent(in) :: type_id       ! Compound datatype identifier
    integer(HID_T), intent(inout) :: file_id

    integer(HID_T) :: filespace     ! Filespace identifier in file
    integer(HID_T) :: memspace      ! Dataspace identifier in memory
    integer(HID_T) :: plist_id
    integer(HID_T) :: dset_id       ! Dataset identifier

    integer(HSIZE_T),  dimension(1)  :: sizeinc    ! Dataset dimensions global.
    integer(HSIZE_T),  dimension(1)  :: sizeincloc ! Dataset dimensions local.
    integer(HSSIZE_T), dimension(1)  :: offset

    integer :: nbdims
    integer(HSIZE_T), dimension(1) :: numdims, numdimsmax

    integer :: i, Ierr

    ! Setup file access property list with parallel I/O access.
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! Create the file collectively.
    call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, Ierr, access_prp = plist_id)

    call h5pclose_f(plist_id, Ierr)

    ! open dataset
    call h5dopen_f(file_id, fieldname, dset_id, Ierr)

    ! sizeinc, sizeincloc
    sizeinc = 0
    do i=1, commSize
       sizeinc = sizeinc + local_Ncpus(i)
    end do
    sizeincloc = local_Ncpus(commRank+1)

    ! offset
    offset = 0
    do i=1, commRank
       offset = offset + local_Ncpus(i)
    end do

    ! get filespace from dataset
    call h5dget_space_f(dset_id, filespace, Ierr)
    call h5sget_simple_extent_ndims_f(filespace, nbdims, Ierr)
    call h5sget_simple_extent_dims_f(filespace, numdims, numdimsmax, Ierr)

    ! check if size is same
    if((nbdims /= 1) .or. (numdims(1) /= sizeinc(1))) then
       write(*,*) "Error in reading cell solution from file at time step ", TimeIter
       write(*,*) "  sizes of not same "
    end if

    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! memspace
    call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    call h5dread_f(dset_id, type_id, targetptr, Ierr, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    call h5dclose_f(dset_id, Ierr)

  end subroutine h5_open_and_read_dset

  !> \brief Write solution (IncCell, IncNode and IncFrac) to file (hdf5).
  !!
  !! \param[in]    dirname             Name of the output directory
  !! \param[in]    TimeIter            Time step, current time iteration
  !! \param[in]    TimeCurrent         Current time (in days)
  !! \param[in]    Delta_t             Current time discretization
  !! \param[in]    TimeOutput          ???
  !! \param[in]    NewtonNiterTotal    Total number of Newton iterations (from the beginning of the simulation)
  !! \param[in]    NewtonNbFailure     Total number of Newton failures (from the beginning of the simulation)
  !! \param[in]    KspNiterTotal       Total number of Ksp iterations (from the beginning of the simulation)
  !! \param[in]    KspNbFailure        Total number of Ksp failures (from the beginning of the simulation)
  !! \param[in]    comptime_total      Total Computation time (from the beginning of the simulation)
  !! \param[in]    comptime_timestep   Computation time of this time step
  !!
  ! TODO write visu infos to hdf5,
  !       the visu starts always from timestep 0 actuellement
  subroutine IncCV_WriteSolToFile(dirname, &
       TimeIter, TimeCurrent, Delta_t, TimeOutput, &
       NewtonNiterTotal, NewtonNbFailure, KspNiterTotal, KspNbFailure, &
       comptime_total, comptime_timestep)

    character(*), intent(in) :: dirname

    integer, intent(in) :: TimeIter

    double precision, intent(in) :: &
         TimeCurrent, Delta_t, TimeOutput, &
         comptime_total, comptime_timestep

    integer, intent(in) :: &
         NewtonNiterTotal, NewtonNbFailure, &
         KspNiterTotal, KspNbFailure

    character(len=300) :: filename  ! filename

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: plist_id      ! Property list identifier
    integer(HID_T) :: inctype_id, perfotype_id    ! Compound datatype identifier

    ! tmp varibale used to attributes
    integer(HSIZE_T), parameter :: nbattr = 1
    integer, dimension(1) :: attrint
    double precision, dimension(1) :: attrdble

    integer :: Ierr

    ! Initialize FORTRAN interface
    call h5open_f(Ierr)

    ! 1. new data type
    call h5_init_datatype(inctype_id, perfotype_id)

    ! 2. write to file: cell, node, frac, well
    ! 2.1 cell
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_cell.h5"
    call h5_create_and_write_dset(filename, "IncCell", NbCellLocal_Ncpus, inctype_id, C_loc(IncCell(1)), file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr) ! Stores MPI IO communicator information to the file access property list.

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_cell.h5"
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbCellLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbCellLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbCellLocal_Ncpus(i)
    ! end do

    ! ! Create the data space for the dataset: filespace 
    ! call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! ! Create the dataset with default properties.
    ! call h5dcreate_f(file_id, "IncCell", inctype_id, filespace, &
    !      dset_id, Ierr)
    ! call h5sclose_f(filespace, Ierr)

    ! ! Select hyperslab in the file: filespace
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr) ! here, we only specify start and count, not stride nor block (which then defaults to 1)

    ! ! memspace, Creates a new simple dataspace and opens it for access
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! ! Create property list for collective dataset read: plist_id
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! ! Write the dataset collectively.
    ! targetptr = C_LOC(IncCell(1))
    ! call h5dwrite_f(dset_id,  inctype_id, targetptr, Ierr, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! attributes
    attrdble(1) = TimeCurrent
    call H5LTset_attribute_double_f(file_id, "IncCell", "TimeCurrent", attrdble, nbattr, Ierr)
    attrdble(1) = Delta_t
    call H5LTset_attribute_double_f(file_id, "IncCell", "Delta_t", attrdble, nbattr, Ierr)
    attrdble(1) = TimeOutput
    call H5LTset_attribute_double_f(file_id, "IncCell", "TimeOutput", attrdble, nbattr, Ierr)

    attrint(1) = NewtonNiterTotal
    call H5LTset_attribute_int_f(file_id, "IncCell", "NewtonNitertotal", attrint, nbattr, Ierr)
    attrint(1) = NewtonNbFailure
    call H5LTset_attribute_int_f(file_id, "IncCell", "NewtonNbFailure", attrint, nbattr, Ierr)
    attrint(1) = KspNiterTotal
    call H5LTset_attribute_int_f(file_id, "IncCell", "KspNiterTotal", attrint, nbattr, Ierr)
    attrint(1) = KspNbFailure
    call H5LTset_attribute_int_f(file_id, "IncCell", "KspNbFailure", attrint, nbattr, Ierr)

    attrdble(1) = comptime_total
    call H5LTset_attribute_double_f(file_id, "IncCell", "comptime_total", attrdble, nbattr, Ierr)
    attrdble(1) = comptime_timestep
    call H5LTset_attribute_double_f(file_id, "IncCell", "comptime_timestep", attrdble, nbattr, Ierr)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.2 node
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_node.h5"
    call h5_create_and_write_dset(filename, "IncNode", NbNodeLocal_Ncpus, inctype_id, C_loc(IncNode(1)), file_id)
    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_node.h5"
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbNodeLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbNodeLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbNodeLocal_Ncpus(i)
    ! end do

    ! ! Create the data space for the dataset: filespace 
    ! call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! ! Create the dataset with default properties.
    ! call h5dcreate_f(file_id, "IncNode", inctype_id, filespace, &
    !      dset_id, Ierr)
    ! call h5sclose_f(filespace, Ierr)

    ! ! Select hyperslab in the file: filespace
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! ! Create property list for collective dataset read: plist_id
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! ! Write the dataset collectively.
    ! targetptr = C_LOC(IncNode(1))
    ! call h5dwrite_f(dset_id,  inctype_id, targetptr, Ierr, &
    !      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.3 frac
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_frac.h5"
    call h5_create_and_write_dset(filename, "IncFrac", NbFracLocal_Ncpus, inctype_id, C_loc(IncFrac(1)), file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_frac.h5"
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbFracLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbFracLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbFracLocal_Ncpus(i)
    ! end do

    ! ! Create the data space for the dataset: filespace 
    ! call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! ! Create the dataset with default properties.
    ! call h5dcreate_f(file_id, "IncFrac", inctype_id, filespace, &
    !      dset_id, Ierr)
    ! call h5sclose_f(filespace, Ierr)

    ! ! Select hyperslab in the file: filespace
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! ! Create property list for collective dataset read: plist_id
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! ! Write the dataset collectively.
    ! targetptr = C_LOC(IncFrac(1))
    ! call h5dwrite_f(dset_id,  inctype_id, targetptr, Ierr, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.4 well inj
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_well_inj.h5"
    call h5_create_and_write_dset(filename, "PerfoWellInj", NbWellInjLocal_Ncpus, perfotype_id, C_loc(PerfoWellInj(1)), file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_well_inj.h5"
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbWellInjLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbWellInjLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbWellInjLocal_Ncpus(i)
    ! end do

    ! ! Create the data space for the dataset: filespace 
    ! call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! ! Create the dataset with default properties.
    ! call h5dcreate_f(file_id, "PerfoWellInj", perfotype_id, filespace, dset_id, Ierr)
    ! call h5sclose_f(filespace, Ierr)

    ! ! Select hyperslab in the file: filespace
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! ! Create property list for collective dataset read: plist_id
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! ! Write the dataset collectively.
    ! targetptr = C_LOC(PerfoWellInj(1))
    ! call h5dwrite_f(dset_id,  inctype_id, targetptr, Ierr, &
    !      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    ! call h5fclose_f(file_id, Ierr)

    ! ! 2.5 well prod

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_well_prod.h5"
    ! call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbWellProdLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbWellProdLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbWellProdLocal_Ncpus(i)
    ! end do

    ! ! Create the data space for the dataset: filespace 
    ! call h5screate_simple_f(1, sizeinc, filespace, Ierr)

    ! ! Create the dataset with default properties.
    ! call h5dcreate_f(file_id, "PerfoWellInj", perfotype_id, filespace, dset_id, Ierr)
    ! call h5sclose_f(filespace, Ierr)

    ! ! Select hyperslab in the file: filespace
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! ! Create property list for collective dataset read: plist_id
    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! ! Write the dataset collectively.
    ! targetptr = C_LOC(PerfoWellInj(1))
    ! call h5dwrite_f(dset_id,  inctype_id, targetptr, Ierr, file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! close 
    call H5Tclose_f(inctype_id, Ierr)
    call H5Tclose_f(perfotype_id, Ierr)
    call h5close_f(Ierr)

  end subroutine IncCV_WriteSolToFile


  !> \brief Read solution (IncCell, IncNode and IncFrac) from file (hdf5).
  !!
  !! \param[in]    dirname             Name of the output directory
  !! \param[in]    TimeIter            Time step, current time iteration
  !! \param[out]   TimeCurrent         Current time (in days)
  !! \param[out]   Delta_t             Current time discretization
  !! \param[out]   TimeOutput          ???
  !! \param[out]   NewtonNiterTotal    Total number of Newton iterations (from the beginning of the simulation)
  !! \param[out]   NewtonNbFailure     Total number of Newton failures (from the beginning of the simulation)
  !! \param[out]   KspNiterTotal       Total number of Ksp iterations (from the beginning of the simulation)
  !! \param[out]   KspNbFailure        Total number of Ksp failures (from the beginning of the simulation)
  !! \param[out]   comptime_total      Total Computation time (from the beginning of the simulation)
  !! \param[out]   comptime_timestep   Computation time of this time step
  subroutine IncCV_ReadSolFromFile(dirname,           &
       TimeIter, TimeCurrent, Delta_t, TimeOutput,    &
       NewtonNiterTotal, NewtonNbFailure,             &
       KspNiterTotal, KspNbFailure,                   &
       comptime_total, comptime_timestep)

    character(*), intent(in) :: dirname
    integer, intent(in) :: TimeIter

    double precision, intent(out) :: &
         TimeCurrent, Delta_t, TimeOutput, &
         comptime_total, comptime_timestep

    integer, intent(out) :: &
         NewtonNiterTotal, NewtonNbFailure, &
         KspNiterTotal, KspNbFailure

    character(len=300) :: filename  ! filename

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: plist_id      ! Property list identifier
    integer(HID_T) :: inctype_id, perfotype_id    ! Compound datatype identifier

    ! tmp varibale used to attributes
    integer, dimension(1) :: attrint
    double precision, dimension(1) :: attrdble

    type(C_PTR) :: targetptr ! buf in h5dread_f in intent inout ...

    integer :: Ierr

    ! Initialize FORTRAN interface
    call h5open_f(Ierr)

    call h5_init_datatype(inctype_id, perfotype_id)

    ! 2. read from file: cell, node, frac, well

    ! 2.1 cell
    targetptr = C_loc(IncCell(1))
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_cell.h5"
    call h5_open_and_read_dset(filename, "IncCell", NbCellLocal_Ncpus, inctype_id, targetptr, TimeIter, file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_cell.h5"
    ! call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! open dataset
    ! call h5dopen_f(file_id, "IncCell", dset_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbCellLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbCellLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbCellLocal_Ncpus(i)
    ! end do

    ! ! get filespace from dataset
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sget_simple_extent_ndims_f(filespace, nbdims, Ierr)
    ! call h5sget_simple_extent_dims_f(filespace, numdims, numdimsmax, Ierr)

    ! ! check if size is same
    ! if((nbdims/=1) .or. (numdims(1)/=sizeinc(1))) then
    !    write(*,*) "Error in reading cell solution from file at time step ", TimeIter
    !    write(*,*) "  sizes of not same "
    ! end if

    ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! targetptr = C_LOC(IncCell(1))
    ! CALL h5dread_f(dset_id, inctype_id, targetptr, Ierr, &
    !      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! attributes
    call H5LTget_attribute_double_f(file_id, "IncCell", "TimeCurrent", attrdble, Ierr)
    TimeCurrent = attrdble(1)
    call H5LTget_attribute_double_f(file_id, "IncCell", "Delta_t", attrdble, Ierr)
    Delta_t = attrdble(1)
    call H5LTget_attribute_double_f(file_id, "IncCell", "TimeOutput", attrdble, Ierr)
    TimeOutput = attrdble(1)

    call H5LTget_attribute_int_f(file_id, "IncCell", "NewtonNitertotal", attrint, Ierr)
    NewtonNiterTotal = attrint(1)
    call H5LTget_attribute_int_f(file_id, "IncCell", "NewtonNbFailure", attrint, Ierr)
    NewtonNbFailure = attrint(1)
    call H5LTget_attribute_int_f(file_id, "IncCell", "KspNiterTotal", attrint, Ierr)
    KspNiterTotal = attrint(1)
    call H5LTget_attribute_int_f(file_id, "IncCell", "KspNbFailure", attrint, Ierr)
    KspNbFailure = attrint(1)

    call H5LTget_attribute_double_f(file_id, "IncCell", "comptime_total", attrdble, Ierr)
    comptime_total = attrdble(1)
    call H5LTget_attribute_double_f(file_id, "IncCell", "comptime_timestep", attrdble, Ierr)
    comptime_timestep = attrdble(1)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.2 node
    targetptr = C_loc(IncNode(1))
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_node.h5"
    call h5_open_and_read_dset(filename, "IncNode", NbNodeLocal_Ncpus, inctype_id, targetptr, TimeIter, file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_node.h5"
    ! call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! open dataset
    ! call h5dopen_f(file_id, "IncNode", dset_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbNodeLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbNodeLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbNodeLocal_Ncpus(i)
    ! end do

    ! ! get filespace from dataset
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sget_simple_extent_ndims_f(filespace, nbdims, Ierr)
    ! call h5sget_simple_extent_dims_f(filespace, numdims, numdimsmax, Ierr)

    ! ! check if size is same
    ! if((nbdims/=1) .or. (numdims(1)/=sizeinc(1))) then
    !    write(*,*) "Error in reading node solution from file at time step ", TimeIter
    !    write(*,*) "  sizes of not same "
    ! end if

    ! call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! targetptr = C_LOC(IncNode(1))
    ! CALL h5dread_f(dset_id, inctype_id, targetptr, Ierr, &
    !      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.3 frac
    targetptr = C_loc(IncFrac(1))
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_frac.h5"
    call h5_open_and_read_dset(filename, "IncFrac", NbFracLocal_Ncpus, inctype_id, targetptr, TimeIter, file_id)

    ! ! Setup file access property list with parallel I/O access.
    ! call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, Ierr)
    ! call h5pset_fapl_mpio_f(plist_id, ComPASS_COMM_WORLD, MPI_INFO_NULL, Ierr)

    ! ! Create the file collectively.
    ! write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_frac.h5"
    ! call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, &
    !      file_id, Ierr, access_prp = plist_id)

    ! call h5pclose_f(plist_id, Ierr)

    ! ! open dataset
    ! call h5dopen_f(file_id, "IncFrac", dset_id, Ierr)

    ! ! sizeinc, sizeincloc
    ! sizeinc = 0
    ! do i=1, commSize
    !    sizeinc = sizeinc + NbFracLocal_Ncpus(i)
    ! end do
    ! sizeincloc = NbFracLocal_Ncpus(commRank+1)

    ! ! offset
    ! offset = 0
    ! do i=1, commRank
    !    offset = offset + NbFracLocal_Ncpus(i)
    ! end do

    ! ! get filespace from dataset
    ! call h5dget_space_f(dset_id, filespace, Ierr)
    ! call h5sget_simple_extent_ndims_f(filespace, nbdims, Ierr)
    ! call h5sget_simple_extent_dims_f(filespace, numdims, numdimsmax, Ierr)

    ! ! check if size is same
    ! if((nbdims/=1) .or. (numdims(1)/=sizeinc(1))) then
    !    write(*,*) "Error in reading frac solution from file at time step ", TimeIter
    !    write(*,*) "  sizes of not same "
    ! end if

    ! call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, sizeincloc, Ierr)

    ! ! memspace
    ! call h5screate_simple_f(1, sizeincloc, memspace, Ierr) 

    ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, Ierr) 
    ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, Ierr)

    ! targetptr = C_LOC(IncFrac(1))
    ! CALL h5dread_f(dset_id, inctype_id, targetptr, Ierr, &
    !      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! call h5dclose_f(dset_id, Ierr)
    call h5fclose_f(file_id, Ierr)

    ! 2.4 well inj
    targetptr = C_loc(PerfoWellInj(1))
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_inj.h5"
    call h5_open_and_read_dset(filename, "PerfoWellInj", NbWellInjLocal_Ncpus, & 
         perfotype_id, targetptr, TimeIter, file_id)

    ! 2.5 well prod
    targetptr = C_loc(PerfoWellProd(1))
    write(filename, '(A,A,I0,A)') trim(dirname), "/TimeStep_", TimeIter, "_prod.h5"
    call h5_open_and_read_dset(filename, "PerfoWellProd", NbWellProdLocal_Ncpus, &
         perfotype_id, targetptr, TimeIter, file_id)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! close 
    call H5Tclose_f(inctype_id, Ierr)
    call H5Tclose_f(perfotype_id, Ierr)
    call h5close_f(Ierr)

  end subroutine IncCV_ReadSolFromFile

#endif

end module IncCV
