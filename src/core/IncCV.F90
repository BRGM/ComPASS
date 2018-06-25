    !
    ! This file is part of ComPASS.
    !
    ! ComPASS is free software: you can redistribute it and/or modify it under both the terms
    ! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
    ! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
    !

    module IncCV

    use IncCVReservoir
    use IncCVWells

    use iso_c_binding

    implicit none


    public :: &
        IncCV_allocate, &
        IncCV_NewtonIncrement, &
        IncCV_LoadIncPreviousTimeStep, &
        IncCV_SaveIncPreviousTimeStep, &
        IncCV_free

    contains

    subroutine IncCV_allocate

    call IncCVReservoir_allocate
    call IncCVWells_allocate
    
    end subroutine IncCV_allocate

    subroutine IncCV_free

    call IncCVReservoir_free
    call IncCVWells_free

    end subroutine IncCV_free

    !> \brief Newton increment of Nodes, Fracture Faces, Cells and Wells unknows.
    subroutine IncCV_NewtonIncrement( &
        NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
        NewtonIncreWellInj, NewtonIncreWellProd, relax) &
        bind(C, name="IncCV_NewtonIncrement")

    real(c_double), dimension(:, :), intent(in) :: &
        NewtonIncreNode, &
        NewtonIncreFrac, &
        NewtonIncreCell
    real(c_double), dimension(:), intent(in) :: &
        NewtonIncreWellInj, &
        NewtonIncreWellProd

    real(c_double), intent(in), value :: relax

    call IncCVReservoir_NewtonIncrement(NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax)
    call IncCVWells_NewtonIncrement(NewtonIncreWellInj, NewtonIncreWellProd, relax)

    end subroutine IncCV_NewtonIncrement

    !> \brief Compute relaxation in Newton.
    !!
    !! relax = min(1, IncreObj/NewtonIncreObjMax)                   <br>
    !! where IncreObj is set by the user in DefModel.F90            <br>
    !! and NewtonIncreObjMax is the maximum of the Nemton increment
    !! in current iteration
    subroutine IncCV_NewtonRelax( &
        NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax)

    real(c_double), dimension(:, :), intent(in) :: &
        NewtonIncreNode, &
        NewtonIncreFrac, &
        NewtonIncreCell

    real(c_double), intent(out) :: relax

    real(c_double) :: &
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
    incremaxlocal_C(:, :) = 0.d0
    incremaxlocal_S(:) = 0.d0

    ! max Newton increment node
    do k = 1, NbNodeOwn_Ncpus(commRank + 1)

        if (IdNodeLocal(k)%P /= "d") then

            ic = IncNode(k)%ic
            NbIncPTC = NbIncPTC_ctx(ic)

            incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreNode(1, k)))

            do i = 2 + IndThermique, NbIncPTC
                icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
                iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

                incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(NewtonIncreNode(i, k)))
            enddo

            do i = 1, NbPhasePresente_ctx(ic)
                iph = NumPhasePresente_ctx(i, ic)

                incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreNode(iph + NbIncPTC, k)))
            end do
        end if

#ifdef _THERMIQUE_
        if (IdNodeLocal(k)%T /= "d") then
            incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreNode(2, k)))
        end if
#endif

    end do

    ! max Newton increment fracture face
    do k = 1, NbFracOwn_Ncpus(commRank + 1)

        ic = IncFrac(k)%ic
        NbIncPTC = NbIncPTC_ctx(ic)

        incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreFrac(1, k)))

#ifdef _THERMIQUE_
        incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreFrac(2, k)))
#endif

        do i = 2 + IndThermique, NbIncPTC
            icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
            iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

            incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(NewtonIncreFrac(i, k)))
        enddo

        do i = 1, NbPhasePresente_ctx(ic)
            iph = NumPhasePresente_ctx(i, ic)

            incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreFrac(iph + NbIncPTC, k)))
        end do
    end do

    ! max Newton increment cell
    do k = 1, NbCellOwn_Ncpus(commRank + 1)

        ic = IncCell(k)%ic
        NbIncPTC = NbIncPTC_ctx(ic)

        incremaxlocal_P = max(incremaxlocal_P, abs(NewtonIncreCell(1, k)))

#ifdef _THERMIQUE_
        incremaxlocal_T = max(incremaxlocal_T, abs(NewtonIncreCell(2, k)))
#endif

        do i = 2 + IndThermique, NbIncPTC
            icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
            iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)

            incremaxlocal_C(icp, iph) = max(incremaxlocal_C(icp, iph), abs(NewtonIncreCell(i, k)))
        enddo

        do i = 1, NbPhasePresente_ctx(ic)
            iph = NumPhasePresente_ctx(i, ic)

            incremaxlocal_S(iph) = max(incremaxlocal_S(iph), abs(NewtonIncreCell(iph + NbIncPTC, k)))
        end do
    end do

    ! relax local = min(1, increobj/incremax)
    relaxlocal = min(1.d0, NewtonIncreObj_P/incremaxlocal_P) ! P

#ifdef _THERMIQUE_
    relaxlocal = min(relaxlocal, NewtonIncreObj_T/incremaxlocal_T) ! T
#endif

    do i = 1, NbPhase ! C_i^alpha
        do j = 1, NbComp

            if (MCP(j, i) == 1 .and. abs(incremaxlocal_C(j, i)) > eps) then
                relaxlocal = min(relaxlocal, NewtonIncreObj_C/incremaxlocal_C(j, i))
            end if
        end do
    end do

    do i = 1, NbPhase ! S^alpha
        if (abs(incremaxlocal_S(i)) > eps) then
            relaxlocal = min(relaxlocal, NewtonIncreObj_S/incremaxlocal_S(i))
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

!    !> \brief Realize Newton increment of each control volume
!    subroutine IncCV_NewtonIncrement_reservoir(inc, incre, relax)
!
!    type(TYPE_IncCV), intent(inout) :: inc
!    double precision, intent(in) :: incre(NbIncPTCSMax), relax
!
!    integer :: i, icp, iph
!    integer :: NbIncPTC, NbIncPTCS
!    integer :: ic
!
!    ic = inc%ic
!    NbIncPTC = NbIncPTC_ctx(ic)
!    NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente_ctx(ic)
!
!    ! increment Pressure
!    inc%Pression = inc%Pression + relax*incre(1)
!
!    !    write(*,*)' increment P ',relax,incre(1)
!
!#ifdef _THERMIQUE_
!
!    ! increment Temperature
!    inc%Temperature = inc%Temperature + relax*incre(2)
!#endif
!
!    ! increment comp
!    do i = 2 + IndThermique, NbIncPTC
!
!        icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
!        iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)
!        inc%Comp(icp, iph) = inc%Comp(icp, iph) + relax*incre(i)
!    enddo
!
!    ! increment saturation
!    do i = 1, NbPhasePresente_ctx(ic)
!        iph = NumPhasePresente_ctx(i, ic)
!
!        inc%Saturation(iph) = inc%Saturation(iph) + relax*incre(iph + NbIncPTC)
!    end do
!
!    ! AccVol
!    do i = 1, NbCompCtilde_ctx(ic)
!        icp = NumCompCtilde_ctx(i, ic)
!        inc%AccVol(icp) = incre(NbIncPTCS + i)
!    end do
!
!    end subroutine IncCV_NewtonIncrement_reservoir

    !> \brief Compute time step (Delta_t) for the next time iteration
    !> DEPRECATED: time step management is to be done in python
    subroutine IncCV_ComputeTimeStep(Delta_t, TimeCurrent)

    double precision, intent(inout) :: Delta_t
    double precision, intent(in) :: TimeCurrent

    double precision :: Delta_tloc, alpha
    integer :: Ierr

    alpha = 1.2d0
    Delta_tloc = Delta_t*alpha

    call MPI_Allreduce(Delta_tloc, Delta_t, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

    Delta_t = min(Delta_t, TimeStepMax)

    end subroutine IncCV_ComputeTimeStep

    !> \brief Save current status if it is necessary to start again current time iteration.
    !!
    !! Copy IncObj to IncObjPreviousTimeStep
    subroutine IncCV_SaveIncPreviousTimeStep() &
        bind(C, name="IncCV_SaveIncPreviousTimeStep")

    call IncCVReservoir_SaveIncPreviousTimeStep
    call IncCVWells_SaveIncPreviousTimeStep

    end subroutine IncCV_SaveIncPreviousTimeStep

    !> \brief Load previous status to start again current time iteration.
    !!
    !! Copy IncObjPreviousTimeStep to IncObj
    subroutine IncCV_LoadIncPreviousTimeStep() &
        bind(C, name="IncCV_LoadIncPreviousTimeStep")

    call IncCVReservoir_LoadIncPreviousTimeStep
    call IncCVWells_LoadIncPreviousTimeStep

    end subroutine IncCV_LoadIncPreviousTimeStep

!    !> \brief Transform IncCV format to a vector,
!    !! necessary for visualization.
!    !!
!    !! Transform Inc into output vector datavisu
!    !! The structure of the vector is                                                 <br>
!    !!   (Pressure, Temperature,                                                      <br>
!    !!        Comp(1), ... , Comp(n),                                                 <br>
!    !!            Saturation(1), ..., Saturation(n))
!    SUBROUTINE IncCV_ToVec_cv(NbIncOwn, Inc, datavisu)
!
!    INTEGER, INTENT(IN) :: NbIncOwn
!    TYPE(Type_IncCV), DIMENSION(:), INTENT(IN) :: Inc
!    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: datavisu
!
!    integer :: i_data, j_data
!    integer :: i, j
!
!    ! Pressure
!    i_data = 1
!    j_data = NbIncOwn
!
!    datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Pression
!
!    ! Temperature
!#ifdef _THERMIQUE_
!    i_data = j_data + 1
!    j_data = j_data + NbIncOwn
!    datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Temperature
!#endif
!
!    ! Comp
!    DO j = 1, NbPhase
!        DO i = 1, NbComp
!
!            IF (MCP(i, j) == 1) THEN
!                i_data = j_data + 1
!                j_data = j_data + NbIncOwn
!                datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Comp(i, j)
!            ENDIF
!        ENDDO
!    ENDDO
!
!    ! Saturation
!    DO i = 1, NbPhase
!        i_data = j_data + 1
!        j_data = j_data + NbIncOwn
!        datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Saturation(i)
!    ENDDO
!    ENDSUBROUTINE IncCV_ToVec_cv
!
!    SUBROUTINE IncCV_ToVec( &
!        datavisucell, &
!        datavisufrac, &
!        datavisunode, &
!        datavisuwellinj, &
!        datavisuwellprod)
!
!    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisucell
!    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisufrac
!    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisunode
!    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellinj
!    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellprod
!
!    INTEGER :: NbCellOwn, NbFracOwn, NbNodeOwn
!
!    NbCellOwn = NbCellOwn_Ncpus(commRank + 1)
!    NbFracOwn = NbFracOwn_Ncpus(commRank + 1)
!    NbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)
!
!    CALL IncCV_ToVec_cv(NbCellOwn, IncCell, datavisucell)
!    CALL IncCV_ToVec_cv(NbFracOwn, IncFrac, datavisufrac)
!    CALL IncCV_ToVec_cv(NbNodeOwn, IncNode, datavisunode)
!
!    datavisuwellinj = 1.d0 ! not implemented
!    datavisuwellprod = 1.d0 ! not implemented
!    ENDSUBROUTINE IncCV_ToVec

    end module IncCV
