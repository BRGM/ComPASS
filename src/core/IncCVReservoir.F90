    !
    ! This file is part of ComPASS.
    !
    ! ComPASS is free software: you can redistribute it and/or modify it under both the terms
    ! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
    ! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
    !

    module IncCVReservoir

    use MeshSchema
    use DefModel
    use Thermodynamics

    use NumbyContext
    use CommonMPI
    use Physics
    use SchemeParameters

    use iso_c_binding

    implicit none

    !> Unknown for Degree Of Freedom (including thermal). DOF can be Cell, Fracture Face or Node.
    TYPE TYPE_IncCVReservoir

        integer(c_int) :: ic !< context (???)

        real(c_double) :: & ! values of Inc
        Pression, & !< Pressure of the element
        Temperature, & !< Temperature of the element
        Comp(NbComp, NbPhase), & !< Molar composition of the element
        Saturation(NbPhase), & !< Saturation of the element
        AccVol(NbCompThermique) !< ??? of the element

    end TYPE TYPE_IncCVReservoir

    !> to allow = between two TYPE_IncCVReservoir
    interface assignment(=)
    module procedure assign_type_inccv
    end interface assignment(=)

    ! Inc for current time step: cell, fracture faces and nodes
    TYPE(TYPE_IncCVReservoir), allocatable, dimension(:), target, public :: &
        IncCell, & !< Cell unknowns for current time step
    IncFrac, & !< Fracture Face unknowns for current time step
    IncNode !< Node unknowns for current time step

    ! Inc for previous time step: current time step - 1
    TYPE(TYPE_IncCVReservoir), allocatable, dimension(:), public :: &
        IncCellPreviousTimeStep, & !< Cell unknowns for previous time step
    IncFracPreviousTimeStep, & !< Fracture Face unknowns for previous time step
    IncNodePreviousTimeStep !< Node unknowns for previous time step

    public :: &
        IncCVReservoir_allocate, &
        IncCVReservoir_NewtonRelax, &
        IncCVReservoir_NewtonIncrement, &
        IncCVReservoir_LoadIncPreviousTimeStep, &
        IncCVReservoir_SaveIncPreviousTimeStep, &
        IncCVReservoir_free

private :: &
        IncCVReservoir_NewtonIncrement_reservoir

    contains

    !> \brief Define operator = between two TYPE_IncCV:  inc2=inc1
    subroutine assign_type_inccv(inc2, inc1)

    type(TYPE_IncCVReservoir), intent(in) :: inc1
    type(TYPE_IncCVReservoir), intent(out) :: inc2

    inc2%ic = inc1%ic

    inc2%Pression = inc1%Pression
#ifdef _THERMIQUE_
    inc2%Temperature = inc1%Temperature
#endif
    inc2%Comp(:, :) = inc1%Comp(:, :)
    inc2%Saturation(:) = inc1%Saturation(:)
    inc2%AccVol(:) = inc1%AccVol(:)

    end subroutine assign_type_inccv

    subroutine IncCVReservoir_allocate

    integer :: Nb, Nnz

    allocate (IncCell(NbCellLocal_Ncpus(commRank + 1)))
    allocate (IncFrac(NbFracLocal_Ncpus(commRank + 1)))
    allocate (IncNode(NbNodeLocal_Ncpus(commRank + 1)))

    allocate (IncCellPreviousTimeStep(NbCellLocal_Ncpus(commRank + 1)))
    allocate (IncFracPreviousTimeStep(NbFracLocal_Ncpus(commRank + 1)))
    allocate (IncNodePreviousTimeStep(NbNodeLocal_Ncpus(commRank + 1)))

    end subroutine IncCVReservoir_allocate

    !> \brief Deallocate unknowns vectors
    subroutine IncCVReservoir_free

    deallocate (IncCell)
    deallocate (IncFrac)
    deallocate (IncNode)

    deallocate (IncCellPreviousTimeStep)
    deallocate (IncFracPreviousTimeStep)
    deallocate (IncNodePreviousTimeStep)

    end subroutine IncCVReservoir_free

    subroutine IncCVReservoir_NewtonIncrement( &
        NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
        relax)

    double precision, dimension(:, :), intent(in) :: &
        NewtonIncreNode, &
        NewtonIncreFrac, &
        NewtonIncreCell

    double precision, intent(in) :: relax

    integer :: k

    ! nodes
    do k = 1, NbNodeLocal_Ncpus(commRank + 1)
        call IncCVReservoir_NewtonIncrement_reservoir(IncNode(k), NewtonIncreNode(:, k), relax)
    end do

    ! fracture faces
    do k = 1, NbFracLocal_Ncpus(commRank + 1)
        call IncCVReservoir_NewtonIncrement_reservoir(IncFrac(k), NewtonIncreFrac(:, k), relax)
    end do

    ! cells
    do k = 1, NbCellLocal_Ncpus(commRank + 1)
        call IncCVReservoir_NewtonIncrement_reservoir(IncCell(k), NewtonIncreCell(:, k), relax)
        !     write(*,*) ' increment cell ',k,NewtonIncreCell(:,k)
    end do

    end subroutine IncCVReservoir_NewtonIncrement

    !> \brief Compute relaxation in Newton.
    !!
    !! relax = min(1, IncreObj/NewtonIncreObjMax)                   <br>
    !! where IncreObj is set by the user in DefModel.F90            <br>
    !! and NewtonIncreObjMax is the maximum of the Nemton increment
    !! in current iteration
    subroutine IncCVReservoir_NewtonRelax( &
        NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, relax)

    double precision, dimension(:, :), intent(in) :: &
        NewtonIncreNode, &
        NewtonIncreFrac, &
        NewtonIncreCell

    double precision, intent(out) :: relax

    double precision :: &
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

    end subroutine IncCVReservoir_NewtonRelax

    !> \brief Realize Newton increment of each control volume
    subroutine IncCVReservoir_NewtonIncrement_reservoir(inc, incre, relax)

    type(TYPE_IncCVReservoir), intent(inout) :: inc
    double precision, intent(in) :: incre(NbIncPTCSMax), relax

    integer :: i, icp, iph
    integer :: NbIncPTC, NbIncPTCS
    integer :: ic

    ic = inc%ic
    NbIncPTC = NbIncPTC_ctx(ic)
    NbIncPTCS = NbIncPTC_ctx(ic) + NbPhasePresente_ctx(ic)

    ! increment Pressure
    inc%Pression = inc%Pression + relax*incre(1)

    !    write(*,*)' increment P ',relax,incre(1)

#ifdef _THERMIQUE_

    ! increment Temperature
    inc%Temperature = inc%Temperature + relax*incre(2)
#endif

    ! increment comp
    do i = 2 + IndThermique, NbIncPTC

        icp = NumIncPTC2NumIncComp_comp_ctx(i, ic)
        iph = NumIncPTC2NumIncComp_phase_ctx(i, ic)
        inc%Comp(icp, iph) = inc%Comp(icp, iph) + relax*incre(i)
    enddo

    ! increment saturation
    do i = 1, NbPhasePresente_ctx(ic)
        iph = NumPhasePresente_ctx(i, ic)

        inc%Saturation(iph) = inc%Saturation(iph) + relax*incre(iph + NbIncPTC)
    end do

    ! AccVol
    do i = 1, NbCompCtilde_ctx(ic)
        icp = NumCompCtilde_ctx(i, ic)
        inc%AccVol(icp) = incre(NbIncPTCS + i)
    end do

    end subroutine IncCVReservoir_NewtonIncrement_reservoir

    !> \brief Save current status if it is necessary to start again current time iteration.
    !!
    !! Copy IncObj to IncObjPreviousTimeStep
    subroutine IncCVReservoir_SaveIncPreviousTimeStep

    integer :: k

    ! save current status
    do k = 1, NbNodeLocal_Ncpus(commRank + 1)
        IncNodePreviousTimeStep(k) = IncNode(k)
    end do
    do k = 1, NbFracLocal_Ncpus(commRank + 1)
        IncFracPreviousTimeStep(k) = IncFrac(k)
    end do
    do k = 1, NbCellLocal_Ncpus(commRank + 1)
        IncCellPreviousTimeStep(k) = IncCell(k)
    end do

    end subroutine IncCVReservoir_SaveIncPreviousTimeStep

    !> \brief Load previous status to start again current time iteration.
    !!
    !! Copy IncObjPreviousTimeStep to IncObj
    subroutine IncCVReservoir_LoadIncPreviousTimeStep

    integer :: k

    do k = 1, NbNodeLocal_Ncpus(commRank + 1)
        IncNode(k) = IncNodePreviousTimeStep(k)
    end do
    do k = 1, NbFracLocal_Ncpus(commRank + 1)
        IncFrac(k) = IncFracPreviousTimeStep(k)
    end do
    do k = 1, NbCellLocal_Ncpus(commRank + 1)
        IncCell(k) = IncCellPreviousTimeStep(k)
    end do

    end subroutine IncCVReservoir_LoadIncPreviousTimeStep

    !> \brief Transform IncCVReservoir format to a vector,
    !! necessary for visualization.
    !!
    !! Transform Inc into output vector datavisu
    !! The structure of the vector is                                                 <br>
    !!   (Pressure, Temperature,                                                      <br>
    !!        Comp(1), ... , Comp(n),                                                 <br>
    !!            Saturation(1), ..., Saturation(n))
    SUBROUTINE IncCVReservoir_ToVec_cv(NbIncOwn, Inc, datavisu)

    INTEGER, INTENT(IN) :: NbIncOwn
    TYPE(Type_IncCVReservoir), DIMENSION(:), INTENT(IN) :: Inc
    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: datavisu

    integer :: i_data, j_data
    integer :: i, j

    ! Pressure
    i_data = 1
    j_data = NbIncOwn

    datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Pression

    ! Temperature
#ifdef _THERMIQUE_
    i_data = j_data + 1
    j_data = j_data + NbIncOwn
    datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Temperature
#endif

    ! Comp
    DO j = 1, NbPhase
        DO i = 1, NbComp

            IF (MCP(i, j) == 1) THEN
                i_data = j_data + 1
                j_data = j_data + NbIncOwn
                datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Comp(i, j)
            ENDIF
        ENDDO
    ENDDO

    ! Saturation
    DO i = 1, NbPhase
        i_data = j_data + 1
        j_data = j_data + NbIncOwn
        datavisu(i_data:j_data) = Inc(1:NbIncOwn)%Saturation(i)
    ENDDO
    ENDSUBROUTINE IncCVReservoir_ToVec_cv

    SUBROUTINE IncCVReservoir_ToVec( &
        datavisucell, &
        datavisufrac, &
        datavisunode, &
        datavisuwellinj, &
        datavisuwellprod)

    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisucell
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisufrac
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisunode
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellinj
    DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: datavisuwellprod

    INTEGER :: NbCellOwn, NbFracOwn, NbNodeOwn

    NbCellOwn = NbCellOwn_Ncpus(commRank + 1)
    NbFracOwn = NbFracOwn_Ncpus(commRank + 1)
    NbNodeOwn = NbNodeOwn_Ncpus(commRank + 1)

    CALL IncCVReservoir_ToVec_cv(NbCellOwn, IncCell, datavisucell)
    CALL IncCVReservoir_ToVec_cv(NbFracOwn, IncFrac, datavisufrac)
    CALL IncCVReservoir_ToVec_cv(NbNodeOwn, IncNode, datavisunode)

    datavisuwellinj = 1.d0 ! not implemented
    datavisuwellprod = 1.d0 ! not implemented
    ENDSUBROUTINE IncCVReservoir_ToVec

    end module IncCVReservoir