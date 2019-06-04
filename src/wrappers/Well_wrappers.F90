!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module WellWrapper

    use, intrinsic :: iso_c_binding
    use mpi, only: MPI_Abort
    use CommonMPI, only: ComPASS_COMM_WORLD, Ncpus, CommonMPI_abort
    use CommonTypesWrapper, only: cpp_COC
    use InteroperabilityStructures, only: cpp_array_wrapper
    use DefModel, only: NbComp
    use DefWell, only: DataWellProd, DataWellInj
    use GlobalMesh, only: &
      NbEdgebyWellInj, NbEdgebyWellProd, &
      NbWellInj, NbWellProd, &
      NumNodebyEdgebyWellInj, NumNodebyEdgebyWellProd

    implicit none

    type, bind(C) :: Producer_data
        real(c_double) :: radius
        real(c_double) :: minimum_pressure
        real(c_double) :: imposed_flowrate
        character(c_char) :: operating_code  ! 'p' for pressure mode ; 'f' for flowrate mode
    end type Producer_data

    type, bind(C) :: Injector_data
        real(c_double) :: radius
        real(c_double) :: temperature
        ! FIXME: introduce composition (CompTotal(NbComp))
        real(c_double) :: maximum_pressure
        real(c_double) :: imposed_flowrate
        character(c_char) :: operating_code  ! 'p' for pressure mode ; 'f' for flowrate mode
    end type Injector_data

    public :: &
        Well_allocate_specific_well_geometries, &
        Well_set_producers_data, &
        Well_set_injectors_data, &
        Well_allocate_well_geometries_from_C, &
        Well_set_wells_data_from_C

    contains

    subroutine Well_allocate_specific_well_geometries(geometries, NbEdgebyWell, NumNodebyEdgebyWell)

    type(cpp_COC), intent(in) :: geometries
    integer, allocatable, dimension(:), intent(inout) :: NbEdgebyWell
    integer, allocatable, dimension(:,:,:), intent(inout) :: NumNodebyEdgebyWell
    integer(c_int), pointer :: offsets(:), edges(:)
    integer i, j, k, nb_wells, nb_nodes, nb_edges

    ! FIXME: Memory management should not be here
    if(allocated(NbEdgebyWell)) then
        deallocate(NbEdgebyWell)
    end if
    ! FIXME: Memory management should not be here
    if(allocated(NumNodebyEdgebyWell)) then
        deallocate(NumNodebyEdgebyWell)
    end if

    nb_wells = geometries%nb_containers
    allocate(NbEdgebyWell(nb_wells))
    call c_f_pointer(geometries%container_offset, offsets, shape=[nb_wells + 1])
    do i=1, nb_wells
        nb_nodes = offsets(i+1) - offsets(i)
        if(mod(nb_nodes, 2)/=0) &
            call CommonMPI_abort('inconsistent well geometry')
        nb_edges = nb_nodes / 2
        NbEdgebyWell(i) = nb_edges
    end do
    allocate(NumNodebyEdgebyWell(2, maxval(NbEdgebyWell), nb_wells))
    call c_f_pointer(geometries%container_content, edges, shape=[offsets(nb_wells + 1)])
    k = 1
    do i=1, nb_wells
        do j=1,  NbEdgebyWell(i)
            NumNodebyEdgebyWell(1:2, j, i) = edges(k:k+1)
            k = k + 2
        end do
    end do
    if(k/=offsets(nb_wells + 1)+1) &
        call CommonMPI_abort('inconsistent well geometry')

    end subroutine Well_allocate_specific_well_geometries

    subroutine Well_allocate_well_geometries_from_C(producers_geometries, injectors_geometries) &
        bind(C, name="Well_allocate_well_geometries")

    type(cpp_COC), intent(in) :: producers_geometries
    type(cpp_COC), intent(in) :: injectors_geometries

    !write(*,*) "FORTRAN - Producers"
    call Well_allocate_specific_well_geometries(producers_geometries, NbEdgebyWellProd, NumNodebyEdgebyWellProd)
    ! FIXME: Setting here a global variable
    NbWellProd = size(NbEdgebyWellProd)
    !write(*,*) "NbEdgebyWellProd:", NbEdgebyWellProd
    !write(*,*) "NumNodebyEdgebyWellProd:", NumNodebyEdgebyWellProd

    !write(*,*) "FORTRAN - Injectors"
    call Well_allocate_specific_well_geometries(injectors_geometries, NbEdgebyWellInj, NumNodebyEdgebyWellInj)
    ! FIXME: Setting here a global variable
    NbWellInj = size(NbEdgebyWellInj)
    !write(*,*) "NbEdgebyWellInj:", NbEdgebyWellInj
    !write(*,*) "NumNodebyEdgebyWellInj:", NumNodebyEdgebyWellInj

    end subroutine Well_allocate_well_geometries_from_C

    subroutine Well_set_producers_data(c_producers_data)

    type(cpp_array_wrapper), intent(in) :: c_producers_data
    integer(c_size_t) :: k, nb_producers
    type(Producer_data), dimension(:), pointer :: producers_data

    nb_producers = c_producers_data%n
    allocate(DataWellProd(nb_producers))
    call c_f_pointer(c_producers_data%p, producers_data, shape=[nb_producers])
    !write(*,*) "decode", nb_producers, "producers"
    do k=1, nb_producers
        DataWellProd(k)%Radius = producers_data(k)%radius
        DataWellProd(k)%PressionMin = producers_data(k)%minimum_pressure
        DataWellProd(k)%ImposedFlowrate = producers_data(k)%imposed_flowrate
        DataWellProd(k)%IndWell = producers_data(k)%operating_code
        !write(*,*) "Setting producer:", &
        !    DataWellProd(k)%Radius, &
        !    DataWellProd(k)%PressionMin, DataWellProd(k)%ImposedFlowrate, &
        !    DataWellProd(k)%IndWell
    end do

    end subroutine Well_set_producers_data

    subroutine Well_set_injectors_data(c_injectors_data)

    type(cpp_array_wrapper), intent(in) :: c_injectors_data
    integer(c_size_t) :: k, nb_injectors
    type(Injector_data), dimension(:), pointer :: injectors_data

    nb_injectors = c_injectors_data%n
    allocate(DataWellInj(nb_injectors))
    call c_f_pointer(c_injectors_data%p, injectors_data, shape=[nb_injectors])
    !write(*,*) "decode", nb_injectors, "injectors"
    do k=1, nb_injectors
        DataWellInj(k)%Radius = injectors_data(k)%radius
        DataWellInj(k)%CompTotal(:) = 1.d0 ! FIXME... for multi component injection
        DataWellInj(k)%Temperature = injectors_data(k)%temperature
        DataWellInj(k)%PressionMax = injectors_data(k)%maximum_pressure
        DataWellInj(k)%ImposedFlowrate = injectors_data(k)%imposed_flowrate
        DataWellInj(k)%IndWell = injectors_data(k)%operating_code
        !write(*,*) "Setting injector:", &
        !    DataWellInj(k)%Radius, DataWellInj(k)%Temperature, &
        !    DataWellInj(k)%PressionMax, DataWellInj(k)%ImposedFlowrate, &
        !    DataWellInj(k)%IndWell
    end do

    end subroutine Well_set_injectors_data

    subroutine Well_set_wells_data_from_C(c_producers_data, c_injectors_data) &
        bind(C, name="Well_set_wells_data")

    type(cpp_array_wrapper), intent(in) :: c_producers_data, c_injectors_data

    call Well_set_producers_data(c_producers_data)
    call Well_set_injectors_data(c_injectors_data)

    end subroutine Well_set_wells_data_from_C

    end module WellWrapper
