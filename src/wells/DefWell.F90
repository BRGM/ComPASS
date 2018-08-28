!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefWell

    use, intrinsic :: iso_c_binding

    use CommonMPI
    use DefModel
    use Physics

implicit none

  type, bind(c) :: WellData_type
     character(c_char) :: &
          IndWell  ! both well types 'p' for pressure mode ; 'f' for flowrate mode
     real(c_double) :: &
          Radius, & ! both well types
          PressionMax, & ! injector only
          PressionMin, & ! producer only
          ImposedFlowrate, & ! both well types (>=0 for producer <0 for injector)
          CompTotal(NbComp), & ! injector only
          Temperature ! injector only
  end type WellData_type

  !> Store data of one Node Well about parent and well index
  type TYPE_DataNodeWell
     integer :: Parent       !< num of parent; -1 if head node
     integer :: PtParent     !< pt of parent; -1 if head node
     double precision :: WID !< Well Index Darcy
     double precision :: WIF !< Well Index Fourier
  end type TYPE_DataNodeWell

  !> CSR type with Pt, Num, and Val (1d DataNodeWell type)
  type TYPE_CSRDataNodeWell
     integer :: Nb
     integer, allocatable, dimension (:) :: Pt
     integer, allocatable, dimension(:) :: Num
     type(TYPE_DataNodeWell), allocatable, dimension (:) :: Val !< Parent and Well Indexes
  end type TYPE_CSRDataNodeWell

  !> to allow = between two DataNodeWell_type
  interface assignment(=)
     module procedure assign_DataNodeWell_equal
  end interface assignment(=)

  !> to allow = between two DataWellInj
  interface assignment(=)
     module procedure assign_DataWell_equal
  end interface assignment(=)

  type(WellData_type), allocatable, target, dimension(:), public :: &
      DataWellInj, &
      DataWellProd

  TYPE(TYPE_CSRDataNodeWell), public :: &
       NodeDatabyWellInj, & !< CSR store data about Parent and Well index of nodes of each injection Well
       NodeDatabyWellProd   !< CSR store data about Parent and Well index of nodes of each production Well

  public :: &
       DefWell_SetDataWellInj,  &
       DefWell_SetDataWellProd, &
       DefWell_Make, &
       DefWell_WellIndex,       &    ! Compute Well index
       DefWell_csrdatawellcopy, &    ! copy CSRDatawell
       DefWell_deallocCSRDataWell, &    ! free CSRdataWell
       get_injectors_data, nb_injectors, &
       get_producers_data, nb_producers

contains

    function get_injectors_data() result(p) &
        bind(C, name="get_injectors_data")

    type(c_ptr) :: p

    if(allocated(DataWellInj)) then
        p = c_loc(DataWellInj(1))
    else
        p = c_null_ptr
    end if

    end function get_injectors_data


    function nb_injectors() result(n) &
        bind(C, name="nb_injectors")

    integer(c_int) :: n

    if(allocated(DataWellInj)) then
        n = size(DataWellInj, 1)
    else
        n = 0
    end if

    end function nb_injectors

    function get_producers_data() result(p) &
        bind(C, name="get_producers_data")

    type(c_ptr) :: p

    if(allocated(DataWellInj)) then
        p = c_loc(DataWellInj(1))
    else
        p = c_null_ptr
    end if

    end function get_producers_data


    function nb_producers() result(n) &
        bind(C, name="nb_producers")

    integer(c_int) :: n

    if(allocated(DataWellInj)) then
        n = size(DataWellInj, 1)
    else
        n = 0
    end if

    end function nb_producers


  subroutine DefWell_print_WellData(datawell)

  type(WellData_type), intent(in) :: datawell

  write(*,*) "%%", "injector data", datawell%Radius, &
	datawell%Temperature, datawell%compTotal(:), &
	datawell%PressionMax, datawell%ImposedFlowrate, &
	datawell%IndWell

  write(*,*) "%%", "producer data", datawell%Radius, &
  datawell%PressionMin, datawell%ImposedFlowrate, &
  datawell%IndWell

  end subroutine DefWell_print_WellData

  ! allocate DataWellInj and set Radius
  subroutine DefWell_SetDataWellInj(NbWell)

    integer, intent(in) :: NbWell
    integer :: k

    allocate(DataWellInj(NbWell))

    do k=1, NbWell

       DataWellInj(k)%Radius = 0.1d0
       DataWellInj(k)%Temperature = 60.d0 + 273.d0

       DataWellInj(k)%CompTotal(:) = 1.d0 ! here NbComp=1
       DataWellInj(k)%PressionMax = 3.d7
       DataWellInj(k)%ImposedFlowrate = - 1.d5/3600.d0

       DataWellInj(k)%IndWell = 'p'
    end do

  end subroutine DefWell_SetDataWellInj


  ! allocate DataWellProd and set Radius
  subroutine DefWell_SetDataWellProd(NbWell)

    integer, intent(in) :: NbWell
    integer :: k

    allocate(DataWellProd(NbWell))

    do k=1, NbWell

       DataWellProd(k)%Radius = 0.1d0
       DataWellProd(k)%PressionMin = 1.d7
       DataWellProd(k)%ImposedFlowrate = 1.d5/3600.d0

       DataWellProd(k)%IndWell = 'p'
    end do

  end subroutine DefWell_SetDataWellProd


  !! ----------------------------------------------------!!

  subroutine DefWell_Make_SetDataWell(NbWellInj, NbWellProd)
    integer, intent(in) :: NbWellInj, NbWellProd

    call DefWell_SetDataWellInj(NbWellInj)   ! allocate DataWellInj and set Radius
    call DefWell_SetDataWellProd(NbWellProd) ! allocate DataWellProd and set Radius

  end subroutine DefWell_Make_SetDataWell

  subroutine DefWell_Make_ComputeWellIndex( &
       NbNode, XNode, CellbyNode, NodebyCell, &
	   FracbyNode, NodebyFace, PermCell, PermFrac)

    integer, intent(in) :: NbNode
    double precision, allocatable, dimension(:,:), intent(in) :: XNode
    type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
    !type(FractureInfoCOC), intent(in) :: FracbyNode
    type(CSR), intent(in) :: FracbyNode

    double precision, allocatable, dimension(:,:,:), intent(in) :: PermCell
    double precision, allocatable, dimension(:), intent(in) :: PermFrac

    integer :: NbWellInj, NbWellProd
    double precision, allocatable, dimension(:) :: WellRadius
    ! FIXME: set consistent values to error codes
    integer :: errcode, Ierr

    if(.NOT.allocated(DataWellProd)) then
        !CHECKME: MPI_Abort is supposed to end all MPI processes
		write(*,*) "ERROR DataWellProd is not allocated."
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if
	NbWellProd = size(DataWellProd)
    if(.NOT.allocated(DataWellInj)) then
        !CHECKME: MPI_Abort is supposed to end all MPI processes
		write(*,*) "ERROR DataWellInj is not allocated."
        call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if
	NbWellInj = size(DataWellInj)

    allocate(WellRadius(max(NbWellInj,NbWellProd)))

    WellRadius(:) = 0
    WellRadius(1:NbWellInj) = DataWellInj(:)%Radius
    call DefWell_WellIndex(NodeDatabyWellInj,NbWellInj,WellRadius, &
         NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
         PermCell, PermFrac)

    WellRadius(:) = 0
    WellRadius(1:NbWellProd) = DataWellProd(:)%Radius
    call DefWell_WellIndex(NodeDatabyWellProd,NbWellProd,WellRadius, &
         NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
         PermCell, PermFrac)

    deallocate(WellRadius)

	end subroutine DefWell_Make_ComputeWellIndex

	subroutine DefWell_Make(NbWellInj, NbWellProd, &
       NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
       PermCell, PermFrac)

    integer, intent(in) :: NbWellInj, NbWellProd, NbNode
    double precision, allocatable, dimension(:,:), intent(in) :: XNode
    type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
    !type(FractureInfoCOC), intent(in) :: FracbyNode
    type(CSR), intent(in) :: FracbyNode

    double precision, allocatable, dimension(:,:,:), intent(in) :: PermCell
    double precision, allocatable, dimension(:), intent(in) :: PermFrac

    double precision, allocatable, dimension(:) :: WellRadius

	call DefWell_Make_SetDataWell(NbWellInj, NbWellProd)

	call DefWell_Make_ComputeWellIndex( &
       NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
       PermCell, PermFrac)

  end subroutine DefWell_Make

  ! Output:
  !  NodeDatabyWell%Val%WI
  ! Use:
  !  NodeDatabyWell%Val%Parent, CellbyNode, FracbyNode, PermCell, PermFrac, Thickness (of frac)
  !  NodebyCell, NodebyFace
  !> \brief Compute the Well Index of every injection and production well using Peaceman formula.
  !!
  !! This model assums that it is derived for a vertical well in a uniform Cartesian grid,
  !! fully penetrating the grid block, with single-phase radial flow and
  !! no interaction with boundaries or other wells.
  !! \WARNING This computation supposes that two wells cannot share a node.
  subroutine DefWell_WellIndex(NodeDatabyWell, NbWell, WellRadius, &
       NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
       PermCell, PermFrac)

    integer, intent(in) :: NbWell, NbNode
    double precision, dimension(:), intent(in) :: WellRadius
    type(TYPE_CSRDataNodeWell), intent(inout) :: NodeDatabyWell

    double precision, allocatable, dimension(:,:), intent(in) :: XNode
    type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
    !type(FractureInfoCOC), intent(in) :: FracbyNode
    type(CSR), intent(in) :: FracbyNode

    double precision, allocatable, dimension(:,:,:), intent(in) :: PermCell
    double precision, allocatable, dimension(:), intent(in) :: PermFrac

    integer :: i, j, k, kp, m, ind, num_node, num_parent, &
         num_cell, num_face, comptCell, comptFrac
    logical :: cell_edge
    double precision :: meanDist, meanPerm, meanThickness, dr0, de, wi, length
    double precision, dimension(3) :: xn1, xn2, xk
    double precision, allocatable, dimension(:) :: WI_global

    double precision, parameter :: Pi = 3.14159265359d0

    allocate(WI_global(NbNode))
    WI_global(:) = 0.d0

    ! computation of Well Index Darcy

    !! MATRIX
    do i=1,NbWell
       ! loop over the EDGES: (node,Parent(node)) for node=1,head_node-1
       do j=NodeDatabyWell%Pt(i)+1,NodeDatabyWell%Pt(i+1)-1
          num_node = NodeDatabyWell%Num(j)
          num_parent = NodeDatabyWell%Val(j)%Parent

          xn1(1:3) = XNode(1:3,num_node)
          xn2(1:3) = XNode(1:3,num_parent)

          comptCell = 0
          meanDist = 0.d0
          meanPerm = 0.d0
          do k=CellbyNode%Pt(num_node)+1,CellbyNode%Pt(num_node+1)  ! cells containing node
             cell_edge = .false.
             num_cell = CellbyNode%Num(k)
             do kp=CellbyNode%Pt(num_parent)+1,CellbyNode%Pt(num_parent+1)  ! cells containing parent
                if(num_cell == CellbyNode%Num(kp)) then   ! cells containing Node and Parent, then containing the edge
                   cell_edge = .true.
                   exit
                endif
             enddo
             if(cell_edge) then
                comptCell = comptCell + 1

                ! find coordinates of center of cell
                xk(:) = 0.d0
                do m = NodebyCell%Pt(num_cell)+1, NodebyCell%Pt(num_cell+1)
                   xk(:) = xk(:) + XNode(:, NodebyCell%Num(m))
                enddo
                xk(:) = xk(:)/dble(NodebyCell%Pt(num_cell+1) - NodebyCell%Pt(num_cell))
                call DistNodetoLine(xk, xn1, xn2, length) ! dist of cell center to the edge (num,parent)
                meanDist = meanDist + length

                meanPerm = meanPerm + PermCell(1,1,num_cell)  ! this formula is true if perm iso !
             endif
          enddo ! loop over cells


          meanPerm = meanPerm / dble(comptCell)
          meanDist = meanDist / dble(comptCell)

          de = dsqrt(dot_product( xn1-xn2 , xn1-xn2 ))   ! length of edge
          dr0 = 0.14036d0*dsqrt(2.d0)*meanDist   ! Peaceman radius
          wi = de * Pi * meanPerm / log(dr0/WellRadius(i))

          ! contribution of the edge to each nodes: node and parent
          WI_global(num_node) = WI_global(num_node) + wi
          WI_global(num_parent) = WI_global(num_parent) + wi

       enddo ! loop over edges
    enddo ! loop over wells

    !! FRACTURES
    do i=1, NbWell

       ! loop over the nodes
       do j=NodeDatabyWell%Pt(i)+1,NodeDatabyWell%Pt(i+1)

          num_node = NodeDatabyWell%Num(j)
          xn1(1:3) = XNode(1:3,num_node)

          meanDist = 0.d0
          meanPerm = 0.d0
          meanThickness = 0.d0
          comptFrac = 0

          ! loop over frac of node
          do k=FracbyNode%Pt(num_node)+1, FracbyNode%Pt(num_node+1)

             comptFrac = comptFrac + 1
             num_face = FracbyNode%Num(k)

             ! find coordinates of center of frac
             xk(:) = 0.d0
             do m = NodebyFace%Pt(num_face)+1, NodebyFace%Pt(num_face+1)
                xk(:) = xk(:) + XNode(:, NodebyFace%Num(m))
             enddo
             xk(:) = xk(:)/dble(NodebyFace%Pt(num_face+1) - NodebyFace%Pt(num_face))
             length = dsqrt( dot_product(xk-xn1, xk-xn1) ) ! dist of cell frac to the node
             meanDist = meanDist + length

             meanPerm = meanPerm + PermFrac(num_face)  ! this formula is true if perm iso !
             meanThickness = meanThickness + Thickness
          enddo
          if(comptFrac>0)then
             ! there is at least one frac in this node, compute wi
             meanPerm = meanPerm / dble(comptFrac)
             meanDist = meanDist / dble(comptFrac)
             meanThickness = meanThickness / dble(comptFrac)

             dr0 = 0.14036d0*dsqrt(2.d0)*meanDist   ! Peaceman radius
             wi = meanThickness * 2.d0 * Pi * meanPerm / log(dr0/WellRadius(i))

             ! contribution of the fracs to the node
             WI_global(num_node) = WI_global(num_node) + wi
          endif
       enddo
    enddo

    ! fill %WIF and %WID
    NodeDatabyWell%Val(:)%WIF = 0.d0  ! not used, negligeable compared to WID

    ! loop over all nodes of all wells
    do j=1,NodeDatabyWell%Pt(NbWell+1)
       num_node = NodeDatabyWell%Num(j)
       NodeDatabyWell%Val(j)%WID = WI_global(num_node)
    enddo

    deallocate(WI_global)

  end subroutine DefWell_WellIndex


  !> \brief Compute the distance from A to the Line (BC)
  !! if B and C are too close ( BC<E-10) compute distance AB.
  subroutine DistNodetoLine(NodeA, NodeB, NodeC, length)

    double precision, dimension(:), intent(in) :: NodeA, NodeB, NodeC
    double precision, intent(out) :: length

    double precision, dimension(3) :: cross_product, vectBA, vectBC

    vectBA = NodeA - NodeB
    vectBC = NodeC - NodeB

    length = dot_product(vectBC, vectBC)
    if(length<1.0d-5) then
       write(*,*)"pb with DistNodetoLine, line(BC) is not a line because B=C"
       length = dsqrt( dot_product(vectBA, vectBA) )
       return
    endif

    ! cross product(BA, BC)
    cross_product(1) = vectBA(2)*vectBC(3)-vectBA(3)*vectBC(2)
    cross_product(2) = vectBA(3)*vectBC(1)-vectBA(1)*vectBC(3)
    cross_product(3) = vectBA(1)*vectBC(2)-vectBA(2)*vectBC(1)

    ! dist(A,BC) = norm( cross product(BA, BC) ) / norm(BC)
    length = dsqrt( dot_product(cross_product,cross_product) / length )

  end subroutine DistNodetoLine


  !> \brief Copy CSRDataWell1 to CSRDataWell2
  subroutine DefWell_csrdatawellcopy(CSR1, CSR2)

    type(TYPE_CSRDataNodeWell), intent(in)  :: CSR1
    type(TYPE_CSRDataNodeWell), intent(out) :: CSR2

    integer :: Nnz, i

    CSR2%Nb = CSR1%Nb

    allocate(CSR2%Pt(CSR1%Nb+1))
    CSR2%Pt(:) = CSR1%Pt(:)

    Nnz = CSR1%Pt(CSR1%Nb+1)

    allocate( CSR2%Num(Nnz))
    do i=1, Nnz
       CSR2%Num(i) = CSR1%Num(i)
    end do

    allocate(CSR2%Val(Nnz))
    do i=1, Nnz
       CSR2%Val(i) = CSR1%Val(i)
    end do

  end subroutine DefWell_csrdatawellcopy

  !> \brief Deallocate CSRDataWell (\%Pt, \%Num and \%Val)
  subroutine DefWell_deallocCSRDataWell(CSR1)

    type(TYPE_CSRDataNodeWell), intent(inout) :: CSR1

    deallocate(CSR1%Pt)

    if(allocated(CSR1%Num)) then
       deallocate(CSR1%Num)
    end if

    if (allocated(CSR1%Val)) then
       deallocate(CSR1%Val)
    end if

  end subroutine DefWell_deallocCSRDataWell


  !> \brief Define operator = between two DataNodeWell_type:  x2 = x1
  subroutine assign_DataNodeWell_equal(x2, x1)

    type(TYPE_DataNodeWell), intent(in) :: x1
    type(TYPE_DataNodeWell), intent(out) :: x2

    x2%Parent = x1%Parent
    x2%PtParent = x1%PtParent
    x2%WID = x1%WID
    x2%WIF = x1%WIF

  end subroutine assign_DataNodeWell_equal

  !> \brief Define operator = between two DataWellInj:  x2 = x1
  subroutine assign_DataWell_equal(x2, x1)

    TYPE(WellData_type), intent(in) :: x1
    TYPE(WellData_type), intent(out) :: x2

    x2%IndWell = x1%IndWell
    x2%Radius = x1%Radius
    x2%PressionMax = x1%PressionMax
    x2%PressionMin = x1%PressionMin
    x2%ImposedFlowrate = x1%ImposedFlowrate
    x2%Temperature = x1%Temperature
    x2%CompTotal = x1%CompTotal
    
  end subroutine assign_DataWell_equal

end module DefWell