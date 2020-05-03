!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module DefWell

   use, intrinsic :: iso_c_binding, only: &
     c_ptr, c_size_t, c_null_ptr, c_loc, c_int, c_double, c_char, c_bool

   use mpi, only: &
     MPI_CHARACTER, &
     MPI_INT, &
     MPI_DOUBLE, &
     MPI_Type_Create_Struct, &
     MPI_ADDRESS_KIND, &
     MPI_Abort

   use mpi_f08, only: &
     MPI_Get_address

   use CommonType, only: CSR
   use CommonMPI, only: ComPASS_COMM_WORLD, CommonMPI_abort
   use DefModel, only: NbComp
   use Physics, only: thickness

   implicit none

   type, bind(c) :: WellData_type
      integer(c_int) :: Id
      real(c_double) :: &
         Radius, & ! both well types
         PressionMax, & ! injector only
         PressionMin, & ! producer only
         ImposedFlowrate, & ! both well types (>=0 for producer <0 for injector)
         CompTotal(NbComp), & ! injector only
         InjectionTemperature, & ! injector only
         actual_mass_flowrate, &
         actual_energy_flowrate, &
         actual_pressure, &
         actual_temperature
      ! WARNING: we put character at the end of the structure
      ! because of "memory padding" when creating mpi well data structure
      ! cf. DefWell_mpi_register_well_data_description
      character(c_char) :: &
         IndWell ! both well types 'p' for pressure mode ; 'f' for flowrate mode; 'c' for closed
   end type WellData_type

   !> Store data of one Node Well about parent and well index
   type, bind(C) :: TYPE_DataNodeWell
      integer(c_int) :: Parent !< num of parent; -1 if head node
      integer(c_int) :: PtParent !< pt of parent; -1 if head node ! FIXME: improve doc !!!
      real(c_double) :: WID !< Well Index Darcy
      real(c_double) :: WIF !< Well Index Fourier
   end type TYPE_DataNodeWell

   !> CSR type with Pt, Num, and Val (1d DataNodeWell type)
   type TYPE_CSRDataNodeWell
      integer(c_int) :: Nb
      integer(c_int), pointer, dimension(:) :: Pt
      integer(c_int), pointer, dimension(:) :: Num
      type(TYPE_DataNodeWell), pointer, dimension(:) :: Val !< Parent and Well Indexes
   end type TYPE_CSRDataNodeWell

   type, bind(C) :: PerforationDataCSR_wrapper
      integer(c_int) :: nb_wells
      type(c_ptr) :: well_offset
      type(c_ptr) :: node_vertex
      type(c_ptr) :: data
   end type PerforationDataCSR_wrapper   

   !> to allow = between two DataNodeWell_type
   interface assignment(=)
      module procedure assign_DataNodeWell_equal
   end interface assignment(=)
      
   !> to allow = between two DataWellInj
   interface assignment(=)
      module procedure assign_DataWell_equal
   end interface assignment(=)

   logical(c_bool) :: DefWell_has_WI_threshold = .false.
   real(c_double) :: DefWell_WI_threshold = 0.d0

   ! FIXME: switch to pointer
   type(WellData_type), allocatable, public, target, dimension(:) :: &
      DataWellInj, &
      DataWellProd

   TYPE(TYPE_CSRDataNodeWell), public :: &
      NodeDatabyWellInj, & !< CSR store data about Parent and Well index of nodes of each injection Well
      NodeDatabyWellProd !< CSR store data about Parent and Well index of nodes of each production Well

   public :: &
      DefWell_WellIndex, & ! Compute Well index
      DefWell_csrdatawellcopy, & ! copy CSRDatawell
      DefWell_deallocCSRDataWell, & ! free CSRdataWell
      get_global_injectors_data, nb_global_injectors, &
      get_global_producers_data, nb_global_producers, &
      DefWell_mpi_register_well_data_description, &
      DefWell_set_WI_threshold, &
      DefWell_unset_WI_threshold

contains

   subroutine DefWell_set_WI_threshold(threshold) &
        bind(C, name="set_Peaceman_WI_threshold")
        
        real(c_double), intent(in), value :: threshold

        if(threshold<=0.d0) then
                call CommonMPI_abort("peaceaman index cannot be negative") 
        endif 
        DefWell_has_WI_threshold = .true.
        DefWell_WI_threshold = threshold
   
   end subroutine DefWell_set_WI_threshold

   subroutine DefWell_unset_WI_threshold() &
        bind(C, name="unset_Peaceman_WI_threshold")
   
        DefWell_has_WI_threshold = .false.

   end subroutine DefWell_unset_WI_threshold


   function get_global_injectors_data() result(p) &
      bind(C, name="get_global_injectors_data")

      type(c_ptr) :: p

      if (allocated(DataWellInj)) then
         p = c_loc(DataWellInj(1))
      else
         p = c_null_ptr
      end if

   end function get_global_injectors_data

   function nb_global_injectors() result(n) &
      bind(C, name="nb_global_injectors")

      integer(c_size_t) :: n

      if (allocated(DataWellInj)) then
         n = size(DataWellInj, 1)
      else
         n = 0
      end if

   end function nb_global_injectors

   function get_global_producers_data() result(p) &
      bind(C, name="get_global_producers_data")

      type(c_ptr) :: p

      if (allocated(DataWellProd)) then
         p = c_loc(DataWellProd(1))
      else
         p = c_null_ptr
      end if

   end function get_global_producers_data

   function nb_global_producers() result(n) &
      bind(C, name="nb_global_producers")

      integer(c_size_t) :: n

      if (allocated(DataWellProd)) then
         n = size(DataWellProd, 1)
      else
         n = 0
      end if

   end function nb_global_producers

   subroutine DefWell_Make_ComputeWellIndex( &
      NbNode, XNode, CellbyNode, NodebyCell, &
      FracbyNode, NodebyFace, PermCell, PermFrac)

      integer, intent(in) :: NbNode
      double precision, allocatable, dimension(:, :), intent(in) :: XNode
      type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
      !type(FractureInfoCOC), intent(in) :: FracbyNode
      type(CSR), intent(in) :: FracbyNode

      double precision, allocatable, dimension(:, :, :), intent(in) :: PermCell
      double precision, allocatable, dimension(:), intent(in) :: PermFrac

      integer :: NbWellInj, NbWellProd
      double precision, allocatable, dimension(:) :: WellRadius
      ! FIXME: set consistent values to error codes
      integer :: errcode, Ierr

      if (.NOT. allocated(DataWellProd)) then
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         write (*, *) "ERROR DataWellProd is not allocated."
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
      NbWellProd = size(DataWellProd)
      if (.NOT. allocated(DataWellInj)) then
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         write (*, *) "ERROR DataWellInj is not allocated."
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
      NbWellInj = size(DataWellInj)

      allocate (WellRadius(max(NbWellInj, NbWellProd)))

      WellRadius(:) = 0
      WellRadius(1:NbWellInj) = DataWellInj(:)%Radius
      call DefWell_WellIndex(NodeDatabyWellInj, NbWellInj, WellRadius, &
                             NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
                             PermCell, PermFrac)

      WellRadius(:) = 0
      WellRadius(1:NbWellProd) = DataWellProd(:)%Radius
      call DefWell_WellIndex(NodeDatabyWellProd, NbWellProd, WellRadius, &
                             NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
                             PermCell, PermFrac)

      deallocate (WellRadius)

   end subroutine DefWell_Make_ComputeWellIndex

   ! FIXME: this is a convenience function that should be elsewhere
   subroutine element_center(vertices, connectivity, element, xc)
      double precision, intent(in) :: vertices(:, :)
      type(CSR), intent(in) :: connectivity
      integer, intent(in) :: element
      double precision, intent(inout) :: xc(3)
      
      integer :: m

      xc = 0.d0
      do m = connectivity%Pt(element) + 1, connectivity%Pt(element + 1)
         xc = xc + vertices(:, connectivity%Num(m))
      enddo
      xc = xc / dble( connectivity%Pt(element + 1) - (connectivity%Pt(element)))
   
   end subroutine element_center

   subroutine compute_peaceman_indices(                   &
      sum_distances, sum_permeabilities, nb_contributors, &
      thickness, well_radius, well_index                  &
   )
      double precision, intent(in) :: sum_distances, sum_permeabilities
      integer, intent(in) :: nb_contributors
      double precision, intent(in) :: well_radius, thickness
      double precision, intent(out) :: well_index
      
      double precision, parameter :: Pi = 3.14159265359d0
      double precision :: d, k, peaceman_radius
      double precision :: wi_max

      d = sum_distances / dble( nb_contributors )
      k = sum_permeabilities / dble( nb_contributors )
#ifndef NDEBUG
      if(d<=0) call CommonMPI_abort('negative distance')
      if(k<=0) call CommonMPI_abort('negative permeability')
#endif
      peaceman_radius = 0.14036d0 * dsqrt(2.d0) * d

      well_index = thickness * Pi * k / log(peaceman_radius / well_radius)
            
      if(peaceman_radius < well_radius) &
         call CommonMPI_abort('Peaceman Well radius is smaller than effective well radius (negative well index).')
      
      if(DefWell_has_WI_threshold) then
         if(peaceman_radius < 2 * well_radius) then
             wi_max = thickness * Pi * k / log(2d0)
             write(*,*) 'WARNING'
             write(*,*) 'WARNING'
             write(*,*) ''
             write(*,*) 'Applying threshold on Peaceman Well Index'
             write(*,*) 'well radius', well_radius, 'vs. Peaceman radius', peaceman_radius
             write(*,*) 'well index', well_index, 'well index limit', wi_max
             write(*,*) ''
             write(*,*) 'WARNING'
             write(*,*) 'WARNING'
             well_index = wi_max
         end if
       end if

   end subroutine compute_peaceman_indices

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

      double precision, allocatable, dimension(:, :), intent(in) :: XNode
      type(CSR), intent(in) :: CellbyNode, NodebyCell, NodebyFace
      !type(FractureInfoCOC), intent(in) :: FracbyNode
      type(CSR), intent(in) :: FracbyNode

      double precision, allocatable, dimension(:, :, :), intent(in) :: PermCell
      double precision, allocatable, dimension(:), intent(in) :: PermFrac

      integer :: i, j, k, kp, num_node, num_parent, &
                 num_cell, num_face, comptCell, comptFrac
      logical :: cell_edge
      double precision :: meanDist, meanPerm, meanThickness, edge_length, wi, length
      double precision, dimension(3) :: xn1, xn2, xk
      double precision, allocatable, dimension(:) :: WI_global

      allocate (WI_global(NbNode))
      WI_global(:) = 0.d0
      !! MATRIX
      do i = 1, NbWell
         ! loop over the EDGES: (node,Parent(node)) for node=1,head_node-1
         do j = NodeDatabyWell%Pt(i) + 1, NodeDatabyWell%Pt(i + 1) - 1
            num_node = NodeDatabyWell%Num(j)
            num_parent = NodeDatabyWell%Val(j)%Parent
            xn1 = XNode(:, num_node)
            xn2 = XNode(:, num_parent)
            comptCell = 0
            meanDist = 0.d0
            meanPerm = 0.d0
            do k = CellbyNode%Pt(num_node) + 1, CellbyNode%Pt(num_node + 1) ! cells containing node
               cell_edge = .false.
               num_cell = CellbyNode%Num(k)
               do kp = CellbyNode%Pt(num_parent) + 1, CellbyNode%Pt(num_parent + 1) ! cells containing parent
                  if (num_cell == CellbyNode%Num(kp)) then ! cells containing Node and Parent, then containing the edge
                     cell_edge = .true.
                     exit
                  endif
               enddo
               if (cell_edge) then
                  comptCell = comptCell + 1
                  call element_center(XNode, NodebyCell, num_cell, xk)
                  call DistNodetoLine(xk, xn1, xn2, length) ! dist of cell center to the edge (num,parent)
                  meanDist = meanDist + length
                  meanPerm = meanPerm + PermCell(1, 1, num_cell) ! this formula is true for istropic permeability
               endif
            enddo ! loop over cells
            edge_length = dsqrt(dot_product(xn1 - xn2, xn1 - xn2)) ! length of edge
            call compute_peaceman_indices(meanDist, meanPerm, comptCell, edge_length, WellRadius(i), wi)
            ! contribution of the edge to each nodes: node and parent
            WI_global(num_node) = WI_global(num_node) + wi
            WI_global(num_parent) = WI_global(num_parent) + wi
         enddo ! loop over edges
      enddo ! loop over wells

      !! FRACTURES
      do i = 1, NbWell
         ! loop over the nodes
         do j = NodeDatabyWell%Pt(i) + 1, NodeDatabyWell%Pt(i + 1)
            num_node = NodeDatabyWell%Num(j)
            xn1(1:3) = XNode(1:3, num_node)
            meanDist = 0.d0
            meanPerm = 0.d0
            meanThickness = 0.d0
            comptFrac = 0
            ! loop over frac of node
            do k = FracbyNode%Pt(num_node) + 1, FracbyNode%Pt(num_node + 1)
               num_face = FracbyNode%Num(k)
               call element_center(XNode, NodebyFace, num_face, xk)
               length = dsqrt(dot_product(xk - xn1, xk - xn1)) ! dist of cell frac to the node
               meanDist = meanDist + length
               meanPerm = meanPerm + PermFrac(num_face) ! this formula is true for isotropic permeability
               meanThickness = meanThickness + Thickness
               comptFrac = comptFrac + 1
            enddo
            ! there is at least one frac in this node, compute wi
            if (comptFrac > 0) then
               meanThickness = meanThickness / dble( comptFrac )
               !Compute peaceman-indices, but for fractures  multiplying by a factor of 2 is necessary 
               call compute_peaceman_indices(meanDist, meanPerm, comptFrac, 2*meanThickness, WellRadius(i), wi)
               ! contribution of the fracs to the node
                WI_global(num_node) = WI_global(num_node) + wi
                              
            endif
         enddo
      enddo

      ! FIXME: unused Fourier Well Index
      NodeDatabyWell%Val(:)%WIF = 0.d0

      ! loop over all nodes of all wells
      do j = 1, NodeDatabyWell%Pt(NbWell + 1)
         num_node = NodeDatabyWell%Num(j)
         NodeDatabyWell%Val(j)%WID = WI_global(num_node)
      enddo

      deallocate (WI_global)

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
      if (length < 1.0d-5) then
         write (*, *) "WARNING - pb with DistNodetoLine, line(BC) is not a line because B=C"
         write (*, *) "          A:", NodeA, "(cell center ?)"
         write (*, *) "          B:", NodeB
         write (*, *) "          C:", NodeA
         length = dsqrt(dot_product(vectBA, vectBA))
         return
      endif

      ! cross product(BA, BC)
      cross_product(1) = vectBA(2)*vectBC(3) - vectBA(3)*vectBC(2)
      cross_product(2) = vectBA(3)*vectBC(1) - vectBA(1)*vectBC(3)
      cross_product(3) = vectBA(1)*vectBC(2) - vectBA(2)*vectBC(1)

      ! dist(A,BC) = norm( cross product(BA, BC) ) / norm(BC)
      length = dsqrt(dot_product(cross_product, cross_product)/length)

   end subroutine DistNodetoLine

   !> \brief Copy CSRDataWell1 to CSRDataWell2
   subroutine DefWell_csrdatawellcopy(CSR1, CSR2)

      type(TYPE_CSRDataNodeWell), intent(in)  :: CSR1
      type(TYPE_CSRDataNodeWell), intent(out) :: CSR2

      integer :: Nnz, i

      CSR2%Nb = CSR1%Nb

      allocate (CSR2%Pt(CSR1%Nb + 1))
      CSR2%Pt(:) = CSR1%Pt(:)

      Nnz = CSR1%Pt(CSR1%Nb + 1)

      allocate (CSR2%Num(Nnz))
      do i = 1, Nnz
         CSR2%Num(i) = CSR1%Num(i)
      end do

      allocate (CSR2%Val(Nnz))
      do i = 1, Nnz
         CSR2%Val(i) = CSR1%Val(i)
      end do

   end subroutine DefWell_csrdatawellcopy

   !> \brief Deallocate CSRDataWell (\%Pt, \%Num and \%Val)
   subroutine DefWell_deallocCSRDataWell(CSR1)

      type(TYPE_CSRDataNodeWell), intent(inout) :: CSR1

      deallocate (CSR1%Pt)

      if (associated(CSR1%Num)) then
         deallocate (CSR1%Num)
      end if

      if (associated(CSR1%Val)) then
         deallocate (CSR1%Val)
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

      x2%Id = x1%Id
      x2%Radius = x1%Radius
      x2%PressionMax = x1%PressionMax
      x2%PressionMin = x1%PressionMin
      x2%ImposedFlowrate = x1%ImposedFlowrate
      x2%CompTotal = x1%CompTotal
      x2%InjectionTemperature = x1%InjectionTemperature
      x2%actual_mass_flowrate = x1%actual_mass_flowrate
      x2%actual_energy_flowrate = x1%actual_energy_flowrate
      x2%actual_pressure = x1%actual_pressure
      x2%actual_temperature = x1%actual_temperature
      x2%IndWell = x1%IndWell

   end subroutine assign_DataWell_equal

   subroutine DefWell_mpi_register_well_data_description(mpi_id)

      integer, intent(out) :: mpi_id

      integer, parameter :: count = 12
      integer :: blocklengths(count)
      integer(kind=MPI_ADDRESS_KIND) :: begin, offset, displacements(count)
      integer :: types(count)
      type(WellData_type) :: dummy
      integer :: Ierr

      call MPI_Get_address(dummy, begin, Ierr)
      call MPI_Get_address(dummy%Id, offset, Ierr)
      displacements(1) = offset - begin
      call MPI_Get_address(dummy%Radius, offset, Ierr)
      displacements(2) = offset - begin
      call MPI_Get_address(dummy%PressionMax, offset, Ierr)
      displacements(3) = offset - begin
      call MPI_Get_address(dummy%PressionMin, offset, Ierr)
      displacements(4) = offset - begin
      call MPI_Get_address(dummy%ImposedFlowrate, offset, Ierr)
      displacements(5) = offset - begin
      call MPI_Get_address(dummy%CompTotal, offset, Ierr)
      displacements(6) = offset - begin
      call MPI_Get_address(dummy%InjectionTemperature, offset, Ierr)
      displacements(7) = offset - begin
      call MPI_Get_address(dummy%actual_mass_flowrate, offset, Ierr)
      displacements(8) = offset - begin
      call MPI_Get_address(dummy%actual_energy_flowrate, offset, Ierr)
      displacements(9) = offset - begin
      call MPI_Get_address(dummy%actual_pressure, offset, Ierr)
      displacements(10) = offset - begin
      call MPI_Get_address(dummy%actual_temperature, offset, Ierr)
      displacements(11) = offset - begin
      call MPI_Get_address(dummy%IndWell, offset, Ierr)
      displacements(12) = offset - begin

      types(:) = MPI_DOUBLE
      types(1) = MPI_INT
      types(11) = MPI_CHARACTER

      blocklengths(:) = 1    
      blocklengths(6) = NbComp

      call MPI_Type_Create_Struct(count, blocklengths, displacements, types, mpi_id, Ierr)
      if (Ierr /= 0) call CommonMPI_abort('Couldt not create well data MPI structure.')
      call MPI_Type_commit(mpi_id, Ierr)
      if (Ierr /= 0) call CommonMPI_abort('Couldt not commit well data MPI structure.')

   end subroutine DefWell_mpi_register_well_data_description

end module DefWell
