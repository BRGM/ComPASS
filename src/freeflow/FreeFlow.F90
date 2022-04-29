!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module FreeFlow

   use iso_c_binding, only: c_int, c_double, c_size_t, c_f_pointer, c_loc, c_null_ptr
   use CommonMPI, only: CommonMPI_abort
   use InteroperabilityStructures, only: cpp_array_wrapper
   use MeshSchema, only: &
      IdNodeLocal, AtmState, &
      IsFreeflowNode, FreeflowFaces, &
      MeshSchema_local_face_surface, &
      SurfFreeFlowLocal, &
      number_of_nodes, &
      NodebyFaceLocal
   use DefModel, only: LIQUID_PHASE
#ifdef ComPASS_WITH_diphasic_PHYSICS
   use DefModel, only: GAS_PHASE
#endif
   use Physics, only: atm_pressure, atm_comp, atm_temperature, rain_temperature, rain_flux, Hm, HT

contains

   subroutine FreeFlow_reset_faces() &
      bind(C, name="clear_freeflow_faces")
      integer(c_size_t) ::  nn

      nn = number_of_nodes()

      if (.not. allocated(IsFreeflowNode)) then
         allocate (IsFreeflowNode(nn))
         if (allocated(FreeflowFaces)) & ! is allocated in FreeFlow_set_faces
            call CommonMPI_abort("FreeFlow_set_freeflow_faces: FreeflowFaces should not be already allocated.")
         if (allocated(SurfFreeFlowLocal)) &
            call CommonMPI_abort("FreeFlow_set_freeflow_faces: SurfFreeFlowLocal should not be already allocated.")
         allocate (SurfFreeFlowLocal(nn))
         if (allocated(AtmState)) &
            call CommonMPI_abort("FreeFlow_set_freeflow_faces: AtmState should not be already allocated.")
         allocate (AtmState(nn))
      endif

      if (size(IsFreeflowNode) /= nn) &
         call CommonMPI_abort("FreeFlow_set_freeflow_faces: inconsistent size for node flags.")
      if (size(SurfFreeFlowLocal) /= nn) &
         call CommonMPI_abort("FreeFlow_set_freeflow_faces: inconsistent size for freeflow surface array.")

      IsFreeflowNode = .false.
      SurfFreeFlowLocal = 0.d0
      if (allocated(FreeflowFaces)) deallocate (FreeflowFaces)

   end subroutine FreeFlow_reset_faces

   subroutine FreeFlow_set_area_distribution(faces, Freeflow_area_distribution)
      integer(c_int), intent(in) :: faces(:)
      real(c_double), dimension(:), allocatable, intent(inout) :: Freeflow_area_distribution
      integer(c_int) :: k, fk
      integer(c_int) :: s, p, nb_freeflow_nodes
      real(c_double) :: surface_fraction

      do k = 1, size(faces)
         fk = faces(k)
         nb_freeflow_nodes = 0
         do p = NodebyFaceLocal%Pt(fk) + 1, NodebyFaceLocal%Pt(fk + 1)
            s = NodebyFaceLocal%Num(p)
            if (IsFreeflowNode(s)) & ! in FF (not a Dirichlet node)
               nb_freeflow_nodes = nb_freeflow_nodes + 1
         enddo
         surface_fraction = MeshSchema_local_face_surface(fk)/nb_freeflow_nodes
         do p = NodebyFaceLocal%Pt(fk) + 1, NodebyFaceLocal%Pt(fk + 1)
            s = NodebyFaceLocal%Num(p)
            if (IsFreeflowNode(s)) &
               Freeflow_area_distribution(s) = Freeflow_area_distribution(s) + surface_fraction
         enddo
      enddo

   end subroutine FreeFlow_set_area_distribution

   subroutine FreeFlow_set_faces(faces)
      integer(c_int), intent(in) :: faces(:)
      integer(c_int) :: k, fk
      integer(c_int) :: s, p

      ! store the freeflow faces index
#ifndef NDEBUG
      if (allocated(FreeflowFaces)) &
         call CommonMPI_abort("FreeFlow_set_faces: FreeflowFaces should not be already allocated.")
#endif
      allocate (FreeflowFaces(size(faces)))
      FreeflowFaces = faces

      do k = 1, size(faces)
         fk = faces(k)
         do p = NodebyFaceLocal%Pt(fk) + 1, NodebyFaceLocal%Pt(fk + 1)
            s = NodebyFaceLocal%Num(p)
            if (IdNodeLocal(s)%T .ne. "d" .AND. IdNodeLocal(s)%P .ne. "d") then
               IsFreeflowNode(s) = .true.
               AtmState(s)%Pressure = atm_pressure
               AtmState(s)%Temperature(LIQUID_PHASE) = rain_temperature
#ifdef ComPASS_WITH_diphasic_PHYSICS
               AtmState(s)%Temperature(GAS_PHASE) = atm_temperature
#endif
               AtmState(s)%Comp = atm_comp
               AtmState(s)%Imposed_flux = rain_flux
               AtmState(s)%Hm = Hm
               AtmState(s)%HT = HT
            endif
         enddo
      enddo

      ! fill SurfFreeFlowLocal with the distribution of the face area
      ! over the freeflow nodes (no Dirichlet node)
      call FreeFlow_set_area_distribution(faces, SurfFreeFlowLocal)

   end subroutine FreeFlow_set_faces

   subroutine FreeFlow_reset_freeflow_nodes() &
      bind(C, name="reset_freeflow_nodes")
      integer :: k, fk
      integer :: s, p

      ! reset the old values (the Dirichlet nodes may have changed,
      ! the freeflow faces no, so use FreeflowFaces)
      ! Do not call FreeFlow_reset_faces to keep FreeflowFaces
      IsFreeflowNode = .false.
      SurfFreeFlowLocal = 0.d0

      ! identify again the FF nodes
      do k = 1, size(FreeflowFaces)
         fk = FreeflowFaces(k)
         do p = NodebyFaceLocal%Pt(fk) + 1, NodebyFaceLocal%Pt(fk + 1)
            s = NodebyFaceLocal%Num(p)
            if (IdNodeLocal(s)%T .ne. "d" .AND. IdNodeLocal(s)%P .ne. "d") then
               IsFreeflowNode(s) = .true.
            endif
         enddo
      enddo

      ! fill SurfFreeFlowLocal with the distribution of the FF face area
      ! over the Freeflow nodes
      call FreeFlow_set_area_distribution(FreeflowFaces, SurfFreeFlowLocal)

   end subroutine FreeFlow_reset_freeflow_nodes

   subroutine FreeFlow_set_faces_C(faces_wrapper) &
      bind(C, name="set_freeflow_faces")
      type(cpp_array_wrapper), intent(in) :: faces_wrapper
      integer(c_int), pointer :: faces(:)

      call c_f_pointer(faces_wrapper%p, faces, [faces_wrapper%n])
      call FreeFlow_set_faces(faces)

   end subroutine FreeFlow_set_faces_C

   subroutine set_atm_temperature(T) &
      bind(C, name="set_atm_temperature")
      real(c_double), value, intent(in) :: T
      integer(c_size_t) ::  nn, s

      ! if not allocated, AtmState will later be initialized with atm_temperature
      atm_temperature = T
      if (allocated(AtmState)) then
         nn = number_of_nodes()
#ifdef ComPASS_WITH_diphasic_PHYSICS
         do s = 1, nn
            if (IsFreeflowNode(s)) then
               AtmState(s)%Temperature(GAS_PHASE) = atm_temperature
            endif
         enddo
#endif
      endif
#ifndef NDEBUG
      write (*, *) " ****************** "
      write (*, *) &
         "WARNING: rain_temperature is not modified, use set_rain_temperature"
      write (*, *) " ****************** "
#endif

   end subroutine set_atm_temperature

   subroutine set_rain_temperature(T) &
      bind(C, name="set_rain_temperature")
      real(c_double), value, intent(in) :: T
      integer(c_size_t) ::  nn, s

      ! if not allocated, AtmState will later be initialized with rain_temperature
      rain_temperature = T
      if (allocated(AtmState)) then
         nn = number_of_nodes()
         do s = 1, nn
            if (IsFreeflowNode(s)) then
               AtmState(s)%Temperature(LIQUID_PHASE) = rain_temperature
            endif
         enddo
      endif

   end subroutine set_rain_temperature

   subroutine set_atm_pressure(p) &
      bind(C, name="set_atm_pressure")
      real(c_double), value, intent(in) :: p
      integer(c_size_t) ::  nn, s

      ! if not allocated, AtmState will later be initialized with atm_pressure
      atm_pressure = p
      if (allocated(AtmState)) then
         nn = number_of_nodes()
         do s = 1, nn
            if (IsFreeflowNode(s)) then
               AtmState(s)%Pressure = atm_pressure
            endif
         enddo
      endif

   end subroutine set_atm_pressure

   subroutine set_atm_rain_flux(q_rain) &
      bind(C, name="set_atm_rain_flux")
      real(c_double), value, intent(in) :: q_rain
      integer(c_size_t) ::  nn, s

      ! if not allocated, AtmState will later be initialized with rain_flux
      rain_flux(LIQUID_PHASE) = q_rain
      if (allocated(AtmState)) then
         nn = number_of_nodes()
         do s = 1, nn
            if (IsFreeflowNode(s)) then
               AtmState(s)%Imposed_flux(LIQUID_PHASE) = q_rain
            endif
         enddo
      endif

   end subroutine set_atm_rain_flux

   subroutine retrieve_freeflow_nodes_mask(cpp_array) &
      bind(C, name="retrieve_freeflow_nodes_mask")

      type(cpp_array_wrapper), intent(inout) :: cpp_array

      cpp_array%p = c_loc(IsFreeflowNode(1))
      cpp_array%n = size(IsFreeflowNode)

   end subroutine retrieve_freeflow_nodes_mask

   subroutine retrieve_freeflow_nodes_area(cpp_array) &
      bind(C, name="retrieve_freeflow_nodes_area")

      type(cpp_array_wrapper), intent(inout) :: cpp_array

      cpp_array%p = c_loc(SurfFreeFlowLocal(1))
      cpp_array%n = size(SurfFreeFlowLocal)

   end subroutine retrieve_freeflow_nodes_area

   subroutine retrieve_freeflow_node_states(cpp_array) &
      bind(C, name="retrieve_freeflow_node_states")

      type(cpp_array_wrapper), intent(out) :: cpp_array
      integer(c_size_t) :: n

      if (.not. allocated(AtmState)) then
         cpp_array%p = c_null_ptr
         cpp_array%n = 0
      else
         n = size(AtmState)
         cpp_array%n = n
         if (n == 0) then
#ifdef TRACK_ZERO_SIZE_ARRAY
            ! FIXME: Remove comment
            write (*, *) '!!!!!!!!!!!!!!!!!!!!!!! Zero size array'
#endif
            cpp_array%p = c_null_ptr
         else
            cpp_array%p = c_loc(AtmState(1))
         end if
      end if

   end subroutine retrieve_freeflow_node_states

end module FreeFlow
