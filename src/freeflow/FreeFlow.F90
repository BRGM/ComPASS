!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module FreeFlow

   use iso_c_binding, only: c_int, c_double, c_size_t, c_f_pointer
   use CommonMPI, only: CommonMPI_abort
   use InteroperabilityStructures, only: cpp_array_wrapper
   use MeshSchema, only: &
      IdFFNodeLocal, &
      MeshSchema_local_face_surface, &
      SurfFreeFlowLocal, &
      number_of_nodes, &
      NodebyFaceLocal

contains

   subroutine FreeFlow_reset_faces() &
      bind(C, name="clear_freeflow_faces")
      integer(c_size_t) ::  nn

      nn = number_of_nodes()

      if (.not. allocated(IdFFNodeLocal)) then
         allocate (IdFFNodeLocal(nn))
         if (allocated(SurfFreeFlowLocal)) &
            call CommonMPI_abort("FreeFlow_set_freeflow_faces: SurfFreeFlowLocal should not be already allocated.")
         allocate (SurfFreeFlowLocal(nn))
      endif

      if (size(IdFFNodeLocal) /= nn) &
         call CommonMPI_abort("FreeFlow_set_freeflow_faces: inconsistent size for node flags.")
      if (size(SurfFreeFlowLocal) /= nn) &
         call CommonMPI_abort("FreeFlow_set_freeflow_faces: inconsistent size for freeflow surface array.")

      IdFFNodeLocal = .false.
      SurfFreeFlowLocal = 0.d0

   end subroutine FreeFlow_reset_faces

   subroutine FreeFlow_set_faces(faces)
      integer(c_int), intent(in) :: faces(:)
      integer(c_int) :: k, fk
      real(c_double) :: surface_fraction
      integer(c_int) :: s, p, nb_facenodes

      do k = 1, size(faces)
         fk = faces(k)
         nb_facenodes = NodebyFaceLocal%Pt(fk + 1) - NodebyFaceLocal%Pt(fk)
         surface_fraction = MeshSchema_local_face_surface(fk)/nb_facenodes
         do p = NodebyFaceLocal%Pt(fk) + 1, NodebyFaceLocal%Pt(fk + 1)
            s = NodebyFaceLocal%Num(p)
            IdFFNodeLocal(s) = .true.
            SurfFreeFlowLocal(s) = SurfFreeFlowLocal(s) + surface_fraction
         enddo
      enddo

   end subroutine FreeFlow_set_faces

   subroutine FreeFlow_set_faces_C(faces_wrapper) &
      bind(C, name="set_freeflow_faces")
      type(cpp_array_wrapper), intent(in) :: faces_wrapper
      integer(c_int), pointer :: faces(:)

      call c_f_pointer(faces_wrapper%p, faces, [faces_wrapper%n])
      call FreeFlow_set_faces(faces)

   end subroutine FreeFlow_set_faces_C

end module FreeFlow
