!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module NeumannContribution

   use iso_c_binding, only: c_int, c_double, c_bool
   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD

   use DefModel, only: NbComp
   use Physics, only: Thickness

   use MeshSchema, only: &
      XNodeLocal, XFaceLocal, NodebyFaceLocal, &
      NbNodeLocal_Ncpus, IdNodeLocal, &
      MeshSchema_local_face_surface, MeshSchema_triangle_area

   implicit none

   !> Unknown for Degree Of Freedom (including thermal). DOF can be Cell, Fracture Face or Node.
   type, bind(C) :: TYPE_NeumannBC
      real(c_double) :: molar_flux(NbComp) !< component molar flux
#ifdef _THERMIQUE_
      real(c_double) :: heat_flux !< heat flux
      logical(c_bool) :: compute_heat_flux
      real(c_double) :: nz
#endif
   end type TYPE_NeumannBC

   TYPE(TYPE_NeumannBC), allocatable, dimension(:), target, public :: &
      NodeNeumannBC !< Neumann contributions (size NbNodeLocal)

   public :: &
      NeumannContribution_allocate, &
      NeumannContribution_free, &
      NeumannContribution_clear, &
      NeumannContribution_set_fracture_edge_with_contribution, &
      NeumannContribution_set_fracture_edges_with_contribution, &
      NeumannContribution_set_fracture_edges, &
      NeumannContribution_set_face_with_contribution, &
      NeumannContribution_set_faces_with_contribution, &
      NeumannContribution_set_face

contains

   subroutine NeumannContribution_clear() &
      bind(C, name="clear_all_neumann_contributions")

      integer :: k

      do k = 1, size(NodeNeumannBC)
         NodeNeumannBC(k)%molar_flux(:) = 0.d0
#ifdef _THERMIQUE_
         NodeNeumannBC(k)%heat_flux = 0.d0
         NodeNeumannBC(k)%compute_heat_flux = .false.
         NodeNeumannBC(k)%nz = 0.d0
#endif
      end do

   end subroutine NeumannContribution_clear

   subroutine NeumannContribution_set_fracture_edge_with_contribution(edge, fluxes) &
      bind(C, name="set_fracture_edge_with_neumann_contribution")

      integer(c_int) :: edge(2)
      type(TYPE_NeumannBC), intent(in) :: fluxes

      integer :: k, s
      integer :: Ierr, errcode ! used for MPI_Abort
      double precision :: half_length

      if (IdNodeLocal(edge(1))%Frac /= 'y' .or. IdNodeLocal(edge(2))%Frac /= 'y') then
         print *, 'ERROR: edge', edge, 'does not seem to be a fracture edge on proc', commRank + 1
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      half_length = 0.5*norm2(XNodeLocal(:, edge(2)) - XNodeLocal(:, edge(1)))
      do k = 1, 2
         s = edge(k)
         ! FIXME: This should take into account node fractions (and variable thicknesses...)
         NodeNeumannBC(s)%molar_flux = NodeNeumannBC(s)%molar_flux &
                                       + half_length*Thickness*fluxes%molar_flux
#ifdef _THERMIQUE_
         NodeNeumannBC(s)%heat_flux = NodeNeumannBC(s)%heat_flux &
                                      + half_length*Thickness*fluxes%heat_flux
         if (NodeNeumannBC(s)%compute_heat_flux .and. .not. fluxes%compute_heat_flux) &
            write (*, *) "WARNING Unsetting heat flux computation"
         NodeNeumannBC(s)%compute_heat_flux = fluxes%compute_heat_flux
         ! FIXME: only valid for planar faces
         NodeNeumannBC(s)%nz = fluxes%nz
#endif
      end do

   end subroutine NeumannContribution_set_fracture_edge_with_contribution

   subroutine NeumannContribution_set_fracture_edges_with_contribution(nbedges, edges, fluxes) &
      bind(C, name="set_fracture_edges_with_neumann_contribution")

      integer(c_int), value :: nbedges
      integer(c_int) :: edges(2, nbedges)
      type(TYPE_NeumannBC), intent(in) :: fluxes

      integer :: k

      do k = 1, nbedges
         call NeumannContribution_set_fracture_edge_with_contribution(edges(:, k), fluxes)
      end do

   end subroutine NeumannContribution_set_fracture_edges_with_contribution

   subroutine NeumannContribution_set_fracture_edges(nbcont, edges, fluxes) &
      bind(C, name="set_fracture_edge_neumann_contributions")

      integer(c_int), value :: nbcont
      integer(c_int) :: edges(2, nbcont)
      type(TYPE_NeumannBC), intent(in) :: fluxes(nbcont)

      integer :: k

      do k = 1, nbcont
         call NeumannContribution_set_fracture_edge_with_contribution(edges(:, k), fluxes(k))
      end do

   end subroutine NeumannContribution_set_fracture_edges

   subroutine NeumannContribution_set_face_with_contribution(face, fluxes) &
      bind(C, name="set_face_with_neumann_contribution")

      integer(c_int), value :: face
      type(TYPE_NeumannBC), intent(in) :: fluxes

      integer :: k, nbnodes
      integer, allocatable :: nodes(:)
      double precision :: barycenter(3)
      double precision :: alpha

      barycenter(:) = XFaceLocal(:, face)
      nbnodes = NodebyFaceLocal%Pt(face + 1) - NodebyFaceLocal%Pt(face)
      allocate (nodes(nbnodes + 2))
      nodes(1:nbnodes) = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(face) + 1:NodebyFaceLocal%Pt(face + 1))
      nodes(nbnodes + 1:nbnodes + 2) = nodes(1:2)
      do k = 2, nbnodes + 1
         alpha = MeshSchema_local_face_surface(face)/nbnodes
         alpha = alpha + MeshSchema_triangle_area(XNodeLocal(:, nodes(k - 1)), XNodeLocal(:, nodes(k)), barycenter)
         alpha = alpha + MeshSchema_triangle_area(XNodeLocal(:, nodes(k)), XNodeLocal(:, nodes(k + 1)), barycenter)
         alpha = alpha/3.d0
         NodeNeumannBC(nodes(k))%molar_flux = NodeNeumannBC(nodes(k))%molar_flux + alpha*fluxes%molar_flux
#ifdef _THERMIQUE_
         NodeNeumannBC(nodes(k))%heat_flux = NodeNeumannBC(nodes(k))%heat_flux + alpha*fluxes%heat_flux
         if (NodeNeumannBC(nodes(k))%compute_heat_flux .and. .not. fluxes%compute_heat_flux) &
            write (*, *) "WARNING Unsetting heat flux computation"
         NodeNeumannBC(nodes(k))%compute_heat_flux = fluxes%compute_heat_flux
         ! FIXME: only valid for planar faces
         NodeNeumannBC(nodes(k))%nz = fluxes%nz
#endif
      end do
      deallocate (nodes)

   end subroutine NeumannContribution_set_face_with_contribution

   subroutine NeumannContribution_set_faces_with_contribution(nbfaces, faces, fluxes) &
      bind(C, name="set_faces_with_neumann_contribution")

      integer(c_int), value :: nbfaces
      integer(c_int) :: faces(nbfaces)
      type(TYPE_NeumannBC), intent(in) :: fluxes

      integer :: k

      do k = 1, nbfaces
         call NeumannContribution_set_face_with_contribution(faces(k), fluxes)
      end do

   end subroutine NeumannContribution_set_faces_with_contribution

   subroutine NeumannContribution_set_face(nbcont, faces, fluxes) &
      bind(C, name="set_face_neumann_contributions")

      integer(c_int), value :: nbcont
      integer(c_int) :: faces(nbcont)
      type(TYPE_NeumannBC), intent(in) :: fluxes(nbcont)

      integer :: k

      do k = 1, nbcont
         call NeumannContribution_set_face_with_contribution(faces(k), fluxes(k))
      end do

   end subroutine NeumannContribution_set_face

   subroutine NeumannContribution_allocate

      allocate (NodeNeumannBC(NbNodeLocal_Ncpus(commRank + 1)))
      call NeumannContribution_clear

   end subroutine NeumannContribution_allocate

   !> \brief Deallocate unknowns vectors
   subroutine NeumannContribution_free

      deallocate (NodeNeumannBC)

   end subroutine NeumannContribution_free

end module NeumannContribution
