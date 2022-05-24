module MSWellsData
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use DefModel
   use mpi
   use CommonMPI
   use IncCVReservoir
   use MeshSchema
   use IncCVMSWells
#else
   use iso_c_binding, only: c_double
   use DefModel, only: NbPhase
#ifdef ComPASS_DIPHASIC_CONTEXT
   use DefModel, only: DIPHASIC_CONTEXT
#endif
   use CommonMPI, only: commRank
   use MeshSchema, only: &
      XNodeLocal, NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus
#endif

   implicit none
   !Edge data
   type, bind(C) :: TYPE_MSWellDataEdge
      real(c_double) :: cross_section
      real(c_double) :: orientation
   end type TYPE_MSWellDataEdge

   type, bind(C) ::TYPE_MSWellDataNode

      !Node data
      type(TYPE_MSWellDataEdge)   :: edata    !node son has the edge data
      real(c_double)              :: vol      !node vol: average of  the edge volumes which touch the node
#ifdef _THERMIQUE_
      real(c_double)              :: thcond(NbPhase)   !Thermal conductivity
#endif
   end type TYPE_MSWellDataNode

   ! like CSR format without stocking %Nb and %Pt (idem that NodebyWellLocal)
   TYPE(TYPE_MSWellDataNode), allocatable, dimension(:), target, public :: &
      MSWellDataNodebyMSWell  !<  Additional Data for MSWellNode

   public :: &
      MSWellsData_allocate, &
      MSWellsData_free
contains

   !> \brief Allocate well unknowns vectors

   subroutine MSWellsData_allocate()

      integer :: Nb, Nnz

      Nb = NodebyMSWellLocal%Nb
      Nnz = NodebyMSWellLocal%Pt(Nb + 1)
      allocate (MSWellDataNodebyMSWell(Nnz))
   end subroutine MSWellsData_allocate

   !> \brief Deallocate well unknowns vectors
   subroutine MSWellsData_free()

      deallocate (MSWellDataNodebyMSWell)

   end subroutine MSWellsData_free

   subroutine MSWellsData_init() &
      bind(C, name="MSWellsData_init")

      double precision :: well_radius
      integer :: s, k, nbwells, num_s, num_sp, wnum_sp
      double precision :: xk_s, xk_sp, yk_s, yk_sp, zk_s, zk_sp, dz, ds
      double precision :: theta, factor, vol_half_edge
      double precision, parameter :: pi = datan(1.0d0)*4.0d0
#ifdef _THERMIQUE_
      double precision, parameter :: user_thcond = 2.0 !Thermal conductivity
#endif

      nbwells = NbMSWellLocal_Ncpus(commRank + 1)

      ! For producers and injectors
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Init node data: volume
      MSWellDataNodebyMSWell(:)%vol = 0.d0 !Init vol (see below)

      do k = 1, nbwells

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Init edge data: cross section and oritentation
         !User setting
         well_radius = DataMSWellLocal(k)%Radius
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! looping for edges, i.e., from  queue to the node before the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            MSWellDataNodebyMSWell(s)%edata%cross_section = pi*well_radius**2

            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the parent in the reservoir

#ifdef _THERMIQUE_
            MSWellDataNodebyMSWell(s)%thcond = user_thcond      !Thermal conductivity
#endif
            xk_s = XNodeLocal(1, num_s)    ! x-cordinate of node s
            xk_sp = XNodeLocal(1, num_sp)  ! x-cordinate of parent of s
            yk_s = XNodeLocal(2, num_s)    ! y-cordinate of node s
            yk_sp = XNodeLocal(2, num_sp)  ! y-cordinate of parent of s
            zk_s = XNodeLocal(3, num_s)    ! z-cordinate of node s
            zk_sp = XNodeLocal(3, num_sp)  ! z-cordinate of parent of s

            ds = dsqrt((xk_s - xk_sp)**2 + (yk_s - yk_sp)**2 + (zk_s - zk_sp)**2)
            dz = dabs(zk_sp - zk_s)
            theta = dacos(dz/ds)
            factor = dsqrt(dcos(theta))*(1 + dsin(theta))**2

            if (zk_sp .gt. zk_s) then
               MSWellDataNodebyMSWell(s)%edata%orientation = factor
            else
               MSWellDataNodebyMSWell(s)%edata%orientation = -factor
            end if

         end do

         !For the head we take the previous one
         s = NodebyMSWellLocal%Pt(k + 1)!head
         MSWellDataNodebyMSWell(s)%edata%cross_section = MSWellDataNodebyMSWell(s - 1)%edata%cross_section

#ifdef _THERMIQUE_
         MSWellDataNodebyMSWell(s)%thcond = user_thcond      !Thermal conductivity
#endif

         ! looping for edges, i.e., from  queue to the node before the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1) - 1

            num_s = NodebyMSWellLocal%Num(s)    ! index of s in the reservoir (local mesh)
            num_sp = NodeDatabyMSWellLocal%Val(s)%Parent ! index of the parent in the reservoir
            wnum_sp = NodeDatabyMSWellLocal%Val(s)%PtParent ! index of the parent in the well (global numbering)

            xk_s = XNodeLocal(1, num_s)    ! x-cordinate of node s
            xk_sp = XNodeLocal(1, num_sp)  ! x-cordinate of parent of s
            yk_s = XNodeLocal(2, num_s)    ! y-cordinate of node s
            yk_sp = XNodeLocal(2, num_sp)  ! y-cordinate of parent of s
            zk_s = XNodeLocal(3, num_s)    ! z-cordinate of node s
            zk_sp = XNodeLocal(3, num_sp)  ! z-cordinate of parent of s

            ds = dsqrt((xk_s - xk_sp)**2 + (yk_s - yk_sp)**2 + (zk_s - zk_sp)**2)

            vol_half_edge = 0.5*ds*MSWellDataNodebyMSWell(s)%edata%cross_section
            MSWellDataNodebyMSWell(s)%vol = MSWellDataNodebyMSWell(s)%vol + vol_half_edge
            MSWellDataNodebyMSWell(wnum_sp)%vol = MSWellDataNodebyMSWell(wnum_sp)%vol + vol_half_edge
         end do

      end do!k

   end subroutine MSWellsData_init

!
end module MSWellsData
