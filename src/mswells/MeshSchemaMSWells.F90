module MeshSchemaMSWells
#ifdef _COMPASS_FORTRAN_DO_NOT_USE_ONLY_
   use iso_c_binding
   use mpi
   use CommonMPI
   use CommonType
   use MeshSchema
#else
   use iso_c_binding, only: c_double, c_bool

   use mpi, only: MPI_Abort
   use CommonMPI, only: commRank, ComPASS_COMM_WORLD, CommonMPI_abort

   use CommonType, only: &
      CSRArray2dble, CSR, CommonType_deallocCSR

   use MeshSchema, only: &
      XNodeLocal, &
      IdNodeLocal, &
      NodebyMSWellLocal, NodeDatabyMSWellLocal, &
      DataMSWellLocal, NbMSWellLocal_Ncpus, &
      NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
      NbMSWellNodeOwn_Ncpus, NbMSWellNodeLocal_Ncpus, &
      MSWellNodePart
#endif
   implicit none
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !This modules creates a new variable unknown ordering for each node in each mswells to be used
   !in the assembly of the linear system, i.e., Jacobian matrix and rhs.
   !Variables indexing stuff
   integer :: &
      NbMSWellNodeOwn, &     ! number of nodes own which come from local mswells
      NbMSWellNodeLocal, &   ! number of nodes (own and ghosts) which come from local mswells
      NbMSWellLocal         ! number of local mswells
   type(CSR), protected :: &
      IncIdxNodebyMSWellLocal  !New unkowns indexing  for local mswells
   !Num: contains the new indexing
   !Val: contains the reservoir node which is paired with the mswell_node.

   public :: &
      MeshSchemaMSWells_make, &    !< init variables indexing
      MeshSchemaMSWells_free

contains

   subroutine MeshSchemaMSWells_make

      integer :: s, k, num_s, Nnz

      !Let us set the new  unknowns ordering for producers to be used for the Jacobian
      NbMSWellLocal = NbMSWellLocal_Ncpus(commRank + 1)    !nblocal mswells
      NbMSWellNodeLocal = NbMSWellNodeLocal_Ncpus(commRank + 1)  ! number of nodes (own and ghosts) which come from local mswells

      !Mapping from mswell node idx to new unkown idx
      IncIdxNodebyMSWellLocal%Nb = NodebyMSWellLocal%Nb !Number of local mswells
      Nnz = NodebyMSWellLocal%Pt(IncIdxNodebyMSWellLocal%Nb + 1)
      allocate (IncIdxNodebyMSWellLocal%Pt(IncIdxNodebyMSWellLocal%Nb + 1))
      allocate (IncIdxNodebyMSWellLocal%Num(Nnz))
      allocate (IncIdxNodebyMSWellLocal%Val(Nnz))
      IncIdxNodebyMSWellLocal%Pt(:) = NodebyMSWellLocal%Pt(:)
      NbMSWellNodeOwn = NbMSWellNodeOwn_Ncpus(commRank + 1)! number of nodes own which come from local mswells

#ifndef NDEBUG
      if (NbMSWellNodeLocal /= MSWellNodePart%Pt(MSWellNodePart%Nb + 1)) &
         call CommonMPI_abort("Number of mswell nodes are not equal in module MeshSchemaMSWells")
#endif

      !We first enumerate the mswell nodes own using indices from 1...NbMSWellNodeOwn
      !Then  the we enumerate the other mswell nodes using indices from NbMSWellNodeOwn+1...NbMSWellNodeLocal,
      !but  respecting the procs index-order.
      !This is almost  done in LocalMesh%MSWellNodebyProc(n_proc)
      !which has been stored in MeshSchema%MSWellNodePart

      do s = 1, MSWellNodePart%Pt(MSWellNodePart%Nb + 1)
         IncIdxNodebyMSWellLocal%Num(MSWellNodePart%Num(s)) = s
      enddo

      !Setting the reservoir_node which is paired with mswell_node
      do k = 1, NbMSWellLocal
         ! looping from  queue to  the head
         do s = NodebyMSWellLocal%Pt(k) + 1, NodebyMSWellLocal%Pt(k + 1)
            num_s = NodebyMSWellLocal%Num(s)
            IncIdxNodebyMSWellLocal%Val(s) = num_s
         end do
      end do

   end subroutine MeshSchemaMSWells_make

   subroutine MeshSchemaMSWells_free

      call CommonType_deallocCSR(IncIdxNodebyMSWellLocal)

   end subroutine MeshSchemaMSWells_free

end module MeshSchemaMSWells
