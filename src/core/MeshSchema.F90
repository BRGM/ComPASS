module MeshSchema

  use CommonType
  use CommonMPI
  use LocalMesh
  use DefModel
  use DefWell
  use mpi

  implicit none

  ! 1. Mesh Info
  integer, allocatable, dimension(:), protected :: &
       NbCellLocal_Ncpus, NbCellOwn_Ncpus, &
       NbFaceLocal_Ncpus, NbFaceOwn_Ncpus, &
       NbNodeLocal_Ncpus, NbNodeOwn_Ncpus, &
       NbFracLocal_Ncpus, NbFracOwn_Ncpus

  double precision, protected :: &
       meshSize_xmax, & !< Global size of mesh xmax, known by all CPUs
       meshSize_xmin, & !< Global size of mesh xmin, known by all CPUs
       meshSize_ymax, & !< Global size of mesh ymax, known by all CPUs
       meshSize_ymin, & !< Global size of mesh ymin, known by all CPUs
       meshSize_zmax, & !< Global size of mesh zmax, known by all CPUs
       meshSize_zmin    !< Global size of mesh zmin, known by all CPUs

  ! 2. Mesh and connectivity
  !    Used to (1) calcul Transmissivities: TkLocal and TkFracLocal
  !            (2) Assembly
  type(CSR), protected :: & 
       NodebyNodeOwn, &   ! (1)
       FracbyNodeOwn, &   ! (1)
       CellbyNodeOwn, &   ! (1)
                                !
       NodebyFracOwn, &   ! (1)
       CellbyFracOwn, &   ! (1)
       FracbyFracOwn, &   ! (1)
                                ! 
       FacebyCellLocal, & ! (1,2)
       FracbyCellLocal, & ! (1,2)
       NodebyCellLocal, & ! (1,2)
                                !
       NodebyFaceLocal    ! (1,2)

  ! Number of Edges by Well
  integer, allocatable, dimension(:), protected :: &
       NbEdgebyWellInjLocal, &
       NbEdgebyWellProdLocal

  ! Num (local) of Nodes by Edge by Well local (own+ghost)
  integer, allocatable, dimension(:,:,:), protected :: &
       NumNodebyEdgebyWellInjLocal, &
       NumNodebyEdgebyWellProdLocal

  ! 3. X node
  double precision, allocatable, dimension(:,:), target :: &
       XNodeLocal
  integer(c_int), allocatable, dimension(:), target :: &
       NodeFlagsLocal
  
  ! 4. IdCell/IdFace/IdNode
  integer, allocatable, dimension(:), protected :: &
       IdCellLocal, &
       IdFaceLocal

  type(Type_IdNode), allocatable, dimension(:), target :: &
       IdNodeLocal

  ! Well 
  integer, dimension(:), allocatable, protected :: &
       NbWellInjLocal_Ncpus, NbWellInjOwn_Ncpus, &
       NbWellProdLocal_Ncpus, NbWellProdOwn_Ncpus

  ! Well connectivity in local 
  type(CSR), protected :: &
       NodebyWellInjLocal, &
       NodebyWellProdLocal

  type(TYPE_DataWellInj), allocatable, dimension(:), public :: &
       DataWellInjLocal !< Data of injection well (Radius,...) for local well
  type(TYPE_DataWellProd), allocatable, dimension(:), public :: &
       DataWellProdLocal !< Data of production well (Radius,...) for local well

  ! Data in nodes of wells
  type(TYPE_CSRDataNodeWell), protected :: &
       NodeDatabyWellInjLocal, &
       NodeDatabyWellProdLocal

  !! The follwoing vectors are used for the strucutre of Jacobian
  type(CSR), protected :: &
       WellInjbyNodeOwn, &  ! numero (local) of well inj connected to this node own
       WellProdbyNodeOwn    ! numero (local) of well prod connected to this node own

  ! 5. Frac to Face, Face to Frac 
  integer, allocatable, dimension(:), protected :: &
       FracToFaceLocal, &
       FaceToFracLocal

  ! 6. NumNodebyProc, NumFracbyProc
  type(CSR), protected:: &
       NumNodebyProc,    &
       NumFracbyProc,    &
       NumWellInjbyProc, &
       NumWellProdbyProc

  ! 7. XCellLocal, XFaceLocal
  double precision, dimension(:,:), allocatable, public :: &
       XCellLocal,&     ! center of cell
       XFaceLocal       ! center of frac

  ! 8. VolCellLocal, SurfFracLocal
  double precision, dimension(:), allocatable, protected :: &
       VolCellLocal, &  ! vol of cell
       SurfFracLocal    ! surf of frac face

  ! 9. max number of nodes/frac in a cell 
  integer, protected :: & 
       NbNodeCellMax, &
       NbFracCellMax, &
       NbNodeFaceMax

  ! 10. Porosity
  double precision, dimension(:), allocatable, protected :: &
       PorositeCellLocal, &
       PorositeFracLocal
  ! permeability
  double precision, dimension(:,:,:), allocatable, public :: &
       PermCellLocal
  double precision, dimension(:), allocatable, public :: &
       PermFracLocal

  ! MPI TYPE for DataNodewell: MPI_DATANODEWELL
  integer, private :: MPI_DATANODEWELL

  ! tmp constant values
  integer, parameter, private :: &
       VALSIZE_ZERO = 0, &
       VALSIZE_NB   = 1, &
       VALSIZE_NNZ  = 2

  private :: &
       MeshSchema_csrsend, &   ! send csr
       MeshSchema_csrrecv, &   ! recv csr
       MeshSchema_csrdatawellsend, & ! send csrdatawell
       MeshSchema_csrdatawellrecv, & ! recv csrdatawell
                                !
       MeshSchema_sendrecv,                 &
       MeshSchema_NumNodebyEdgebyWellLocal, &
       MeshSchema_XCellLocal,               &
       MeshSchema_VolCellLocal,             &
       MeshSchema_XFaceLocal,               &
       MeshSchema_SurfFracLocal,            &
       MeshSchema_Surf12f,                  & ! used by MeshSchema_SurfFraclocal
       MeshSchema_NbNodeCellMax,            &
       MeshSchema_NbNodeFaceMax

  public :: &
       MeshSchema_make, &
       MeshSchema_Free

contains

  subroutine MeshSchema_make

    call MeshSchema_sendrecv

    ! List of Edge (local number of nodes) by local Well (own+ghost)
    call MeshSchema_NumNodebyEdgebyWellLocal

    ! XCellLocal and XFaceLocal
    call MeshSchema_XCellLocal
    call MeshSchema_XFaceLocal

    ! VolCellLocal and SurfFraclocal
    call MeshSchema_VolCellLocal
    call MeshSchema_SurfFracLocal

    call MeshSchema_NbNodeCellMax ! max nb of nodes in a cell 
    call MeshSchema_NbFracCellMax ! max nb of frac in a cell
    call MeshSchema_NbNodeFaceMax ! max nb of nodes in a frac

  end subroutine MeshSchema_make


  ! send/receive (commRank>=1)
  ! or copy (commRank=0) 
  subroutine MeshSchema_sendrecv

    integer :: dest, Ierr, i, j, Nb, Nnnz
    integer stat(MPI_STATUS_SIZE)

    integer :: blen(1), offsets(1), oldtypes(1), MPI_IDNODE
    integer :: blocklen(4), arraytype(4)
    integer(kind=MPI_ADDRESS_KIND) ::disp(4)

    integer :: MPI_DATAWELLINJ
    integer :: blocklen_datawellinj(6), arraytype_datawellinj(6)
    integer(kind=MPI_ADDRESS_KIND) ::disp_datawellinj(6)

    integer :: MPI_DATAWELLPROD
    integer :: blocklen_datawellprod(4), arraytype_datawellprod(4)
    integer(kind=MPI_ADDRESS_KIND) ::disp_datawellprod(4)

    
    ! ************************************* !

    ! Send Nb*

    if (commRank==0) then

       do dest=1,Ncpus-1

          ! Nb*ResS_Ncpus
          call MPI_Send(NbCellResS_Ncpus, Ncpus, MPI_INTEGER, dest, 11, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbFaceResS_Ncpus, Ncpus, MPI_INTEGER, dest, 12, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbNodeResS_Ncpus, Ncpus, MPI_INTEGER, dest, 13, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbFracResS_Ncpus, Ncpus, MPI_INTEGER, dest, 14, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbWellInjResS_Ncpus, Ncpus, MPI_INTEGER, dest, 141, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbWellProdResS_Ncpus, Ncpus, MPI_INTEGER, dest, 142, ComPASS_COMM_WORLD, Ierr)


          ! Nb*OwnS_Ncpus
          call MPI_Send(NbCellOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 15, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbFaceOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 16, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbNodeOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 17, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbFracOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 18, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbWellInjOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 181, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(NbWellProdOwnS_Ncpus, Ncpus, MPI_INTEGER, dest, 182, ComPASS_COMM_WORLD, Ierr)

       end do

       allocate(NbCellLocal_Ncpus(Ncpus))
       allocate(NbFaceLocal_Ncpus(Ncpus))
       allocate(NbNodeLocal_Ncpus(Ncpus))
       allocate(NbFracLocal_Ncpus(Ncpus))
       allocate(NbWellInjLocal_Ncpus(Ncpus))
       allocate(NbWellProdLocal_Ncpus(Ncpus))

       allocate(NbCellOwn_Ncpus(Ncpus))
       allocate(NbFaceOwn_Ncpus(Ncpus))
       allocate(NbNodeOwn_Ncpus(Ncpus))
       allocate(NbFracOwn_Ncpus(Ncpus))
       allocate(NbWellInjOwn_Ncpus(Ncpus))
       allocate(NbWellProdOwn_Ncpus(Ncpus))

       NbCellLocal_Ncpus(:) = NbCellResS_Ncpus(:)
       NbFaceLocal_Ncpus(:) = NbFaceResS_Ncpus(:)
       NbNodeLocal_Ncpus(:) = NbNodeResS_Ncpus(:)
       NbFracLocal_Ncpus(:) = NbFracResS_Ncpus(:)
       NbWellInjLocal_Ncpus(:) = NbWellInjResS_Ncpus(:)
       NbWellProdLocal_Ncpus(:) = NbWellProdResS_Ncpus(:)
       
       NbCellOwn_Ncpus(:) = NbCellOwnS_Ncpus(:)
       NbFaceOwn_Ncpus(:) = NbFaceOwnS_Ncpus(:)
       NbNodeOwn_Ncpus(:) = NbNodeOwnS_Ncpus(:)
       NbFracOwn_Ncpus(:) = NbFracOwnS_Ncpus(:)
       NbWellInjOwn_Ncpus(:) = NbWellInjOwnS_Ncpus(:)
       NbWellProdOwn_Ncpus(:) = NbWellProdOwnS_Ncpus(:)

    else
       allocate(NbCellLocal_Ncpus(Ncpus))
       allocate(NbFaceLocal_Ncpus(Ncpus))
       allocate(NbNodeLocal_Ncpus(Ncpus))
       allocate(NbFracLocal_Ncpus(Ncpus))
       allocate(NbWellInjLocal_Ncpus(Ncpus))
       allocate(NbWellProdLocal_Ncpus(Ncpus))

       allocate(NbCellOwn_Ncpus(Ncpus))
       allocate(NbFaceOwn_Ncpus(Ncpus))
       allocate(NbNodeOwn_Ncpus(Ncpus))
       allocate(NbFracOwn_Ncpus(Ncpus))
       allocate(NbWellInjOwn_Ncpus(Ncpus))
       allocate(NbWellProdOwn_Ncpus(Ncpus))

       ! Nb*ResS_Ncpus
       call MPI_recv(NbCellLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbFaceLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbNodeLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbFracLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbWellInjLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 141, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbWellProdLocal_Ncpus, Ncpus, MPI_INTEGER, 0, 142, ComPASS_COMM_WORLD, stat, Ierr)

       ! Nb*OwnS_Ncpus
       call MPI_Recv(NbCellOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 15, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbFaceOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 16, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbNodeOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 17, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbFracOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 18, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbWellInjOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 181, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_Recv(NbWellProdOwn_Ncpus, Ncpus, MPI_INTEGER, 0, 182, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       deallocate(NbCellResS_Ncpus)
       deallocate(NbFaceResS_Ncpus)
       deallocate(NbNodeResS_Ncpus)
       deallocate(NbFracResS_Ncpus)
       deallocate(NbWellInjResS_Ncpus)
       deallocate(NbWellProdResS_Ncpus)
       deallocate(NbCellOwnS_Ncpus)
       deallocate(NbFaceOwnS_Ncpus)
       deallocate(NbNodeOwnS_Ncpus)
       deallocate(NbFracOwnS_Ncpus)
       deallocate(NbWellInjOwnS_Ncpus)
       deallocate(NbWellProdOwnS_Ncpus)
    end if


    ! Send mesh size
    if (commRank==0) then

       do dest=1,Ncpus-1

          ! meshSize_*max, meshSize_*min
          call MPI_Send(meshSizeS_xmax, 1, MPI_DOUBLE, dest, 11, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(meshSizeS_xmin, 1, MPI_DOUBLE, dest, 12, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(meshSizeS_ymax, 1, MPI_DOUBLE, dest, 13, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(meshSizeS_ymin, 1, MPI_DOUBLE, dest, 14, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(meshSizeS_zmax, 1, MPI_DOUBLE, dest, 15, ComPASS_COMM_WORLD, Ierr)
          call MPI_Send(meshSizeS_zmin, 1, MPI_DOUBLE, dest, 16, ComPASS_COMM_WORLD, Ierr)
       end do

       meshSize_xmax = meshSizeS_xmax
       meshSize_xmin = meshSizeS_xmin
       meshSize_ymax = meshSizeS_ymax
       meshSize_ymin = meshSizeS_ymin
       meshSize_zmax = meshSizeS_zmax
       meshSize_zmin = meshSizeS_zmin
    else

       call MPI_recv(meshSize_xmax, 1, MPI_DOUBLE, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_recv(meshSize_xmin, 1, MPI_DOUBLE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_recv(meshSize_ymax, 1, MPI_DOUBLE, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_recv(meshSize_ymin, 1, MPI_DOUBLE, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_recv(meshSize_zmax, 1, MPI_DOUBLE, 0, 15, ComPASS_COMM_WORLD, stat, Ierr)
       call MPI_recv(meshSize_zmin, 1, MPI_DOUBLE, 0, 16, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    ! **************************************** !

    ! Send XNodeRes_Ncpus to XNodeLocal
    if (commRank==0) then

       do i=1,Ncpus-1
          call MPI_Send(XNodeRes_Ncpus(i+1)%Array2d, NbNodeLocal_Ncpus(i+1)*3, MPI_DOUBLE, i, 21, ComPASS_COMM_WORLD, Ierr)
       end do

       allocate(XNodeLocal(3,NbNodeLocal_Ncpus(1)) )
       do j=1,NbNodeLocal_Ncpus(1)
          XNodeLocal(1,j) = XNodeRes_Ncpus(1)%Array2d(1,j)
          XNodeLocal(2,j) = XNodeRes_Ncpus(1)%Array2d(2,j)
          XNodeLocal(3,j) = XNodeRes_Ncpus(1)%Array2d(3,j)
       end do

    else
       allocate(XNodeLocal(3,NbNodeLocal_Ncpus(commRank+1)) )
       call MPI_Recv(XNodeLocal,  NbNodeLocal_Ncpus(commRank+1)*3, MPI_DOUBLE, 0, 21, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    ! Send node flags    
    if (commRank==0) then

       do i=1,Ncpus-1
          call MPI_Send(NodeFlags_Ncpus(i+1)%Val, NbNodeLocal_Ncpus(i+1), MPI_INTEGER, i, 22, ComPASS_COMM_WORLD, Ierr)
       end do

       allocate(NodeFlagsLocal(NbNodeLocal_Ncpus(1)))
       NodeFlagsLocal = NodeFlags_Ncpus(1)%Val
       
    else
       allocate(NodeFlagsLocal(NbNodeLocal_Ncpus(commRank+1)))
       call MPI_Recv(NodeFlagsLocal,  NbNodeLocal_Ncpus(commRank+1), MPI_INTEGER, 0, 22, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(XNodeRes_Ncpus(i)%Array2d)
          deallocate(NodeFlags_Ncpus(i)%Val)
       end do
       deallocate(XNodeRes_Ncpus)
       deallocate(NodeFlags_Ncpus)
    end if


    ! ************************************** !

    ! Send mesh and connectivities
    if (commRank==0) then

       ! for proc >=1, send
       do i=1,Ncpus-1

          ! NodebyNodeOwn
          call MeshSchema_csrsend(NodebyNodeOwn_Ncpus(i+1), i, 100, VALSIZE_ZERO)
          ! FracbyNodeOwn
          call MeshSchema_csrsend(FracbyNodeOwn_Ncpus(i+1), i, 200, VALSIZE_ZERO)
          ! CellbyNodeOwn
          call MeshSchema_csrsend(CellbyNodeOwn_Ncpus(i+1), i, 300, VALSIZE_ZERO)

          ! NodebyFracOwn
          call MeshSchema_csrsend(NodebyFracOwn_Ncpus(i+1), i, 400, VALSIZE_ZERO)
          ! CellbyFracOwn
          call MeshSchema_csrsend(CellbyFracOwn_Ncpus(i+1), i, 500, VALSIZE_ZERO)
          ! FracbyFracOwn
          call MeshSchema_csrsend(FracbyFracOwn_Ncpus(i+1), i, 600, VALSIZE_ZERO)
          ! WellInjbyNodeOwn
          call MeshSchema_csrsend(WellInjbyNodeOwn_Ncpus(i+1), i, 610, VALSIZE_ZERO)
          ! WellProdbyNodeOwn
          call MeshSchema_csrsend(WellProdbyNodeOwn_Ncpus(i+1), i, 620, VALSIZE_ZERO)

          ! FacebyCellLocal
          call MeshSchema_csrsend(FacebyCellLocal_Ncpus(i+1), i, 700, VALSIZE_ZERO)
          ! FracbyCellLocal
          call MeshSchema_csrsend(FracbyCellLocal_Ncpus(i+1), i, 800, VALSIZE_ZERO)      
          ! NodebyCellLocal
          call MeshSchema_csrsend(NodebyCellLocal_Ncpus(i+1), i, 900, VALSIZE_ZERO)

          ! NodebyFaceLocal
          call MeshSchema_csrsend(NodebyFaceLocal_Ncpus(i+1), i, 1000, VALSIZE_ZERO)

          ! NodebyWellInjLocal
          call MeshSchema_csrsend(NodebyWellInjLocal_Ncpus(i+1), i, 1100, VALSIZE_ZERO)
          ! NodebyWellProdLocal
          call MeshSchema_csrsend(NodebyWellProdLocal_Ncpus(i+1), i, 1200, VALSIZE_ZERO)

       end do

       ! for proc 0, copy
       call CommonType_csrcopy(NodebyNodeOwn_Ncpus(1), NodebyNodeOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(FracbyNodeOwn_Ncpus(1), FracbyNodeOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(CellbyNodeOwn_Ncpus(1), CellbyNodeOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(WellInjbyNodeOwn_Ncpus(1), WellInjbyNodeOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(WellProdbyNodeOwn_Ncpus(1), WellProdbyNodeOwn, VALSIZE_ZERO)

       call CommonType_csrcopy(NodebyFracOwn_Ncpus(1), NodebyFracOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(CellbyFracOwn_Ncpus(1), CellbyFracOwn, VALSIZE_ZERO)
       call CommonType_csrcopy(FracbyFracOwn_Ncpus(1), FracbyFracOwn, VALSIZE_ZERO)

       call CommonType_csrcopy(FacebyCellLocal_Ncpus(1), FacebyCellLocal, VALSIZE_ZERO)
       call CommonType_csrcopy(FracbyCellLocal_Ncpus(1), FracbyCellLocal, VALSIZE_ZERO)
       call CommonType_csrcopy(NodebyCellLocal_Ncpus(1), NodebyCellLocal, VALSIZE_ZERO)
       call CommonType_csrcopy(NodebyFaceLocal_Ncpus(1), NodebyFaceLocal, VALSIZE_ZERO)

       call CommonType_csrcopy(NodebyWellInjLocal_Ncpus(1), NodebyWellInjLocal, VALSIZE_ZERO)
       call CommonType_csrcopy(NodebyWellProdLocal_Ncpus(1), NodebyWellProdLocal, VALSIZE_ZERO)

    else

       ! NodebyNodeOwn
       call MeshSchema_csrrecv(NodebyNodeOwn, 0, 100, VALSIZE_ZERO)
       ! FracbyNodeOwn
       call MeshSchema_csrrecv(FracbyNodeOwn, 0, 200, VALSIZE_ZERO)
       ! CellbyNodeOwn
       call MeshSchema_csrrecv(CellbyNodeOwn, 0, 300, VALSIZE_ZERO)

       ! NodebyFracOwn
       call MeshSchema_csrrecv(NodebyFracOwn, 0, 400, VALSIZE_ZERO)
       ! CellbyFracOwn
       call MeshSchema_csrrecv(CellbyFracOwn, 0, 500, VALSIZE_ZERO)
       ! FracbyFracOwn
       call MeshSchema_csrrecv(FracbyFracOwn, 0, 600, VALSIZE_ZERO)

       ! WellInjbyNodeOwn
       call MeshSchema_csrrecv(WellInjbyNodeOwn, 0, 610, VALSIZE_ZERO)
       ! WellProdbyNodeOwn
       call MeshSchema_csrrecv(WellProdbyNodeOwn, 0, 620, VALSIZE_ZERO)

       ! FacebyCellLocal
       call MeshSchema_csrrecv(FacebyCellLocal, 0, 700, VALSIZE_ZERO)
       ! FracbyCellLocal
       call MeshSchema_csrrecv(FracbyCellLocal, 0, 800, VALSIZE_ZERO)
       ! NodebyCellLocal
       call MeshSchema_csrrecv(NodebyCellLocal, 0, 900, VALSIZE_ZERO)
       ! NodebyFaceLocal
       call MeshSchema_csrrecv(NodebyFaceLocal, 0, 1000, VALSIZE_ZERO)

       ! NodebyWellInjLocal
       call MeshSchema_csrrecv(NodebyWellInjLocal, 0, 1100, VALSIZE_ZERO)
       ! NodebyWellProdLocal
       call MeshSchema_csrrecv(NodebyWellProdLocal, 0, 1200, VALSIZE_ZERO)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          call CommonType_deallocCSR(NodebyNodeOwn_Ncpus(i))
          call CommonType_deallocCSR(FracbyNodeOwn_Ncpus(i))
          call CommonType_deallocCSR(CellbyNodeOwn_Ncpus(i))

          call CommonType_deallocCSR(NodebyFracOwn_Ncpus(i))
          call CommonType_deallocCSR(CellbyFracOwn_Ncpus(i))
          call CommonType_deallocCSR(FracbyFracOwn_Ncpus(i))

          call CommonType_deallocCSR(WellInjbyNodeOwn_Ncpus(i))
          call CommonType_deallocCSR(WellProdbyNodeOwn_Ncpus(i))

          call CommonType_deallocCSR(FacebyCellLocal_Ncpus(i))
          call CommonType_deallocCSR(FracbyCellLocal_Ncpus(i))
          call CommonType_deallocCSR(NodebyCellLocal_Ncpus(i))
          call CommonType_deallocCSR(NodebyFaceLocal_Ncpus(i))

          call CommonType_deallocCSR(NodebyWellInjLocal_Ncpus(i))
          call CommonType_deallocCSR(NodebyWellProdLocal_Ncpus(i))
       end do

       deallocate(NodebyNodeOwn_Ncpus)
       deallocate(FracbyNodeOwn_Ncpus)
       deallocate(CellbyNodeOwn_Ncpus)

       deallocate(NodebyFracOwn_Ncpus)
       deallocate(CellbyFracOwn_Ncpus)
       deallocate(FracbyFracOwn_Ncpus)

       deallocate(WellInjbyNodeOwn_Ncpus)
       deallocate(WellProdbyNodeOwn_Ncpus)

       deallocate(FacebyCellLocal_Ncpus)
       deallocate(FracbyCellLocal_Ncpus)
       deallocate(NodebyCellLocal_Ncpus)
       deallocate(NodebyFaceLocal_Ncpus)

       deallocate(NodebyWellInjLocal_Ncpus)
       deallocate(NodebyWellProdLocal_Ncpus)
    end if

    ! send DataWellInjLocal

    blocklen_datawellinj(1) = 1
    blocklen_datawellinj(2) = 1
    blocklen_datawellinj(3) = NbComp
    blocklen_datawellinj(4) = 1
    blocklen_datawellinj(5) = 1
    blocklen_datawellinj(6) = 1
    arraytype_datawellinj(1) = MPI_DOUBLE
    arraytype_datawellinj(2) = MPI_DOUBLE
    arraytype_datawellinj(3) = MPI_DOUBLE
    arraytype_datawellinj(4) = MPI_DOUBLE
    arraytype_datawellinj(5) = MPI_DOUBLE
    arraytype_datawellinj(6) = MPI_CHARACTER
    disp_datawellinj(1) = 0 !   = 0
    disp_datawellinj(2) = 8 !   + double
    disp_datawellinj(3) = 16 !  + double
    disp_datawellinj(4) = 8*(NbComp+2) ! + double * NbComp
    disp_datawellinj(5) = 8*(NbComp+2) + 8 ! + double 
    disp_datawellinj(6) = 8*(NbComp+2) +16 ! + double

    ! Create and commit
    call MPI_Type_Create_Struct(6, blocklen_datawellinj, disp_datawellinj, arraytype_datawellinj, MPI_DATAWELLINJ, Ierr)
    call MPI_Type_commit(MPI_DATAWELLINJ, Ierr)

    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbWellInjLocal_Ncpus(i+1)
          call MPI_Send(DataWellInjRes_Ncpus(:,i+1), Nb, MPI_DATAWELLINJ, i, 210+i, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbWellInjLocal_Ncpus(1)
       allocate(DataWellInjLocal(Nb))
       DataWellInjLocal(:) = DataWellInjRes_Ncpus(1:Nb,1)
    else

       Nb = NbWellInjLocal_Ncpus(commRank+1)
       allocate(DataWellInjLocal(Nb))
       call MPI_Recv(DataWellInjLocal, Nb, MPI_DATAWELLINJ, 0, 210+commRank, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    ! Free TYPE MPI_DATAWELLINJ
    call MPI_Type_free(MPI_DATAWELLINJ, Ierr)


    ! send DataWellProdLocal
    blocklen_datawellprod(1) = 1
    blocklen_datawellprod(2) = 1
    blocklen_datawellprod(3) = 1
    blocklen_datawellprod(4) = 1
    arraytype_datawellprod(1) = MPI_DOUBLE
    arraytype_datawellprod(2) = MPI_DOUBLE
    arraytype_datawellprod(3) = MPI_DOUBLE
    arraytype_datawellprod(4) = MPI_CHARACTER
    disp_datawellprod(1) = 0 !   = 0
    disp_datawellprod(2) = 8 !   + double
    disp_datawellprod(3) = 16 !  + double
    disp_datawellprod(4) = 24 !  + double 

    ! Create and commit
    call MPI_Type_Create_Struct(4, blocklen_datawellprod, disp_datawellprod, arraytype_datawellprod, MPI_DATAWELLPROD, Ierr)
    call MPI_Type_commit(MPI_DATAWELLPROD, Ierr)

    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbWellProdLocal_Ncpus(i+1)
          call MPI_Send(DataWellProdRes_Ncpus(:,i+1), Nb, MPI_DATAWELLPROD, i, 210+i, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbWellProdLocal_Ncpus(1)
       allocate(DataWellProdLocal(Nb))
       DataWellProdLocal(:) = DataWellProdRes_Ncpus(1:Nb,1)

    else

       Nb = NbWellProdLocal_Ncpus(commRank+1)
       allocate(DataWellProdLocal(Nb))
       call MPI_Recv(DataWellProdLocal, Nb, MPI_DATAWELLPROD, 0, 210+commRank, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    ! Free TYPE MPI_DATAWELLPROD
    call MPI_Type_free(MPI_DATAWELLPROD, Ierr)

    
    ! ************************************* !

    ! Send IdCellLocal
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbCellLocal_Ncpus(i+1)
          call MPI_Send(IdCellRes_Ncpus(i+1)%Val, Nb, MPI_INTEGER, i, 11, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbCellLocal_Ncpus(1)
       allocate(IdCellLocal(Nb))
       IdCellLocal(:) = IdCellRes_Ncpus(1)%Val(:)

    else   
       Nb = NbCellLocal_Ncpus(commRank+1)
       allocate(IdCellLocal(Nb))
       call MPI_Recv(IdCellLocal, Nb, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(IdCellRes_Ncpus(i)%Val)
       end do
       deallocate(IdCellRes_Ncpus)
    end if

    ! Send IdFaceLocal
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbFaceLocal_Ncpus(i+1)
          call MPI_Send(IdFaceRes_Ncpus(i+1)%Val, Nb, MPI_INTEGER, i, 11, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbFaceLocal_Ncpus(1)
       allocate(IdFaceLocal(Nb))
       IdFaceLocal(:) = IdFaceRes_Ncpus(1)%Val(:)

    else   
       Nb = NbFaceLocal_Ncpus(commRank+1)
       allocate(IdFaceLocal(Nb))
       call MPI_Recv(IdFaceLocal, Nb, MPI_INTEGER, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(IdFaceRes_Ncpus(i)%Val)
       end do
       deallocate(IdFaceRes_Ncpus)
    end if

    ! ************************************* !

    ! new MPI type: MPI_IDNODE 
    blen(1) = 4
    offsets(1) = 0
    oldtypes(1) = MPI_CHARACTER

    call MPI_Type_struct(1, blen, offsets, oldtypes, MPI_IDNODE, Ierr)
    call MPI_Type_commit(MPI_IDNODE, Ierr)

    ! Send IdNodeLocal
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbNodeLocal_Ncpus(i+1)
          call MPI_Send(IdNodeRes_Ncpus(i+1)%Val, Nb, MPI_IDNODE, i, 12, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbNodeLocal_Ncpus(1)
       allocate(IdNodeLocal(Nb))
       IdNodeLocal(:) = IdNodeRes_Ncpus(1)%Val(:)

    else   
       Nb = NbNodeLocal_Ncpus(commRank+1)
       allocate(IdNodeLocal(Nb))
       call MPI_Recv(IdNodeLocal, Nb, MPI_IDNODE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    ! free MPI_IDNODE
    call MPI_Type_free(MPI_IDNODE, Ierr)

    ! free IdNode
    if(commRank==0) then
       do i=1, Ncpus
          deallocate(IdNodeRes_Ncpus(i)%Val)
       end do
       deallocate(IdNodeRes_Ncpus)
    end if

    ! ************************************ !

    ! Send FracToFace and FaceToFrac
    if (commRank==0) then

       ! proc>=1, send
       do i=1, Ncpus-1
          Nb = NbFracLocal_Ncpus(i+1)
          call MPI_Send(FracToFaceLocal_Ncpus(i+1)%Val, Nb, MPI_INTEGER, i, 21, ComPASS_COMM_WORLD, Ierr) ! FracToFace

          Nb = NbFaceLocal_Ncpus(i+1)
          call MPI_Send(FaceToFracLocal_Ncpus(i+1)%Val, Nb, MPI_INTEGER, i, 22, ComPASS_COMM_WORLD, Ierr) ! FracToFace
       end do

       ! proc=0, copy
       Nb = NbFracLocal_Ncpus(1)
       allocate(FracToFaceLocal(Nb))
       FracToFaceLocal(:) = FracToFaceLocal_Ncpus(1)%Val(:) ! FracToFace

       Nb = NbFaceLocal_Ncpus(1)
       allocate(FaceToFracLocal(Nb))
       FaceToFracLocal(:) = FaceToFracLocal_Ncpus(1)%Val(:) ! FaceToFrac

    else
       Nb = NbFracLocal_Ncpus(commRank+1)
       allocate(FracToFaceLocal(Nb))
       call MPI_Recv(FracToFaceLocal, Nb, MPI_INTEGER, 0, 21, ComPASS_COMM_WORLD, stat, Ierr) ! FracToFace

       Nb = NbFaceLocal_Ncpus(commRank+1)
       allocate(FaceToFracLocal(Nb))
       call MPI_Recv(FaceToFracLocal, Nb, MPI_INTEGER, 0, 22, ComPASS_COMM_WORLD, stat, Ierr) ! FaceToFrac
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(FracToFaceLocal_Ncpus(i)%Val)
          deallocate(FaceToFracLocal_Ncpus(i)%Val)
       end do
       deallocate(FracToFaceLocal_Ncpus)
       deallocate(FaceToFracLocal_Ncpus)
    end if


    ! ********************************** !

    ! Send NumNodebyProc, NumFracbyProc, NumWellInjbyProc, NumWellProdbyProc
    if(commRank==0) then
       do i=1,Ncpus-1
          call MeshSchema_csrsend(NumNodebyProc_Ncpus(i+1), i, 100, VALSIZE_NNZ)
          call MeshSchema_csrsend(NumFracbyProc_Ncpus(i+1), i, 300, VALSIZE_NNZ)
          call MeshSchema_csrsend(NumWellInjbyProc_Ncpus(i+1), i, 500, VALSIZE_NNZ)
          call MeshSchema_csrsend(NumWellProdbyProc_Ncpus(i+1), i, 700, VALSIZE_NNZ)
       end do

       call CommonType_csrcopy(NumNodebyProc_Ncpus(1), NumNodebyProc, VALSIZE_NNZ)
       call CommonType_csrcopy(NumFracbyProc_Ncpus(1), NumFracbyProc, VALSIZE_NNZ)
       call CommonType_csrcopy(NumWellInjbyProc_Ncpus(1), NumWellInjbyProc, VALSIZE_NNZ)
       call CommonType_csrcopy(NumWellProdbyProc_Ncpus(1), NumWellProdbyProc, VALSIZE_NNZ)

    else
       call MeshSchema_csrrecv(NumNodebyProc, 0, 100, VALSIZE_NNZ)
       call MeshSchema_csrrecv(NumFracbyProc, 0, 300, VALSIZE_NNZ)
       call MeshSchema_csrrecv(NumWellInjbyProc, 0, 500, VALSIZE_NNZ)
       call MeshSchema_csrrecv(NumWellProdbyProc, 0, 700, VALSIZE_NNZ)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          call CommonType_deallocCSR(NumNodebyProc_Ncpus(i))
          call CommonType_deallocCSR(NumFracbyProc_Ncpus(i))
          call CommonType_deallocCSR(NumWellInjbyProc_Ncpus(i))
          call CommonType_deallocCSR(NumWellProdbyProc_Ncpus(i))
       end do
       deallocate(NumNodebyProc_Ncpus)
       deallocate(NumFracbyProc_Ncpus)
       deallocate(NumWellInjbyProc_Ncpus)
       deallocate(NumWellProdbyProc_Ncpus)
    end if

    ! ******************************* !

    ! send PorositeCell
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbCellLocal_Ncpus(i+1)
          call MPI_Send(PorositeCell_Ncpus(i+1)%Val, Nb, MPI_DOUBLE, i, 11, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbCellLocal_Ncpus(1)
       allocate(PorositeCellLocal(Nb))
       PorositeCellLocal(:) = PorositeCell_Ncpus(1)%Val(:)

    else   
       Nb = NbCellLocal_Ncpus(commRank+1)
       allocate(PorositeCellLocal(Nb))
       call MPI_Recv(PorositeCellLocal, Nb, MPI_DOUBLE, 0, 11, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(PorositeCell_Ncpus(i)%Val)
       end do
       deallocate(PorositeCell_Ncpus)
    end if

    ! send PorositeFrac
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbFracLocal_Ncpus(i+1)
          call MPI_Send(PorositeFrac_Ncpus(i+1)%Val, Nb, MPI_DOUBLE, i, 12, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbFracLocal_Ncpus(1)
       allocate(PorositeFracLocal(Nb))
       PorositeFracLocal(:) = PorositeFrac_Ncpus(1)%Val(:)

    else   
       Nb = NbFracLocal_Ncpus(commRank+1)
       allocate(PorositeFracLocal(Nb))
       call MPI_Recv(PorositeFracLocal, Nb, MPI_DOUBLE, 0, 12, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(PorositeFrac_Ncpus(i)%Val)
       end do
       deallocate(PorositeFrac_Ncpus)
    end if

    ! send PermCellLocal
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbCellLocal_Ncpus(i+1)
          call MPI_Send(PermCellLocal_Ncpus(i+1)%Array3d, Nb*9, MPI_DOUBLE, i, 13, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbCellLocal_Ncpus(1)
       allocate(PermCellLocal(3,3,Nb))
       PermCellLocal(:,:,:) = PermCellLocal_Ncpus(1)%Array3d(:,:,:)

    else   
       Nb = NbCellLocal_Ncpus(commRank+1)
       allocate(PermCellLocal(3,3,Nb))
       call MPI_Recv(PermCellLocal, Nb*9, MPI_DOUBLE, 0, 13, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          !          deallocate(PermCellLocal_Ncpus(i+1)%Array3d)    ! ???
       end do
       ! deallocate(PermCellLocal_Ncpus)   ! ???
    end if

    ! send PermFrac
    if (commRank==0) then

       ! proc >=1, send
       do i=1, Ncpus-1
          Nb = NbFracLocal_Ncpus(i+1)
          call MPI_Send(PermFracLocal_Ncpus(i+1)%Val, Nb, MPI_DOUBLE, i, 14, ComPASS_COMM_WORLD, Ierr)
       end do

       ! proc=0, copy
       Nb = NbFracLocal_Ncpus(1)
       allocate(PermFracLocal(Nb))
       PermFracLocal(:) = PermFracLocal_Ncpus(1)%Val(:)

    else   
       Nb = NbFracLocal_Ncpus(commRank+1)
       allocate(PermFracLocal(Nb))
       call MPI_Recv(PermFracLocal, Nb, MPI_DOUBLE, 0, 14, ComPASS_COMM_WORLD, stat, Ierr)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          deallocate(PermFracLocal_Ncpus(i)%Val)
       end do
       !       deallocate(PermFracLocal_Ncpus)   ! ???
    end if



    ! MPI TYPE for DataNodewell: MPI_DATANODEWELL
    blocklen(1) = 1
    blocklen(2) = 1
    blocklen(3) = 1
    blocklen(4) = 1
    arraytype(1) = MPI_INTEGER
    arraytype(2) = MPI_INTEGER
    arraytype(3) = MPI_DOUBLE
    arraytype(4) = MPI_DOUBLE
    disp(1) = 0 !  = 0
    disp(2) = 4 !  + integer
    disp(3) = 8 !  + integer
    disp(4) = 16 ! + double

    ! Create and commit
    call MPI_Type_Create_Struct(4, blocklen, disp, arraytype, MPI_DATANODEWELL, Ierr)
    call MPI_Type_commit(MPI_DATANODEWELL, Ierr)

    ! send NodeDatabyWell
    if(commRank==0) then

       do i=1, Ncpus-1
          call MeshSchema_csrdatawellsend(NodeDatabyWellInjLocal_Ncpus(i+1), i, 100)
          call MeshSchema_csrdatawellsend(NodeDatabyWellProdLocal_Ncpus(i+1), i, 200)
       end do

       call DefWell_csrdatawellcopy(NodeDatabyWellInjLocal_Ncpus(1), NodeDatabyWellInjLocal)
       call DefWell_csrdatawellcopy(NodeDatabyWellProdLocal_Ncpus(1), NodeDatabyWellProdLocal)
    else

       call MeshSchema_csrdatawellrecv(NodeDatabyWellInjLocal, 0, 100)
       call MeshSchema_csrdatawellrecv(NodeDatabyWellProdLocal, 0, 200)
    end if

    if(commRank==0) then
       do i=1, Ncpus
          call DefWell_deallocCSRDataWell(NodeDatabyWellInjLocal_Ncpus(i))
          call DefWell_deallocCSRDataWell(NodeDatabyWellProdLocal_Ncpus(i))
       end do
       deallocate(NodeDatabyWellInjLocal_Ncpus)
       deallocate(NodeDatabyWellProdLocal_Ncpus)
    end if

    ! Free TYPE MPI_DATANODEWELL
    call MPI_Type_free(MPI_DATANODEWELL, Ierr)

	!do i=1, size(DataWellInjLocal)
	!	write(*,*) "%%", "local injection data on proc", commRank
	!	call DefWell_print_DataWellInj(DataWellInjLocal(i))
	!end do
	!do i=1, size(DataWellProdLocal)
	!	write(*,*) "%%", "local production data on proc", commRank
	!	call DefWell_print_DataWellProd(DataWellProdLocal(i))
	!end do
    
  end subroutine MeshSchema_sendrecv


  ! Send csr1 to "dest" with "tag"
  subroutine MeshSchema_csrsend(csr1, dest, tag, valsize)

    type(CSR) :: csr1
    integer, intent(in) :: tag, dest, valsize
    integer :: Nb, Nnz, Ierr
    integer stat(MPI_STATUS_SIZE)

    call MPI_Send(csr1%Nb, 1, MPI_INTEGER, dest, tag+1, ComPASS_COMM_WORLD, Ierr)
    Nb = csr1%Nb
    call MPI_Send(csr1%Pt, Nb+1, MPI_INTEGER, dest, tag+2, ComPASS_COMM_WORLD, Ierr)
    Nnz = csr1%Pt(Nb+1) 
    call MPI_Send(csr1%Num, Nnz, MPI_INTEGER, dest, tag+3, ComPASS_COMM_WORLD, Ierr)
    if(valsize==VALSIZE_NB) then
       call MPI_Send(csr1%Val, Nb, MPI_INTEGER, dest, tag+4, ComPASS_COMM_WORLD, Ierr)
    else if(valsize==VALSIZE_NNZ) then
       call MPI_Send(csr1%Val, Nnz, MPI_INTEGER, dest, tag+4, ComPASS_COMM_WORLD, Ierr)
    end if

  end subroutine MeshSchema_csrsend


  ! Recv csr2 from "source" with "tag"
  subroutine MeshSchema_csrrecv(csr2, source, tag, valsize)

    type(CSR), intent(inout) :: csr2
    integer, intent(in) :: tag, source, valsize

    integer :: Nb, Nnz, Ierr
    integer stat(MPI_STATUS_SIZE)

    call MPI_Recv(csr2%Nb, 1, MPI_INTEGER, source, tag+1, ComPASS_COMM_WORLD, stat, Ierr)

    Nb = csr2%Nb
    allocate(csr2%Pt(Nb+1))
    call MPI_Recv(csr2%Pt, Nb+1, MPI_INTEGER, source, tag+2, ComPASS_COMM_WORLD, stat, Ierr)

    Nnz = csr2%Pt(Nb+1)
    allocate(csr2%Num(Nnz))
    call MPI_Recv(csr2%Num, Nnz, MPI_INTEGER, source, tag+3, ComPASS_COMM_WORLD, stat, Ierr)

    if(valsize==VALSIZE_NB) then
       allocate(csr2%Val(Nb))
       call MPI_Recv(csr2%Val, Nb, MPI_INTEGER, source, tag+4, ComPASS_COMM_WORLD, stat, Ierr)
    else if(valsize==VALSIZE_NNZ) then
       allocate(csr2%Val(Nnz))
       call MPI_Recv(csr2%Val, Nnz, MPI_INTEGER, source, tag+4, ComPASS_COMM_WORLD, stat, Ierr)
    end if

  end subroutine MeshSchema_csrrecv


  ! Send csr1 to "dest" with "tag"
  ! %val type is DataNodeWell
  subroutine MeshSchema_csrdatawellsend(csr1, dest, tag)

    type(TYPE_CSRDataNodeWell) :: csr1
    integer, intent(in) :: tag, dest
    integer :: Nb, Nnz, Ierr
    integer stat(MPI_STATUS_SIZE)

    call MPI_Send(csr1%Nb, 1, MPI_INTEGER, dest, tag+1, ComPASS_COMM_WORLD, Ierr)
    Nb = csr1%Nb
    call MPI_Send(csr1%Pt, Nb+1, MPI_INTEGER, dest, tag+2, ComPASS_COMM_WORLD, Ierr)
    Nnz = csr1%Pt(Nb+1) 
    call MPI_Send(csr1%Num, Nnz, MPI_INTEGER, dest, tag+3, ComPASS_COMM_WORLD, Ierr)

    call MPI_Send(csr1%Val, Nnz, MPI_DATANODEWELL, dest, tag+4, ComPASS_COMM_WORLD, Ierr)

  end subroutine MeshSchema_csrdatawellsend


  ! Recv csr2 from "source" with "tag"
  ! %val type is DataNodeWell
  subroutine MeshSchema_csrdatawellrecv(csr2, source, tag)

    type(TYPE_CSRDataNodeWell), intent(inout) :: csr2
    integer, intent(in) :: tag, source

    integer :: Nb, Nnz, Ierr
    integer stat(MPI_STATUS_SIZE)

    call MPI_Recv(csr2%Nb, 1, MPI_INTEGER, source, tag+1, ComPASS_COMM_WORLD, stat, Ierr)

    Nb = csr2%Nb
    allocate(csr2%Pt(Nb+1))
    call MPI_Recv(csr2%Pt, Nb+1, MPI_INTEGER, source, tag+2, ComPASS_COMM_WORLD, stat, Ierr)

    Nnz = csr2%Pt(Nb+1)
    allocate(csr2%Num(Nnz))
    call MPI_Recv(csr2%Num, Nnz, MPI_INTEGER, source, tag+3, ComPASS_COMM_WORLD, stat, Ierr)

    allocate(csr2%Val(Nnz))
    call MPI_Recv(csr2%Val, Nnz, MPI_DATANODEWELL, source, tag+4, ComPASS_COMM_WORLD, stat, Ierr)

  end subroutine MeshSchema_csrdatawellrecv

  ! Output:
  !  NumNodebyEdgebyWellLocal
  ! Input:
  !  NodeDatabyWellLocal
  ! Local list of Edge by Well (own+ghost)
  ! List is oriented (Id_parent, Id_son)
  subroutine MeshSchema_NumNodebyEdgebyWellLocal

    integer :: NbWellLocal, nnz, i, j, NbEdgemax, comptNode

    !! INJ WELL
    NbWellLocal = NodeDatabyWellInjLocal%Nb ! Number of well inj
    allocate(NbEdgebyWellInjLocal(NbWellLocal))
    NbEdgemax = 0

    do i=1,NbWellLocal
       ! Number of edges = (Number of nodes in local well i) - 1
       NbEdgebyWellInjLocal(i) = NodeDatabyWellInjLocal%Pt(i+1) - NodeDatabyWellInjLocal%Pt(i) - 1
       NbEdgemax = max(NbEdgemax,NbEdgebyWellInjLocal(i))
    enddo

    allocate(NumNodebyEdgebyWellInjLocal(2,NbEdgemax,NbWellLocal))
    NumNodebyEdgebyWellInjLocal(:,:,:) = -1

    do i=1,NbWellLocal

       ! comptNode = 0
       ! ! loop over every nodes of local well, minus the head node of well
       ! do j=NodeDatabyWellInjLocal%Pt(i)+1,NodeDatabyWellInjLocal%Pt(i+1)-1
       !    comptNode = comptNode + 1
       !    NumNodebyEdgebyWellInjLocal(1,comptNode,i) = NodeDatabyWellInjLocal%Val(j)%Parent
       !    NumNodebyEdgebyWellInjLocal(2,comptNode,i) = NodeDatabyWellInjLocal%Num(j)
       ! enddo

       comptNode = 0
       
       ! loop over every nodes of local well, minus the head node of well
       do j=1,NodeDatabyWellInjLocal%Pt(i+1)-NodeDatabyWellInjLocal%Pt(i)-1
          NumNodebyEdgebyWellInjLocal(1,comptNode+j,i) = NodeDatabyWellInjLocal%Val(j+NodeDatabyWellInjLocal%Pt(i))%Parent
          NumNodebyEdgebyWellInjLocal(2,comptNode+j,i) = NodeDatabyWellInjLocal%Num(j+NodeDatabyWellInjLocal%Pt(i))
       enddo
       comptNode = comptNode + NodeDatabyWellInjLocal%Pt(i+1)-NodeDatabyWellInjLocal%Pt(i)-1
    enddo

    !! PROD WELL
    NbWellLocal = NodeDatabyWellProdLocal%Nb ! Number of well Prod
    allocate(NbEdgebyWellProdLocal(NbWellLocal))
    NbEdgemax = 0

    do i=1,NbWellLocal
       ! Number of edges = (Number of nodes in local well i) - 1
       NbEdgebyWellProdLocal(i) = NodeDatabyWellProdLocal%Pt(i+1) - NodeDatabyWellProdLocal%Pt(i) - 1
       NbEdgemax = max(NbEdgemax,NbEdgebyWellProdLocal(i))
    enddo

    allocate(NumNodebyEdgebyWellProdLocal(2,NbEdgemax,NbWellLocal))
    NumNodebyEdgebyWellProdLocal(:,:,:) = -1

    do i=1,NbWellLocal
       comptNode = 0
       ! loop over every nodes of local well, minus the head node of well
       do j=NodeDatabyWellProdLocal%Pt(i)+1,NodeDatabyWellProdLocal%Pt(i+1)-1
          comptNode = comptNode + 1
          NumNodebyEdgebyWellProdLocal(1,comptNode,i) = NodeDatabyWellProdLocal%Val(j)%Parent
          NumNodebyEdgebyWellProdLocal(2,comptNode,i) = NodeDatabyWellProdLocal%Num(j)
       enddo
    enddo

  end subroutine MeshSchema_NumNodebyEdgebyWellLocal

  ! center of Cell
  subroutine MeshSchema_XCellLocal

    integer :: k, m, nbCellLocal
    double precision, dimension(3) :: xk(3)

    nbCellLocal = NbCellLocal_Ncpus(commRank+1)

    allocate(XCellLocal(3,nbCellLocal))

    ! boucle sur les mailles
    do k=1,nbCellLocal

       ! center of cell
       xk(:) = 0.d0
       do m = NodebyCellLocal%Pt(k)+1, NodebyCellLocal%Pt(k+1)
          xk(:) = xk(:) + XNodeLocal(:, NodebyCellLocal%Num(m))
       enddo
       xk(:) = xk(:)/dble(NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k))

       XCellLocal(:,k) = xk(:)       
    enddo

  end subroutine MeshSchema_XCellLocal


  ! Vol of Cell
  subroutine MeshSchema_VolCellLocal

    ! calcul du volume et du centre de gravite
    ! maille polyèdrique quelconque non necessairement convexe 
    ! faces non planes (decoupée en triangles avec un point au centre)  

    integer :: i,j,k,m,n1,n2
    double precision :: volk,volT
    double precision, dimension(3) :: yk,xk,xT,x1,x2,xs,e0,e1,e2,e3

    integer :: Ierr, errcode ! used for MPI_Abort

    ! check if XCellLocal is computed
    if (allocated(XCellLocal) .eqv. .false.) then
       if(commRank==0) then
          print*, "XCellLocal: center of cell not computed"
       end if

       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    allocate(VolCellLocal(NbCellLocal_Ncpus(commRank+1)))
    VolCellLocal(:) = 0.d0

    ! boucle sur les mailles
    do k=1, NbCellLocal_Ncpus(commRank+1)

       volk = 0.d0

       ! center of cell
       yk(:) = XCellLocal(:,k)

       ! boucle sur les faces i de la maille k 
       do j = FacebyCellLocal%Pt(k)+1, FacebyCellLocal%Pt(k+1)
          i = FacebyCellLocal%Num(j)

          ! isobarycentre de la face 
          xs(:) = XFaceLocal(:,i)

          ! boucle sur les nodes n1 de la face i 
          do m = NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1)
             n1 = NodebyFaceLocal%Num(m)
             x1(:) = XNodeLocal(:,n1)

             if (m == NodebyFaceLocal%Pt(i+1)) then 
                n2 = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(i)+1)
             else 
                n2 = NodebyFaceLocal%Num(m+1)
             endif

             x2(:) = XNodeLocal(:,n2)
             xT(:) = (x1(:)+x2(:)+xs(:)+yk(:))/4.d0 

             e0(:) = xT(:)-yk(:)
             e1(:) = x1(:)-yk(:)
             e2(:) = x2(:)-yk(:)
             e3(:) = xs(:)-yk(:)

             volT = e1(1)*e2(2)*e3(3) + e2(1)*e3(2)*e1(3) + e3(1)*e1(2)*e2(3) &
                  -e1(1)*e2(3)*e3(2) - e2(1)*e3(3)*e1(2) - e3(1)*e1(3)*e2(2) 

             volT =  abs(volT)/6.d0 
             volk = volk + volT 
          enddo
       enddo

       VolCellLocal(k) = volk
    enddo

  end subroutine MeshSchema_VolCellLocal


  ! center of frac
  subroutine MeshSchema_XFaceLocal

    integer :: ifrac, i, m
    integer :: nbFaceLocal
    double precision :: xf(3)

    nbFaceLocal = NbFaceLocal_Ncpus(commRank+1)

    allocate(XFaceLocal(3,nbFaceLocal))    ! center of face

    ! boucle sur les face frac     
    do i = 1, nbFaceLocal

       ! isobarycentre de la face 
       xf(:) = 0.d0
       do m = NodebyFaceLocal%Pt(i)+1,NodebyFaceLocal%Pt(i+1)
          xf(:) = xf(:) + XNodeLocal(:,NodebyFaceLocal%Num(m))
       enddo
       XFaceLocal(:,i) = xf(:)/dble(NodebyFaceLocal%Pt(i+1)-NodebyFaceLocal%Pt(i))

    end do ! end of loop frac

  end subroutine MeshSchema_XFaceLocal


  ! surf of frac face
  subroutine MeshSchema_SurfFracLocal

    double precision, dimension(3) :: &
         xk, xf, x1, x2, xt, & ! cordinate
         v                     ! normal directive

    double precision :: &
         SurfFace, & ! surface of a face
         Surf12f     ! surface of a triangle with nodes 1,2 and center of face

    integer :: &
         n1, n2, & ! num (local) of a node
         in1, in2  ! num (face) of a node

    integer :: &
         i, ifrac,    & ! i: loop of face frac, ifrac: num (local) of i
         k, j, m

    integer :: errcode, Ierr

    integer :: nbNodeFace

    ! check if XFaceLocal is computed
    if (allocated(XFaceLocal) .eqv. .false.) then
       if(commRank==0) then
          print*, "XFaceLocal: center of cell not computed"
       end if

       call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
    end if

    allocate(SurfFracLocal(NbFracLocal_Ncpus(commRank+1))) ! surf of face
    SurfFracLocal(:) = 0.d0

    ! boucle sur les face frac     
    do ifrac = 1, NbFracLocal_Ncpus(commRank+1)
       i = FracToFaceLocal(ifrac)

       ! num of nodes in face i
       nbNodeFace = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)

       ! isobarycentre de la face 
       xf(:) = XFaceLocal(:,i)

       ! init SurfFace as zero, surface of face i
       SurfFace = 0.d0

       do m = NodebyFaceLocal%Pt(i)+1, NodebyFaceLocal%Pt(i+1)

          ! edge of nodes n1, n2
          ! num (face) of n1 and n2 are in1 and in2            
          in1 = m - NodebyFaceLocal%Pt(i)
          n1 = NodebyFaceLocal%Num(m)

          if (m==NodebyFaceLocal%Pt(i+1)) then 
             n2 = NodebyFaceLocal%Num(NodebyFaceLocal%Pt(i)+1)
             in2 = 1
          else 
             n2 = NodebyFaceLocal%Num(m+1)
             in2 = in1 + 1
          endif

          x1(:) = XNodeLocal(:,n1)
          x2(:) = XNodeLocal(:,n2)

          ! Surface of triangle (arete n1, n2 and center of face) and normal vector
          xt(:) = (x1(:)+x2(:)+xf(:))/3.d0

          call MeshSchema_Surf12f(x1,x2,xf,xt,Surf12f)

          ! Surface of Face i is sum of surf12f
          SurfFace = SurfFace + Surf12f 

       end do ! end of loop edge in face       

       SurfFracLocal(ifrac) = SurfFace ! area of frac

    end do ! end of loop face

  end subroutine MeshSchema_SurfFracLocal


  subroutine MeshSchema_Surf12f(x1,x2,x3,x,surf)

    ! calcul du vecteur normal unitaire d'un triangle defini par les 
    ! coordonnees de ses trois points sortant par rapport 
    ! au point x,y,z > vx,vy,vz 
    ! + surface du triangle = surf 

    double precision, dimension(3), intent(in) :: x1 ,x2, x3, x
    double precision, intent(out) :: surf
    double precision, dimension(3) :: xt, v
    double precision :: s

    xt(:) = (x1(:)+x2(:)+x3(:))/3.d0

    v(1) = (x1(3)-x3(3))*(x1(2)-x2(2)) - (x1(2)-x3(2))*(x1(3)-x2(3))
    v(2) = (x1(1)-x3(1))*(x1(3)-x2(3)) - (x1(3)-x3(3))*(x1(1)-x2(1))
    v(3) = (x1(2)-x3(2))*(x1(1)-x2(1)) - (x1(1)-x3(1))*(x1(2)-x2(2))

    s = dsqrt(v(1)**2+v(2)**2+v(3)**2)
    surf = s/2.d0

  end subroutine MeshSchema_Surf12f



  ! max number of nodes in a cell
  subroutine MeshSchema_NbNodeCellMax

    integer :: k, Nb

    NbNodeCellMax = 0
    do k=1, NbCellLocal_Ncpus(commRank+1)

       Nb = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
       if(NbNodeCellMax<Nb) then
          NbNodeCellMax = Nb
       end if
    end do

  end subroutine MeshSchema_NbNodeCellMax

  ! max number of frac in a cell
  subroutine MeshSchema_NbFracCellMax

    integer :: k, Nb

    NbFracCellMax = 0
    do k = 1, NbCellLocal_Ncpus(commRank+1)

       Nb = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)
       if(NbFracCellMax<Nb) then
          NbFracCellMax = Nb
       end if
    end do

  end subroutine MeshSchema_NbFracCellMax


  ! max number of nodes in a cell
  subroutine MeshSchema_NbNodeFaceMax

    integer :: i, Nb

    NbNodeFaceMax = 0
    do i=1, NbFaceLocal_Ncpus(commRank+1)

       Nb = NodebyFaceLocal%Pt(i+1) - NodebyFaceLocal%Pt(i)
       if(NbNodeFaceMax<Nb) then
          NbNodeFaceMax = Nb
       end if
    end do

  end subroutine MeshSchema_NbNodeFaceMax

  ! Free mesh and some connecivities
  subroutine MeshSchema_Free

    deallocate(NbCellLocal_Ncpus)
    deallocate(NbCellOwn_Ncpus)
    deallocate(NbFaceLocal_Ncpus)
    deallocate(NbFaceOwn_Ncpus)
    deallocate(NbNodeLocal_Ncpus)
    deallocate(NbNodeOwn_Ncpus)
    deallocate(NbFracLocal_Ncpus)
    deallocate(NbFracOwn_Ncpus)

    call CommonType_deallocCSR(FacebyCellLocal)
    call CommonType_deallocCSR(FracbyCellLocal)
    call CommonType_deallocCSR(NodebyCellLocal)
    call CommonType_deallocCSR(NodebyFaceLocal)

    deallocate(XNodeLocal)
    deallocate(NodeFlagsLocal)

    deallocate(IdCellLocal)
    deallocate(IdFaceLocal)
    deallocate(IdNodeLocal)

    deallocate(FracToFaceLocal)
    deallocate(FaceToFracLocal)

    call CommonType_deallocCSR(NumNodebyProc)
    call CommonType_deallocCSR(NumFracbyProc)

    deallocate(NbEdgebyWellInjLocal)
    deallocate(NbEdgebyWellProdLocal)
    deallocate(NumNodebyEdgebyWellInjLocal)
    deallocate(NumNodebyEdgebyWellProdLocal)

    deallocate(XCellLocal)
    deallocate(XFaceLocal)

    deallocate(VolCellLocal)
    deallocate(SurfFracLocal)

    ! the follwoing free could be put after VAGFrac?
    deallocate(PorositeCellLocal)
    deallocate(PorositeFracLocal)

    ! the folllowing free could be put after Assembly
    call CommonType_deallocCSR(NodebyNodeOwn)
    call CommonType_deallocCSR(FracbyNodeOwn)
    call CommonType_deallocCSR(CellbyNodeOwn)

    call CommonType_deallocCSR(NodebyFracOwn)
    call CommonType_deallocCSR(CellbyFracOwn)
    call CommonType_deallocCSR(FracbyFracOwn)

  end subroutine MeshSchema_Free

end module MeshSchema
