!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module NN

   use PathUtilities

   use GlobalMesh
   use PartitionMesh
   use LocalMesh
   use MeshSchema
   use SchemeParameters

   use DefModel
   use NumbyContext
   use IncCV
   use IncCVReservoir
   use IncCVWells
   use DirichletContribution
   use NeumannContribution
   use VAGFrac

   use LoisThermoHydro
   use Flux
   use Residu
   use Jacobian
   use SolvePetsc

   use DefFlash
   use DefFlashWells

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h> 
#else
#include <petsc/finclude/petsc.h>
#endif

   use petsc

!! FIXME: This is define through CMake options #define _VISU_
#ifdef _VISU_
   use VisuVTK
#endif

   implicit none

   ! Mesh file and Perm file
   character(len=200) :: Wellinfoname

   character(len=200) :: output_path

   integer :: Ierr, errcode
   logical :: file_exists

   integer :: i, j, k, s, iph, fi, numj, head
   integer :: nsf, is, n, fk
   integer, allocatable, dimension(:) :: fd
   integer :: rowk, nz, colk
   double precision, pointer :: ptr(:)

   double precision :: err_cell_L1, err_cell_L2, err_cell_Linf
   double precision :: errlocal_cell_L1, errlocal_cell_L2, errlocal_cell_Linf
   double precision :: err_frac_L1, err_frac_L2, err_frac_Linf
   double precision :: errlocal_frac_L1, errlocal_frac_L2, errlocal_frac_Linf

   double precision :: sol

   double precision :: Tempmaxloc, Tempminloc, Tempmax, Tempmin

   ! Time variables
   double precision :: Delta_t, TimeCurrent, TimeOutput
   integer :: VisuTimeIter = 1
   double precision :: Psat, dTsat

   ! Computation time variables
   double precision :: &
      comptime_total, &
      comptime_timestep, &
      comptime_start, &
      comptime_part

   double precision :: &
      comptime_readmesh !, &
   !comptime_meshmake

   ! Newton variables
   logical :: NewtonConv
   integer :: NewtonIter
   double precision :: NewtonResNormRel(NewtonNiterMax)
   double precision :: NewtonRelax
   integer :: NewtonNiterTotal = 0
   integer :: NewtonNbFailure = 0
   double precision :: & ! Residu relative norm from first Newton iteration
      NewtonResConvInit(NbCompThermique), &
      NewtonResClosInit

   ! ksp variables
   double precision, allocatable, dimension(:) :: KspHistory

   logical :: KspConv
   integer :: KspNiter, KspNiterTimeStep
   integer :: KspNiterTotal = 0
   integer :: KspNbFailure = 0

   ! Newton increment (NbInc,NbNodelocal)
   double precision, dimension(:, :), allocatable :: &
      NewtonIncreNode, &
      NewtonIncreFrac, &
      NewtonIncreCell

   ! as well have only one unknown (pressure), only a vector (NbWellLocal)
   double precision, dimension(:), allocatable :: &
      NewtonIncreWellInj, &
      NewtonIncreWellProd

   !  ! Perm
   !  double precision, dimension(:,:,:), allocatable :: PermCellLocal
   !  double precision, dimension(:), allocatable :: PermFracLocal

#ifdef _VISU_
   ! vectors used for visu
   double precision, dimension(:), allocatable :: datavisucell
   double precision, dimension(:), allocatable :: datavisufrac
   double precision, dimension(:), allocatable :: datavisunode
   double precision, dimension(:), allocatable :: datavisuwellinj
   double precision, dimension(:), allocatable :: datavisuwellprod
   
   
#endif

   double precision :: visutime

   double precision :: f, CC(NbComp), SS(NbPhase), dPf, dTf, dCf(NbComp), dSf(NbPhase)

   ! ! ********************************** ! !

   public :: &
      NN_init, &
      NN_main, &
      NN_main_make_timestep, &
      NN_main_output_visu, &
      NN_main_summarize_timestep, &
      NN_finalize

contains

   subroutine NN_init_output_streams(Logfile)

      character(len=*), intent(in) :: LogFile

      ! Report file
      if (commRank == 0) then
         open (11, file=trim(LogFile), status="unknown")
      end if

      allocate (fd(2)); fd = (/6, 11/) ! stdout: 6, logfile: 11
      ! allocate(fd(1)); fd = (/11/)

   end subroutine NN_init_output_streams

   subroutine NN_init_read_mesh(MeshFile)

      character(len=*), intent(in) :: MeshFile

      !FIXME: This is more of an assertion for consistency, might be removed
      if (.NOT. commRank == 0) then
         print *, "Mesh is supposed to be read by master process."
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      inquire (FILE=MeshFile, EXIST=file_exists)
      if (file_exists .eqv. .false.) then
         print *, " "
         print *, "Mesh does not exist   ", MeshFile
         print *, " "
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if
      print *, "Mesh read from file: ", MeshFile

      ! Read Global Mesh
      call GlobalMesh_Make_read_file(MeshFile)

      call GlobalMesh_Make_post_read()

      call DefWell_Make_SetDataWell(NbWellInj, NbWellProd)
      call DefWell_Make_ComputeWellIndex( &
         NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
         PermCell, PermFrac)

   end subroutine NN_init_read_mesh

   subroutine NN_init_build_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)

      real(kind=c_double), intent(in)  :: Ox, Oy, Oz
      real(kind=c_double), intent(in)  :: lx, ly, lz
      integer(kind=c_int), intent(in)  :: nx, ny, nz

      !FIXME: This is more of an assertion for consistency, might be removed
      if (.NOT. commRank == 0) then
         print *, "Mesh is supposed to be built by master process."
         !CHECKME: MPI_Abort is supposed to end all MPI processes
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      call GlobalMesh_Build_cartesian_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)

      call GlobalMesh_SetWellCar(nx, ny, nz)

      call GlobalMesh_Make_post_read()

      call DefWell_Make_SetDataWell(NbWellInj, NbWellProd)
      call DefWell_Make_ComputeWellIndex( &
         NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
         PermCell, PermFrac)

   end subroutine NN_init_build_grid

   subroutine NN_partition_mesh(status)

      logical, intent(out) :: status

      status = .true.

      !FIXME: This is more of an assertion for consistency, might be removed
      if (.NOT. commRank == 0) then
         print *, "Mesh is supposed to be partitioned by master process."
         status = .false.
         return
      end if

      ! ****
      ! Partition Global Mesh
      ! Output:
      !        ProcbyCell
      call PartitionMesh_Metis(Ncpus)

      ! make local mesh
      call LocalMesh_Make

      ! free global mesh
      call GlobalMesh_free

      !comptime_meshmake = MPI_WTIME() - comptime_meshmake

      !do i=1,size(fd)
      !    write(fd(i),'(A,F16.3)') "Computation time of making mesh:  ", &
      !        comptime_meshmake
      !end do

   end subroutine NN_partition_mesh

   subroutine NN_init_warmup(LogFile)

      character(len=*), intent(in) :: LogFile

      ! initialisation petsc/MPI
      call PetscInitialize(PETSC_NULL_CHARACTER, Ierr)
      CHKERRQ(Ierr)

      ! init mpi, communicator/commRank/commSize
      call CommonMPI_init(PETSC_COMM_WORLD)

      call NN_init_output_streams(Logfile)

      comptime_readmesh = MPI_WTIME()

   end subroutine NN_init_warmup

   subroutine NN_init_warmup_and_read_mesh(MeshFile, LogFile)

      character(len=*), intent(in) :: MeshFile, LogFile

      call NN_init_warmup(Logfile)

      ! *** Global Mesh *** !

      if (commRank == 0) then
         call NN_init_read_mesh(MeshFile)
      end if
      call MPI_Barrier(ComPASS_COMM_WORLD, Ierr)

   end subroutine NN_init_warmup_and_read_mesh

subroutine init_visualization(OutputDir)

      character(len=*), intent(in) :: OutputDir

#ifdef _VISU_

      ! initialize visu
      if (trim(MESH_TYPE) == "cartesian-quad") then
         call VisuVTK_VisuTime_Init(MESH_CAR, OutputDir, TimeFinal, output_frequency, &
                                    NbComp, NbPhase, MCP, IndThermique)

      else if (trim(MESH_TYPE) == "hexahedron-quad") then
         call VisuVTK_VisuTime_Init(MESH_HEX, OutputDir, TimeFinal, output_frequency, &
                                    NbComp, NbPhase, MCP, IndThermique)

      else if (trim(MESH_TYPE) == "tetrahedron-triangle") then
         call VisuVTK_VisuTime_Init(MESH_TET, OutputDir, TimeFinal, output_frequency, &
                                    NbComp, NbPhase, MCP, IndThermique)

      else if (trim(MESH_TYPE) == "wedge") then
         call VisuVTK_VisuTime_Init(MESH_WEDGE, OutputDir, TimeFinal, output_frequency, &
                                    NbComp, NbPhase, MCP, IndThermique)
      else
         write (*, *) ""
         write (*, *) "This mesh type is not supported!"
         call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
      end if

      ! vectors that regroup datas (P, T, C, S) from IncCV(:)
      allocate (datavisucell(NbIncPTCSMax*NbCellOwn_Ncpus(commRank + 1)))
      allocate (datavisufrac(NbIncPTCSMax*NbFracOwn_Ncpus(commRank + 1)))
      allocate (datavisunode(NbIncPTCSMax*NbNodeOwn_Ncpus(commRank + 1)))

      ! pressure at well edges, inj/prod
      n = sum(NbEdgebyWellInjLocal(1:NbWellInjOwn_Ncpus(commRank + 1)))
      allocate (datavisuwellinj(n))
      n = sum(NbEdgebyWellProdLocal(1:NbWellProdOwn_Ncpus(commRank + 1)))
      allocate (datavisuwellprod(n))

      if (commRank == 0) then
         do i = 1, size(fd)
            write (fd(i), *) ""
            write (fd(i), *) " *** Warning : visualization vtk of data of wells has not been implemented ***"
         end do
      end if

#endif

    VisuTimeIter = 0

end subroutine init_visualization

subroutine NN_init_phase2(OutputDir)

      character(len=*), intent(in) :: OutputDir
      logical                      :: ok

      if (commRank == 0) then
         do i = 1, size(fd)
            j = fd(i)
            write (j, *) "  NbCell:      ", NbCell
            write (j, *) "  NbFace:      ", NbFace
            write (j, *) "  NbNode:      ", NbNode
            write (j, *) "  NbFrac:      ", NbFrac
            write (j, *) "  NbWellInj    ", NbWellInj
            write (j, *) "  NbWellProd   ", NbWellProd
            write (j, *) "  NbDirNode P: ", NbDirNodeP
#ifdef _THERMIQUE_
            write (j, *) "  NbDirNode T: ", NbDirNodeT
#endif
            write (j, *) "  Ncpus :      ", commSize
            write (j, *) ""
            write (j, *) "Final time: ", TimeFinal/OneSecond
            write (j, *) ""
         end do

         comptime_readmesh = MPI_WTIME() - comptime_readmesh
         do i = 1, size(fd)
            write (fd(i), '(A,F16.3)') "Computation time warm up and reading mesh: ", &
               comptime_readmesh
         end do
      end if

      !comptime_meshmake = MPI_WTIME()

      if (commRank == 0) then
         call NN_partition_mesh(ok)
         if (.NOT. ok) then
            !CHECKME: MPI_Abort is supposed to end all MPI processes
            call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
         end if
      end if
      call MPI_Barrier(ComPASS_COMM_WORLD, Ierr)

      ! *** Global Mesh -> Local Mesh *** !

      ! Main subroutine of module MeshSchema
      ! This module constains all infos of local mesh
      !   1. Send mesh from proc 0 to other procs
      !   2. compute XCellLocal, XFaceLocal
      !   3. compute VolCelllocal, SurfFracLocal
      !   4. compute NbNodeCellmax, NbFracCellMax, NbNodeFaceMax

      call MeshSchema_make

      ! free some tmps in LocalMesh
      if (commRank == 0) then
         call LocalMesh_Free
      end if

      comptime_start = MPI_WTIME() ! total time start
      comptime_total = 0.d0

      ! if(commRank==0) then
      !    print*, sum(NbNodeOwn_Ncpus)
      !    print*, sum(NbFaceOwn_Ncpus)
      !    print*, sum(NbFracOwn_Ncpus)
      !    print*, sum(NbCellOwn_Ncpus)
      !    ! print*, NbFracLocal_Ncpus
      !    ! print*, NbNodeLocal_Ncpus
      !    ! print*, NbCellLocal_Ncpus
      !    ! print*, NbWellInjOwn_Ncpus
      !    ! print*, NbWellInjLocal_Ncpus
      !    ! print*, NbWellProdOwn_Ncpus
      !    ! print*, NbWellProdLocal_Ncpus
      ! end if

      ! if(commRank==1) then
      !    print*, DataWellInjLocal(:)%Radius
      !    print*, DataWellInjLocal(:)%Temperature
      !    print*, DataWellInjLocal(:)%IndWell
      !    print*, DataWellInjLocal(:)%PressionMax
      !    print*, DataWellInjLocal(:)%Flowrate
      !    print*, ""

      !    print*, DataWellProdLocal(:)%Radius
      !    print*, DataWellProdLocal(:)%IndWell
      !    print*, DataWellProdLocal(:)%PressionMin
      !    print*, DataWellProdLocal(:)%Flowrate
      ! end if

      ! if(commRank==0) then
      !    k = 1
      !    do s=NodebyWellInjLocal%Pt(k)+1, NodebyWellInjLocal%Pt(k+1)
      !       print*, DataofNodebyWellInjLocal%Val(s)%WID
      !    end do
      ! end if

      ! *** Numeratation derived from model *** !

      call NumbyContext_make

      ! *** VAG Transmissivity *** !


      call VAGFrac_TransDarcy(PermCellLocal, PermFracLocal)

#ifdef _THERMIQUE_
      call VAGFrac_TransFourier(CondThermalCellLocal, CondThermalFracLocal)
#endif

      deallocate (PermCellLocal)
      deallocate (PermFracLocal)

#ifdef _THERMIQUE_
      deallocate (CondThermalCellLocal)
      deallocate (CondThermalFracLocal)
#endif

      call VAGFrac_VolsDarcy

#ifdef _THERMIQUE_
      call VAGFrac_VolsFourier
#endif

      comptime_part = MPI_WTIME() - comptime_start
      comptime_start = MPI_WTIME()
      if (commRank == 0) then
         do i = 1, size(fd)
            write (fd(i), '(A,F16.3)') "Computation time VAG init.        :", &
               comptime_part
         end do
      end if

      comptime_total = comptime_total + comptime_part

      ! **** allocate structure **** !

      ! unknowns allocate
      call IncCV_allocate
      call DirichletContribution_allocate
      call NeumannContribution_allocate

      ! allcoate Loisthermohydro
      call LoisThermoHydro_allocate

      ! allocate flux
      call Flux_allocate

      ! allocate Residu
      call Residu_allocate

      ! csr sturcture of Jacobian
      ! allocate memory of Jacobian ans Sm
      call Jacobian_StrucJacBigA

      ! csr sturcture of Jacobian after Schur
      ! allocate memory of Jacobian and Sm after Schur
      call Jacobian_StrucJacA

      ! init solver: allocate mat, vector, etc.
      ! call SolvePetsc_Init(KspNiterMax, KspTol)
      call SolvePetsc_cpramgInit(KspNiterMax, KspTol)

      ! allocate KspHistory
      allocate (KspHistory(KspNiterMax + 1))

      ! sync mat create and set value
      call SolvePetsc_SyncMat

      ! allocate increment
      allocate (NewtonIncreNode &
                (NbIncPTCSMax, NbNodeLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreFrac &
                (NbIncPTCSMax, NbFracLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreCell &
                (NbIncPTCSMax, NbCellLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreWellInj &
                (NbWellInjLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreWellProd &
                (NbWellProdLocal_Ncpus(commRank + 1)))

      call init_visualization(OutputDir)

      comptime_part = MPI_WTIME() - comptime_start
      comptime_start = MPI_WTIME()
      if (commRank == 0) then
         do i = 1, size(fd)
            write (fd(i), '(A,F16.3)') "Computation time of allocation:   ", &
               comptime_part
         end do
      end if

      comptime_total = comptime_total + comptime_part

      !! set dir bc values
      !call IncCV_SetDirBCValue
      !
      !! initial value
      !call IncCV_SetInitialValue

      ! init and sort for flash
      call DefFlashWells_allocate

#ifdef _VISU_

      ! visu initial solution

      ! visutime = MPI_WTIME()

      !call IncCV_ToVec( &
      !   datavisucell, datavisufrac, &
      !   datavisuwellinj, datavisuwellprod)
      !
      !call VisuVTK_VisuTime_writedata(0.d0, &
      !      datavisucell, datavisufrac, datavisunode, &
      !                                datavisuwellinj, datavisuwellprod)
      !
      !! if this proc constains at least one well, then write data to file
      !if (NbWellInjOwn_Ncpus(commRank + 1) > 0 .or. &
      !    NbWellProdOwn_Ncpus(commRank + 1) > 0) then
      !
      !   write (output_path, '(A)') trim(OutputDir)//"/wellinfo"
      !   call make_directory(output_path)
      !
      !   write (Wellinfoname, '(A,I0,A)') &
      !      trim(OutputDir)//"/wellinfo/proc_", commRank, ".txt"
      !
      !   open (12, file=Wellinfoname, status="unknown")
      !
      !   write (12, *) "Nb of mpi procs"
      !   write (12, *) commSize
      !
      !   ! nb of well inj of all procs
      !   write (12, *) "Nb of injection wells of all procs"
      !   do i = 1, commSize
      !      write (12, '(I0,A)', advance="no") NbWellInjOwn_Ncpus(i), " "
      !   end do
      !   write (12, *) ""
      !
      !   ! nb of well prod of all procs
      !   write (12, *) "Nb of production wells of all procs"
      !   do i = 1, commSize
      !      write (12, '(I0,A)', advance="no") NbWellProdOwn_Ncpus(i), " "
      !   end do
      !   write (12, *) ""
      !   write (12, *) ""
      !
      !   ! perforation nodes info of well inj
      !   do i = 1, NbWellInjOwn_Ncpus(commRank + 1)
      !
      !      write (12, '(A,I0)') &
      !         "Nb of perforation nodes of own inj well ", i
      !      write (12, '(I0,A)') &
      !         NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i), " "
      !
      !      ! coordinate of perforation nodes
      !      write (12, '(A,I0)') &
      !         "Coordinate of perforation nodes of own inj well ", i
      !      do j = NodebyWellInjLocal%Pt(i) + 1, NodebyWellInjLocal%Pt(i + 1)
      !         numj = NodebyWellInjLocal%Num(j)
      !         write (12, '(F15.6, F15.6, F15.6)') &
      !            XNodeLocal(1, numj), XNodeLocal(2, numj), XNodeLocal(3, numj)
      !      end do
      !   end do
      !
      !   ! perforation nodes info of well prod
      !   do i = 1, NbWellProdOwn_Ncpus(commRank + 1)
      !
      !      write (12, '(A,I0)') &
      !         "Nb of perforation nodes of own prod well ", i
      !      write (12, '(I0,A)') &
      !         NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i), " "
      !
      !      ! coordinate of perforation nodes
      !      write (12, '(A,I0)') &
      !         "Coordinate of perforation nodes of own prod well ", i
      !      do j = NodebyWellProdLocal%Pt(i) + 1, NodebyWellProdLocal%Pt(i + 1)
      !         numj = NodebyWellProdLocal%Num(j)
      !         write (12, '(F15.6, F15.6, F15.6)') &
      !            XNodeLocal(1, numj), XNodeLocal(2, numj), XNodeLocal(3, numj)
      !      end do
      !   end do
      !
      !   close (12)
      !
      !   VisuTimeIter = 1
      !
      !end if

      ! visutime = MPI_WTIME() - visutime
      ! if(commRank==0) then
      !    print*, "visutime is", visutime
      ! end if
#endif

      comptime_total = comptime_total + (MPI_WTIME() - comptime_start)
      comptime_start = MPI_WTIME()

      ! *** Time steps *** !

      TimeCurrent = 0.d0
      TimeOutput = 0.d0

      Delta_t = TimeStepInit

   end subroutine NN_init_phase2

   subroutine NN_init(MeshFile, LogFile, OutputDir)

      character(len=*), intent(in) :: MeshFile, LogFile, OutputDir

      call NN_init_warmup_and_read_mesh(MeshFile, LogFile)
      call NN_init_phase2(OutputDir)

   end subroutine NN_init

   subroutine NN_main_output_visu(TimeIter, OutputDir)

      integer, intent(in) :: TimeIter
      character(len=*), intent(in) :: OutputDir

#ifdef _VISU_

      call IncCV_ToVec( &
        datavisucell, &
        datavisufrac, &
        datavisunode, &
        datavisuwellinj, &
        datavisuwellprod)

      call VisuVTK_VisuTime_writedata( &
        TimeCurrent/OneDay, &
        datavisucell, &
        datavisufrac, &
        datavisunode, &
        datavisuwellinj, &
        datavisuwellprod)

      ! max and min temperature
      Tempmaxloc = -1.d4
      Tempminloc = 1.d4

      do k = 1, NbCellOwn_Ncpus(commRank + 1)
         Tempmaxloc = max(IncCell(k)%Temperature, Tempmaxloc)
         Tempminloc = min(IncCell(k)%Temperature, Tempminloc)
      end do
      do k = 1, NbNodeOwn_Ncpus(commRank + 1)
         Tempmaxloc = max(IncNode(k)%Temperature, Tempmaxloc)
         Tempminloc = min(IncNode(k)%Temperature, Tempminloc)
      end do
      do k = 1, NbFracOwn_Ncpus(commRank + 1)
         Tempmaxloc = max(IncFrac(k)%Temperature, Tempmaxloc)
         Tempminloc = min(IncFrac(k)%Temperature, Tempminloc)
      end do

      call MPI_AllReduce(Tempmaxloc, Tempmax, 1, MPI_DOUBLE, MPI_MAX, ComPASS_COMM_WORLD, Ierr)
      call MPI_AllReduce(Tempminloc, Tempmin, 1, MPI_DOUBLE, MPI_MIN, ComPASS_COMM_WORLD, Ierr)

      ! if(commRank==0) then

      !    ! print*, ""
      !    ! print*, ""
      !    ! print*, Tempmax, Tempmin

      !    do k=1, NbCellOwn_Ncpus(commRank+1)

      !       if( abs(IncCell(k)%Temperature-Tempmax)<1.d-5) then
      !          print*, "cell ", k, XCellLocal(:,k)
      !       end if
      !    end do

      !    do k=1, NbNodeOwn_Ncpus(commRank+1)

      !       if( abs(IncNode(k)%Temperature-Tempmax)<1.d-5) then
      !          print*, "node ", k, XNodeLocal(:,k)
      !       end if
      !    end do

      !    do k=1, NbFracOwn_Ncpus(commRank+1)

      !       if( abs(IncFrac(k)%Temperature-Tempmax)<1.d-5) then
      !          print*, "frac ", k, XCellLocal(:,FracToFaceLocal(k))
      !       end if
      !    end do
      ! end if

      ! if this proc constains at least one well, then write well data to file
      if (NbWellInjOwn_Ncpus(commRank + 1) > 0 .or. &
          NbWellProdOwn_Ncpus(commRank + 1) > 0) then

         ! write well data to file
         write (output_path, '(A,I0)') trim(OutputDir)//"/wellinfo/time_", VisuTimeIter
         call make_directory(output_path)

         write (Wellinfoname, '(A,I0,A,I0,A)') &
            trim(OutputDir)//"/wellinfo/time_", VisuTimeIter, "/proc_", commRank, ".txt"

         open (12, file=Wellinfoname, status="unknown")

         ! time step info
         write (12, *) "TimeStep"
         write (12, *) TimeIter
         write (12, *) "Time"
         write (12, '(F16.5)') TimeCurrent/OneDay

         ! data of perforation nodes of inj well
         do i = 1, NbWellInjOwn_Ncpus(commRank + 1)

            head = NodebyWellInjLocal%Pt(i + 1) ! head

            write (12, '(A,I0)') &
               "Nb of perforation nodes of own inj well ", i
            write (12, '(I0,A)') &
               NodebyWellInjLocal%Pt(i + 1) - NodebyWellInjLocal%Pt(i), " "

            ! data of head perforation node
            write (12, '(A,I0)') &
               "Data of head perforation node of own inj well ", i

            write (12, '(ES17.7)', advance='no') IncPressionWellInj(i) ! pressure
            write (12, '(ES17.7)', advance='no') PerfoWellInj(head)%Temperature ! temperature
! FIXME: headmolarFluxInj is not defined in all configuration files
!             write(12,'(ES17.7)',advance='no') headmolarFluxInj(i) ! head molar flux
            write (12, *) ""
         end do

         ! data of perforation nodes of prod well
         do i = 1, NbWellProdOwn_Ncpus(commRank + 1)

            write (12, '(A,I0)') &
               "Nb of perforation nodes of own prod well ", i
            write (12, '(I0,A)') &
               NodebyWellProdLocal%Pt(i + 1) - NodebyWellProdLocal%Pt(i), " "

            ! data of head perforation node
            write (12, '(A,I0,A)') &
               "Data of head perforation node of own prod well ", i, &
               ":  Pressure / Temperature / cumulFluxmolar / cumulFluxEnergy"
            head = NodebyWellProdLocal%Pt(i + 1)

            do j = NodebyWellProdLocal%Pt(i) + 1, NodebyWellProdLocal%Pt(i + 1)

               ! write(12,'(ES17.7)',advance='no') IncPressionWellProd(i)     ! pressure
               write (12, '(ES17.7)', advance='no') PerfoWellProd(j)%Pression ! pressure
               write (12, '(ES17.7)', advance='no') PerfoWellProd(j)%Temperature ! temperature
               write (12, '(ES17.7)', advance='no') summolarFluxProd(:, j) ! head molar flux
               write (12, '(ES17.7)', advance='no') sumnrjFluxProd(j) ! head energy flux
               write (12, *) ""
            end do
         end do

         ! max and min temperature
         write (12, '(A)') "max/min temperature"
         write (12, '(ES17.7)', advance='no') Tempmax
         write (12, '(ES17.7)', advance='no') Tempmin
         write (12, *) ""

         close (12)

         VisuTimeIter = VisuTimeIter + 1

      end if
#endif

   end subroutine NN_main_output_visu

   subroutine NN_main_checkpoint(OutputDir)

      character(len=*), intent(in) :: OutputDir

#ifdef _HDF5_
      call IncCV_WriteSolToFile(OutputDir, &
                                TimeIter, TimeCurrent, Delta_t, TimeOutput, &
                                NewtonNiterTotal, NewtonNbFailure, KspNiterTotal, KspNbFailure, &
                                comptime_total, comptime_timestep)
#endif

   end subroutine NN_main_checkpoint

   subroutine NN_main_make_timestep(initial_time_step)

      real(c_double), optional, intent(in) :: initial_time_step

      if(present(initial_time_step)) Delta_t = max(0.d0, initial_time_step)

      ! init start time
      comptime_start = MPI_WTIME()

      NewtonConv = .false.
      KspConv = .false.

      ! save current status, copy Inc* to Inc*PrevisousTimeStep
      call IncCV_SaveIncPreviousTimeStep

      call IncCVWells_PressureDropWellInj ! compute PerfoWellInj%Pression
      call IncCVWells_PressureDropWellProd ! compute PerfoWellDrop%Pression

      ! *** Newton iterations *** !

      ! if Jacobian Ksp solver doesn't converge
      !   restart a new Newton with a smaller time step
      ! else if converge,
      !   KspConv will is set as .true.

      ! if Newton doesn't converge
      !   restart a new Newton with a smaller time step
      ! else if converge
      !   NewtonConv and KspConv are set as .true.
      do while ((NewtonConv .eqv. .false.) .or. (KspConv .eqv. .false.))

         ! init as false
         NewtonConv = .false.
         KspConv = .false.

         KspNiterTimeStep = 0 ! total nb of Ksp of time step

         ! Newton iteration
         ! if Newton converge
         !   NewtonConv will be set as .true.
         !   KspConv will also be set as .true.
         do NewtonIter = 1, NewtonNiterMax

            ! Copy Dir boundary values to Inc
            call DirichletContribution_update

            ! compute pressure of perforations with well pressure
!           IncPressionWellInj(:) = 2.d7
            call IncCVWells_PressureDropWellInj
            call IncCVWells_PressureDropWellProd

            ! LoisThermohydro
            call LoisThermoHydro_compute

            ! compute flux cell/frac
            call Flux_DarcyFlux_Cell
            call Flux_DarcyFlux_Frac

#ifdef _THERMIQUE_
            call Flux_FourierFlux_Cell
            call Flux_FourierFlux_Frac
#endif

            ! compute Residu
            call Residu_compute(Delta_t, NewtonIter)

           !  if(commRank==1) then
            !    ! print*, ResiduNode
            !    ! print*, ResiduCell
            !    ! print*, ResiduWellInj
           !  end if

            ! write(*,*) ""
            ! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)

            ! test Newton converge
            call Residu_RelativeNorm(NewtonIter, Delta_t, &
                                     NewtonResNormRel(NewtonIter), NewtonResConvInit, NewtonResClosInit)

            if (commRank == 0) then
               if (NewtonIter == 1) then
                  do i = 1, size(fd)
                     j = fd(i)
                     write (j, *) ""
                     write (j, '(A)', advance='no') "     *Residu init conv:  "
                     do k = 1, NbCompThermique
                        write (j, '(A,ES12.5)', advance='no') "  ", NewtonResConvInit(k)
                     end do
                     write (j, *) ""
                     write (j, '(A,ES12.5)') "     *Residu init clos:    ", NewtonResClosInit
                     write (j, *) ""
                  end do
               end if

               do i = 1, size(fd)
                  j = fd(i)
                  write (j, '(A,I4)', advance="no") "     *Newton iter:", NewtonIter
                  write (j, '(A,E18.10)', advance="no") "      Newton res norm:", NewtonResNormRel(NewtonIter)
               end do
            end if

            if (NewtonResNormRel(NewtonIter) < NewtonTol) then

               NewtonConv = .true.
               KspConv = .true.
               exit
            end if

            ! Jacobian and second member
            !   inputs:  Residu
            !   outputs: JacA, Sm
            call Jacobian_ComputeJacSm(Delta_t)

            
            ! set values of matrix, vector, Ksp solver
            !   inputs : JacA, Sm
            !   outputs : solver Petsc

            call SolvePetsc_SetUp

            ! write(*,*) ""
            ! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)

            ! solve, return nb of iterations
            ! if KspNiter<0, not converge
            call SolvePetsc_KspSolve(KspNiter, KspHistory)

            if (commRank == 0) then
               do i = 1, size(fd)
                  write (fd(i), '(A,I5)', advance="no") "    Nb of Ksp iter:  ", KspNiter
               end do
            end if

            if (KspNiter < 0) then ! not converge

               ! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)

               KspConv = .false. ! solver not converge
               KspNbFailure = KspNbFailure + 1

               ! restart a new Newton with a smaller time step
               Delta_t = Delta_t*0.5d0

               ! load status
               call IncCV_LoadIncPreviousTimeStep
               call IncCVWells_PressureDropWellInj ! compute PerfoWellInj%Pression
               call IncCVWells_PressureDropWellProd ! compute PerfoWellDrop%Pression

               ! print residu if Ksp failure
               if (commRank == 0) then
                  do i = 1, size(fd)
                     j = fd(i)
                     write (j, *) ""
                     do k = 1, KspNiterMax
                        write (j, '(A,I4,A,ES18.10)') &
                           "             Ksp iter:   ", k, "    res:  ", KspHistory(k)
                     end do

                     write (j, *) ""
                     write (j, '(A)', advance="no") &
                        "   -- Restart a Newton with a smaller time step (Ksp does not converge): "
                     write (j, *) Delta_t/OneDay
                  end do
               end if

               exit
            else
               KspConv = .true.
               KspNiterTimeStep = KspNiterTimeStep + KspNiter

               call SolvePetsc_Sync ! sync for ghost

               ! Get increment prim of node, frac and wells
               ! Increment of prim is NewtonIncre(1:NbComp+IndThermique)
               ! Get values from vector of petsc/trilinos
               call SolvePetsc_GetSolNodeFracWell( &
                  NewtonIncreNode, NewtonIncreFrac, &
                  NewtonIncreWellInj, NewtonIncreWellProd)

               ! print*, NewtonIncreWellInj

               ! Compute increment prim of cell
               ! Inverse Schur
               call Jacobian_GetSolCell(NewtonIncreNode, &
                                        NewtonIncreFrac, NewtonIncreCell)

               ! Compute incremment increment of secd using increment prim
               !   Input : increment prim is NewtonIncre(1:NbComp+IndThermique)
               !   Output: complete increment = (prim and secd, acc)
               call LoisThermoHydro_PrimToSecd( &
                  NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell)

               ! compute Newton relaxation
               call IncCVReservoir_NewtonRelax( &
                  NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, NewtonRelax)

               ! ???
               ! NewtonRelax = 1.d0

               if (commRank == 0) then
                  do i = 1, size(fd)
                     write (fd(i), '(A,F12.7)') "    Relaxation:", NewtonRelax
                  end do
               end if

               ! update Inc with Increment
               call IncCV_NewtonIncrement( &
                  NewtonIncreNode, NewtonIncreFrac, NewtonIncreCell, &
                  NewtonIncreWellInj, NewtonIncreWellProd, NewtonRelax)

               call DirichletContribution_update

               ! call IncCV_ToVec( &
               !      dataviscuell, datavisufrac, datavisunode &
               !      datavisuwellinj, datavisuwellprod)

               ! call VisuVTK_VisuTime_writedata(TimeCurrent/OneDay, &
               !      datavisucell, datavisufrac, datavisunode, &
               !      datavisuwellinj, datavisuwellprod)

               ! Flash
               call DefFlash_Flash

            end if ! end if converge/not converge

         end do ! end of : NewtonIter=1, NewtonNiterMax

         if (NewtonConv .eqv. .true.) then

            NewtonNiterTotal = NewtonNiterTotal + NewtonIter
            KspNiterTotal = KspNiterTotal + KspNiterTimeStep

         else if ((NewtonConv .eqv. .false.) .and. &
                  (KspConv .eqv. .true.)) then

            NewtonNbFailure = NewtonNbFailure + 1

            ! restart a new Newton with a smaller time step
            Delta_t = Delta_t*0.5d0

            if (commRank == 0) then
               write (*, *) ""
               write (*, '(A)', advance="no") &
                  "   -- Restart a Newton with a smaller time step (Newton does not converge): "
               write (*, *) Delta_t/OneSecond

               write (11, *) ""
               write (11, '(A)', advance="no") &
                  "   -- Restart a Newton with a smaller time step (Newton does not converge): "
               write (11, *) Delta_t/OneSecond

               
            end if

            ! load status
            call IncCV_LoadIncPreviousTimeStep
            call IncCVWells_PressureDropWellInj ! compute PerfoWellInj%Pression
            call IncCVWells_PressureDropWellProd ! compute PerfoWellDrop%Pression

            ! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
         end if

      end do ! end of do while( (NewtonConv .eqv. .false.) .or. (KspConv .eqv. .false.))

      call DefFlashWells_TimeFlash

      TimeCurrent = TimeCurrent + Delta_t      

      ! FiXME: What is the policy for time step management
      ! compute Delta_t for the next time step
      call IncCV_ComputeTimeStep(Delta_t, TimeCurrent)
      ! ???
      ! Delta_t = TimeStepInit

      ! total computation time and computation time of this time step
      comptime_timestep = MPI_WTIME() - comptime_start
      comptime_total = comptime_total + comptime_timestep

   end subroutine NN_main_make_timestep

   subroutine NN_main(TimeIter, OutputDir)

      integer, intent(inout) :: TimeIter
      character(len=*), intent(in) :: OutputDir

#ifdef _HDF5_
      if (TimeIter > 0) then
         call IncCV_ReadSolFromFile(OutputDir, &
                                    TimeIter, TimeCurrent, Delta_t, TimeOutput, &
                                    NewtonNiterTotal, NewtonNbFailure, KspNiterTotal, KspNbFailure, &
                                    comptime_total, comptime_timestep)
      end if
#endif

      do while (TimeCurrent < (TimeFinal + eps))

         TimeIter = TimeIter + 1

         if (commRank == 0) then
            do i = 1, size(fd)
               j = fd(i)
               write (j, *) ""
               write (j, *) ""
               write (j, '(A,I0)') "Time Step: ", TimeIter
               write (j, '(A,F16.5)') "Time at previous time step: ", TimeCurrent/OneSecond, "seconds"

               write (j, *)
               write (j, '(A)', advance="no") "   -- Initial time step: "
               write (j, *) Delta_t/OneSecond
            end do
         end if

         call NN_main_make_timestep

         ! checkpoint and visu
         if (TimeCurrent > TimeOutput) then

            call NN_main_output_visu(TimeIter, OutputDir)

            call NN_main_checkpoint(OutputDir)

            TimeOutput = NN_ceiling(TimeCurrent / output_frequency) * output_frequency

         end if

         if (commRank == 0) then
            call NN_main_summarize_timestep
         end if

      end do ! end of time steps

   end subroutine NN_main


   ! personal ceiling function
   ! ceiling(x) returns an integer which unbound when x is huge
   FUNCTION NN_ceiling(x)
     DOUBLE PRECISION, INTENT(IN) :: x

     DOUBLE PRECISION :: NN_ceiling

     NN_ceiling = AINT(x)
     NN_ceiling = NN_ceiling + MERGE(1, 0, x > 0 .AND. x /= NN_ceiling)
   END FUNCTION NN_ceiling


   subroutine NN_main_summarize_timestep()

      do i = 1, size(fd)
         j = fd(i)
         write (j, *)
         write (j, *)
         write (j, '(A,I0)') "     -Total nb of Newton iters:     ", NewtonNiterTotal
         write (j, '(A,I0)') "     -Total nb of Ksp iters:        ", KspNiterTotal
         write (j, '(A,I0)') "     -Total nb of Newton failures:  ", NewtonNbFailure
         write (j, '(A,I0)') "     -Total nb of Ksp failures:     ", KspNbFailure

         write (j, *)
         write (j, '(A,F15.3)') "     -Total Computation time:              ", comptime_total
         write (j, '(A,F15.3)') "     -Computation time of this time step:  ", comptime_timestep
      end do

   end subroutine NN_main_summarize_timestep

   subroutine NN_finalize()

      ! ! compute errors
      ! errlocal_cell_L1 = 0.d0
      ! errlocal_cell_L2 = 0.d0
      ! errlocal_cell_Linf = 0.d0

      ! errlocal_frac_L1 = 0.d0
      ! errlocal_frac_L2 = 0.d0
      ! errlocal_frac_Linf = 0.d0

      ! do k=1, NbCellOwn_Ncpus(commRank+1)

      !    call DefModel_AnalyticSol(XCellLocal(:,k),sol)

      !    errlocal_cell_L2   = errlocal_cell_L2 + (IncCell(k)%Pression - sol)**2 * VolCellLocal(k)
      !    errlocal_cell_L1   = errlocal_cell_L1 + abs(IncCell(k)%Pression - sol) * VolCellLocal(k)
      !    errlocal_cell_Linf = max(abs(IncCell(k)%Pression - sol), errlocal_cell_Linf)
      ! end do

      ! do k=1, NbFracOwn_Ncpus(commRank+1)

      !    call DefModel_AnalyticSol(XFaceLocal(:,FracToFaceLocal(k)),sol)

      !    errlocal_frac_L2   = errlocal_frac_L2 + (IncFrac(k)%Pression - sol)**2 * SurfFracLocal(k)
      !    errlocal_frac_L1   = errlocal_frac_L1 + abs(IncFrac(k)%Pression - sol) * SurfFracLocal(k)
      !    errlocal_frac_Linf = max(abs(IncFrac(k)%Pression - sol), errlocal_frac_Linf)
      ! end do

      ! call MPI_Reduce(errlocal_cell_L2, err_cell_L2, 1, MPI_DOUBLE, MPI_SUM, 0, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Reduce(errlocal_cell_L1, err_cell_L1, 1, MPI_DOUBLE, MPI_SUM, 0, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Reduce(errlocal_cell_Linf, err_cell_Linf, 1, MPI_DOUBLE, MPI_MAX, 0, ComPASS_COMM_WORLD, Ierr)

      ! call MPI_Reduce(errlocal_frac_L2, err_frac_L2, 1, MPI_DOUBLE, MPI_SUM, 0, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Reduce(errlocal_frac_L1, err_frac_L1, 1, MPI_DOUBLE, MPI_SUM, 0, ComPASS_COMM_WORLD, Ierr)
      ! call MPI_Reduce(errlocal_frac_Linf, err_frac_Linf, 1, MPI_DOUBLE, MPI_MAX, 0, ComPASS_COMM_WORLD, Ierr)

      ! if(commRank==0) then
      !    print*, ""
      !    print*, "errl2cell[mm] =", sqrt(err_cell_L2) / (2000.d0)**3
      !    print*, "errl1cell[mm] =", err_cell_L1 / (2000.d0)**3
      !    print*, "errlinfcell[mm] =", err_cell_Linf

      !    print*, ""
      !    print*, "errl2frac[mm] =", sqrt(err_frac_L2) / (2000.d0)**2
      !    print*, "errl1frac[mm] =", err_frac_L1 / (2000.d0)**2
      !    print*, "errlinffrac[mm] =", err_frac_Linf
      ! end if

      ! *** Report *** !
      if (commRank == 0) then
         do i = 1, size(fd)
            j = fd(i)
            write (j, *) ""
            write (j, *) ""
            write (j, *) "Final Report"

            write (j, *) ""
            write (j, *) "    *Final time:  ", TimeFinal/OneSecond

            write (j, *) ""
            write (j, *) "    *Mesh:"
            write (j, '(A,I0)') "        -Nb of cells:  ", NbCell
            write (j, '(A,I0)') "        -Nb of faces:  ", NbFace
            write (j, '(A,I0)') "        -Nb of nodes:  ", NbNode
            write (j, '(A,I0)') "        -Nb of fracs:  ", NbFrac

            write (j, *) ""
            write (j, *) "    *Newton/Ksp Iterations:"
            write (j, '(A,I0)') "        -Total nb of Newton iters:  ", NewtonNiterTotal
            write (j, '(A,I0)') "        -Total nb of Ksp iters:  ", KspNiterTotal

            write (j, *) ""
            write (j, *) "    *Newton/Ksp Failures:"
            write (j, '(A,I0)') "        -Total nb of Newton failures:  ", NewtonNbFailure
            write (j, '(A,I0)') "        -Total nb of Ksp failures:  ", KspNbFailure

            write (j, *) ""
            write (j, *) "    *Total simulation time:  ", comptime_total
         end do
      end if

      ! *** Free *** !
#ifdef _VISU_
      call VisuVTK_VisuTime_pvdwriter ! write pvd file
      call VisuVTK_VisuTime_free
      deallocate (datavisucell)
      deallocate (datavisufrac)
      deallocate (datavisunode)
      deallocate (datavisuwellinj)
      deallocate (datavisuwellprod)
#endif

      deallocate (NewtonIncreNode)
      deallocate (NewtonIncreFrac)
      deallocate (NewtonIncreCell)
      deallocate (NewtonIncreWellInj)
      deallocate (NewtonIncreWellProd)
      deallocate (KspHistory)
      if (allocated(fd)) then
         deallocate (fd)
      end if

      call SolvePetsc_free
      call Jacobian_free
      call Residu_free
      call Flux_free
      call VAGFrac_free
      call LoisThermoHydro_free
      call IncCV_free
      call DirichletContribution_free
      call NeumannContribution_free
      ! call DefFlash_free
      call MeshSchema_free

      call MPI_Barrier(ComPASS_COMM_WORLD, Ierr)

      if (commRank == 0) then
         close (11)
         write (*, *) ""
         write (*, *) "NN Finalize"
      end if
      call PetscFinalize(Ierr)

   end subroutine NN_finalize

! FIXME: Spare comments? Remove them?

! if(commRank==0) then

!    do i=1, NbPhasePresente_ctx(IncCell(1)%ic)
!       iph = NumPhasePresente_ctx(i,IncCell(1)%ic)

! print*, NumIncPTCSPrimCell(:,1)
! print*, NumIncPTCSecondCell(:,1)

! write(*,'(ES22.14)'), DensiteMassiqueCell(i,1)
! print*, ""

! do j=1, NbIncPTCSPrimMax
!    write(*,'(ES22.14)'), divDensiteMassiqueCell(j,i,1)
! end do
! print*, ""

! write(*,'(ES22.14)'), SmDensiteMassiqueCell(i,1)
! print*, ""

!    ! ! *** !

!    do j=1, NbIncPTCSPrimMax
!       write(*,'(ES22.14)'), divDensitemolaireKrViscoCell(j,i,1)
!    end do
!    print*, ""
!    write(*,'(ES22.14)'), SmDensitemolaireKrViscoCell(i,1)
!    print*, ""

!    ! ! *** !

! do j=1, NbIncPTCSPrimMax
!    write(*,'(ES22.14)'), divDensitemolaireKrViscoEnthalpieCell(j,i,1)
! end do
! print*, ""
! write(*,'(ES22.14)'), SmDensitemolaireKrViscoEnthalpieCell(i,1)
! print*, ""

!        print*, NumIncPTCSecondCell(:,1)

! do j=1, NbComp
!    if(MCP(j,iph)==1) then
!       do s=1, NbIncPTCSPrim_ctx(IncCell(1)%ic)
!          write(*,'(ES22.14)'), divDensitemolaireKrViscoCompCell(s,j,i,1)
!       end do
!       print*, ""
!       write(*,'(ES22.14)'), SmDensitemolaireKrViscoCompCell(j,i,1)
!    end if
!    print*, ""
! end do

! do j=1, NbComp
!    if(MCP(j,iph)==1) then
!       do s=1, NbIncPTCSPrim_ctx(IncCell(1)%ic)
!          write(*,'(ES22.14)'), divDensitemolaireSatCompCell(s,j,i,1)
!       end do
!       print*, ""
!       write(*,'(ES22.14)'), SmDensitemolaireSatCompCell(j,i,1)
!    end if
!    print*, ""
! end do

! do j=1, NbIncPTCSPrimMax
!    write(*,'(ES22.14)'), divDensitemolaireKrViscoEnthalpieCell(j,i,1)
! end do
! print*, ""
! write(*,'(ES22.14)'), SmDensitemolaireKrViscoEnthalpieCell(i,1)
! print*, ""

!    ! *** !

!    do j=1, NbIncPTCSPrimMax
!       write(*,'(ES22.14)'), divDensitemolaireEnergieInterneSatCell(j,i,1)
!    end do
!    print*, ""
!    write(*,'(ES22.14)'), SmDensitemolaireEnergieInterneSatCell(i,1)
!    print*, ""
!    end do

! end if

! if(commRank==1) then
! do k=1, NbCellLocal_Ncpus(commRank+1) ! loop of cell

!    ! number of nodes/fracs in cell k
!    NbNodeCell = NodebyCellLocal%Pt(k+1) - NodebyCellLocal%Pt(k)
!    NbFracCell = FracbyCellLocal%Pt(k+1) - FracbyCellLocal%Pt(k)

!    print*, FluxDarcyKI(1,1:NbNodeCell+NbFracCell, k)
!    print*, FluxFourierKI(1:NbNodeCell+NbFracCell, k)
! end do

! do k=1, NbFracLocal_Ncpus(commRank+1)

!    fk = FracToFaceLocal(k) ! fk is face number

!    NbNodeFrac = NodebyFaceLocal%Pt(fk+1) - NodebyFaceLocal%Pt(fk)

!    !print*, FluxDarcyFI(NumPhasePresente_ctx(1,IncFrac(k)%ic ),1:NbNodeFrac, k)
!    print*, FluxFourierFI(1:NbNodeFrac, k)
! end do

! end if

! call Residu_Accmolaire

! allocate( x12( NbComp+IndThermique, NbCellLocal_Ncpus(commRank+1) &
!      + NbNodeLocal_Ncpus(commRank+1)+NbFracLocal_Ncpus(commRank+1) ))

! x12(:,:) = 0.d0

! do i=1, NbNodeLocal_Ncpus(commRank+1)+NbFracLocal_Ncpus(commRank+1)

! ! do i=1, NbCellLocal_Ncpus(commRank+1) &
! !      + NbNodeLocal_Ncpus(commRank+1)+NbFracLocal_Ncpus(commRank+1)

!    do j=1, NbComp+IndThermique
!       x12(j,i) = 10.d0 * dble(i)! + dble(j)
!    end do
! end do

! if(commRank==1) then
!    do i=1, NbCellLocal_Ncpus(commRank+1) &
!         + NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)

! print*, x12(:,1)
!    end do
!    print*, ""
! end if

! bigSm(:,:) = 0.d0

! do i=1, NbCellLocal_Ncpus(commRank+1) &
!      + NbNodeOwn_Ncpus(commRank+1) + NbFracOwn_Ncpus(commRank+1)

!    do j=JacBigA%Pt(i)+1, JacBigA%Pt(i+1)
!       colj = JacBigA%Num(j)

!       do i1=1, 5
!          do j1= 1, 5
!             bigSm(i1,i) = bigSm(i1,i) + &
!                  JacBigA%Val(j1,i1,j) * x12(j1,colj)
!          end do
!       end do

!    end do
! end do

! if(commRank==1) then
!    do i=1, NbCellLocal_Ncpus(commRank+1) &
!         + NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)

!       print*, bigSm(:,i)
!    end do
!    print*, ""
! end if

! do i=1, NbCellLocal_Ncpus(commRank+1)+ &
!      NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)

!    do j=JacbigA%Pt(i)+1, JacbigA%Pt(i+1)
!       colj = JacbigA%Num(j)

!       do i1=1, 5
!          do j1 = 1, 5
!             bigSm(i1,i) = bigSm(i1,i) - &
!                  JacbigA%Val(j1,i1,j) * x12(j1,colj)
!          end do
!       end do
!    end do
! end do

! if(commRank==1) then
!    do i=1, NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)
!       print*, i, bigSm(1:5,i)
!    end do
! end if

! if(commRank==1) then
!    print*, bigSm(:,1)
!    print*, Sm(:,1)
! end if

! do i=1, NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)

!    tmp(:) = 0.d0
!    do j=JacA%Pt(i)+1, JacA%Pt(i+1)
!       colj = JacA%Num(j)

!       do i1=1, 5
!          do j1=1, 5
!             tmp(i1) = tmp(i1) +  &
!                  JacA%Val(j1,i1,j) * x12(j1,colj)
!          end do
!       end do
!    end do

!    if(commRank==1) then
!       print*, ""
!       print*, (tmp(:)-Sm(:,i))!/Sm(:,i)
!    end if
! end do

! if(commRank==1) then
!    do i=1, 1! NbNodeOwn_Ncpus(commRank+1)+NbFracOwn_Ncpus(commRank+1)
!      ! print*, i, Sm(1:5, i)
!    end do
! end if

! if(commRank==0) then

!    ! do i=JacA%Pt(j)+1, JacBigA%Pt(j+1)
!    !    if( JacBigA%Num(i)==colk) then
!    !       nz = i
!    !       exit
!    !    end if
!    ! end do

!    ! do i=1, 5
!    !    print*, JacA%Pt(i+1)-JacA%Pt(i)
!    ! end do

!    ! do i=1, NbComp+1
!    !    do j=1, NbIncPTCSPrim_ctx(IncCell(1)%ic)
!    !       write(*,'(ES22.13)') JacA%Val(j,i,1)
!    !    end do
!    !    print*, ""
!    ! end do
! end if

! if(commRank==1) then
!    do i=21, 25
!       call MatGetRow(A_p, i, ncols, cols, values, Ierr)

!       print*, values(61:65)
!       ! do j=1, ncols
!       !    if(i==cols(j)) then
!       !       print*, values(j)
!       !    end if
!       ! end do

!       call MatRestoreRow(A_p, i, ncols, cols, values, Ierr)
!    end do
! end if

! if(commRank==1) then

!    do i=1,12
!       print*, IncNode(i)%Pression
!       print*, IncNode(i)%Temperature
!       print*, IncNode(i)%Comp(1,:)
!       print*, IncNode(i)%Saturation(:)
!       print*, ""
!    end do

!    do i=1,2
!       print*, IncCell(i)%Pression
!       print*, IncCell(i)%Temperature
!       print*, IncCell(i)%Comp(1,:)
!       print*, IncCell(i)%Saturation(:)
!       print*, ""
!    end do

!    print*, IncFrac(1)%Pression
!    print*, IncFrac(1)%Temperature
!    print*, IncFrac(1)%Comp(1,:)
!    print*, IncFrac(1)%Saturation(:)
!    print*, ""
! end if

! VolDarcyCell(:) = 0.5d3
! PoroVolDarcyCell(:) = 0.5d3
! PoroVolFourierCell(:) = 0.5d3
! Poro_1volFourierCell(:) = 0.5d3

! VolDarcyFrac(:) = 0.5d3
! PoroVolDarcyFrac(:) = 0.5d3
! PoroVolFourierFrac(:) = 0.5d3
! Poro_1volFourierFrac(:) = 0.5d3

! VolDarcyNode(:) = 0.5d3
! PoroVolDarcyNode(:) = 0.5d3
! PoroVolFourierNode(:) = 0.5d3
! Poro_1volFourierNode(:) = 0.5d3

! do i=1, NbCellLocal_Ncpus(commRank+1)
!    TkLocal_Darcy(i)%pt(:,:) = 0.d0
!    TkLocal_Fourier(i)%pt(:,:) = 0.d0
! end do

! do i=1, NbFracLocal_Ncpus(commRank+1)
!    TkFracLocal_Darcy(i)%pt(:,:) = 0.d0
!    TkFracLocal_Fourier(i)%pt(:,:) = 0.d0
! end do

! do k=1, NbCellLocal_Ncpus(commRank+1)

!    j = NodebyCellLocal%pt(k+1)-NodebyCellLocal%pt(k) &
!         + FracbyCellLocal%pt(k+1)-FracbyCellLocal%Pt(k)

!    ! do i=1, j
!    !    TkLocal_Darcy(k)%pt(i,i) = 1.d-11
!    !    TkLocal_Fourier(k)%pt(i,i) = 1.d2
!    ! end do
! end do

! do i=1, NbFracLocal_Ncpus(commRank+1)
!    fi = FracToFaceLocal(i)
!    j = NodebyFaceLocal%Pt(fi+1)-NodebyFaceLocal%Pt(fi)

!    do k=1, j
!       ! TkFracLocal_Darcy(i)%pt(k,k) = 1.d-9
!       TkFracLocal_Fourier(i)%pt(k,k) = 1.d4
!    end do
! end do

! if(commRank==1) then
!    print*, PoroVolDarcyCell(:)
!    print*, PoroVolFourierCell(:)
!    print*, Poro_1volFourierCell(:)

!    print*, PoroVolDarcyFrac
!    print*, PoroVolDarcyNode
! end if

! VolDarcyCell(:) = 0.d0
! PoroVolDarcyCell(:) = 0.d0
! PoroVolFourierCell(:) = 0.d0
! Poro_1volFourierCell(:) = 0.d0

! VolDarcyFrac(:) = 0.d0
! PoroVolDarcyFrac(:) = 0.d0
! PoroVolFourierFrac(:) = 0.d0
! Poro_1volFourierFrac(:) = 0.d0

! VolDarcyNode(:) = 0.d0
! PoroVolDarcyNode(:) = 0.d0
! PoroVolFourierNode(:) = 0.d0
! Poro_1volFourierNode(:) = 0.d0

! open(unit=12, file="res.log", status="unknown")
! do i=1, NewtonIter
!    write(12,*) NewtonResNormRel(i)
! end do
! close(12)

! if(commRank==0) then
!    ! print*, XCellLocal(:,1)
!    ! print*, "T cell", IncCell(1)%Temperature
!    ! print*, "S cell", IncCell(1)%Saturation
!    ! print*, "P cell", IncCell(1)%Pression
!    print*, ""
!    print*, "T frac", IncFrac(1)%Temperature
!    print*, "S frac", IncFrac(1)%Saturation
!    print*, "P frac", IncFrac(1)%Pression
! end if

! print*, ""
! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)

! do i=1, NbNodeLocal_Ncpus(commRank+1)
!    if(IncNode(i)%Temperature>451.d0) then
!       print*, "T>450 node:", commRank, i, IncNode(i)%Temperature, XNodeLocal(:,i)
!    end if
! end do

! do i=1, NbFracLocal_Ncpus(commRank+1)
!    if(IncFrac(i)%Temperature>451.d0) then
!       print*, "T>450 frac:", commRank, i, IncFrac(i)%Temperature, XFaceLocal(:,FracToFaceLocal(i))
!    end if
! end do

! do i=1, NbCellLocal_Ncpus(commRank+1)
!    if(IncCell(i)%Temperature>451.d0) then
!       print*, "T>450 cell:", commRank, i, IncCell(i)%Temperature, XCellLocal(:,i)
!    end if
! end do

end module NN
