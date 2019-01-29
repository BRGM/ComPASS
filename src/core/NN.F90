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
   ! Residu relative norm from first Newton iteration
   type(CTVector) :: NewtonResConvInit
   real(c_double) :: NewtonResClosInit

   ! ksp variables
   double precision, allocatable, dimension(:) :: KspHistory

   logical :: KspConv
   integer :: KspNiter, KspNiterTimeStep
   integer :: KspNiterTotal = 0
   integer :: KspNbFailure = 0

   ! Newton increment (NbInc,NbNodelocal)
   real(c_double), dimension(:, :), allocatable :: &
      NewtonIncreNode, &
      NewtonIncreFrac, &
      NewtonIncreCell

   ! as well have only one unknown (pressure), only a vector (NbWellLocal)
   real(c_double), dimension(:), allocatable :: &
      NewtonIncreWellInj, &
      NewtonIncreWellProd

   !  ! Perm
   !  double precision, dimension(:,:,:), allocatable :: PermCellLocal
   !  double precision, dimension(:), allocatable :: PermFracLocal

   double precision :: visutime

   double precision :: f, CC(NbComp), SS(NbPhase), dPf, dTf, dCf(NbComp), dSf(NbPhase)

   ! ! ********************************** ! !

   public :: &
      NN_main_make_timestep, &
      NN_main_summarize_timestep, &
      NN_flash_all_control_volumes, &
      NN_finalize
   
   private :: &
       NN_flash_control_volumes

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

   !subroutine NN_init_read_mesh(MeshFile)
   !
   !   character(len=*), intent(in) :: MeshFile
   !
   !   !FIXME: This is more of an assertion for consistency, might be removed
   !   if (.NOT. commRank == 0) then
   !      print *, "Mesh is supposed to be read by master process."
   !      !CHECKME: MPI_Abort is supposed to end all MPI processes
   !      call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
   !   end if
   !
   !   inquire (FILE=MeshFile, EXIST=file_exists)
   !   if (file_exists .eqv. .false.) then
   !      print *, " "
   !      print *, "Mesh does not exist   ", MeshFile
   !      print *, " "
   !      !CHECKME: MPI_Abort is supposed to end all MPI processes
   !      call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
   !   end if
   !   print *, "Mesh read from file: ", MeshFile
   !
   !   ! Read Global Mesh
   !   call GlobalMesh_Make_read_file(MeshFile)
   !
   !   call GlobalMesh_Make_post_read()
   !
   !   call DefWell_Make_SetDataWell(NbWellInj, NbWellProd)
   !   call DefWell_Make_ComputeWellIndex( &
   !      NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
   !      PermCell, PermFrac)
   !
   !end subroutine NN_init_read_mesh

   !subroutine NN_init_build_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)
   !
   !   real(kind=c_double), intent(in)  :: Ox, Oy, Oz
   !   real(kind=c_double), intent(in)  :: lx, ly, lz
   !   integer(kind=c_int), intent(in)  :: nx, ny, nz
   !
   !   !FIXME: This is more of an assertion for consistency, might be removed
   !   if (.NOT. commRank == 0) then
   !      print *, "Mesh is supposed to be built by master process."
   !      !CHECKME: MPI_Abort is supposed to end all MPI processes
   !      call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
   !   end if
   !
   !   call GlobalMesh_Build_cartesian_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)
   !
   !   call GlobalMesh_SetWellCar(nx, ny, nz)
   !
   !   call GlobalMesh_Make_post_read()
   !
   !   call DefWell_Make_SetDataWell(NbWellInj, NbWellProd)
   !   call DefWell_Make_ComputeWellIndex( &
   !      NbNode, XNode, CellbyNode, NodebyCell, FracbyNode, NodebyFace, &
   !      PermCell, PermFrac)
   !
   !end subroutine NN_init_build_grid

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
      call PetscInitialize(PETSC_NULL_CHARACTER, Ierr); CHKERRQ(Ierr)
   
      ! cf. https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInitializeFortran.html
      call PetscInitializeFortran(Ierr); CHKERRQ(Ierr)

      ! init mpi, communicator/commRank/commSize
      call CommonMPI_init(PETSC_COMM_WORLD)
   
      call NN_init_output_streams(Logfile)
   
      comptime_readmesh = MPI_WTIME()
   
   end subroutine NN_init_warmup

subroutine NN_init_phase2(OutputDir, activate_cpramg, activate_direct_solver)

      character(len=*), intent(in) :: OutputDir
      logical(c_bool), intent(in) :: activate_cpramg, activate_direct_solver
      logical :: ok

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

      ! csr sturcture of Jacobian
      ! allocate memory of Jacobian ans Sm
      call Jacobian_StrucJacBigA

      ! csr sturcture of Jacobian after Schur
      ! allocate memory of Jacobian and Sm after Schur
      call Jacobian_StrucJacA

      call SolvePetsc_Init(KspNiterMax, KspTol, activate_cpramg, activate_direct_solver)

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

      comptime_part = MPI_WTIME() - comptime_start
      comptime_start = MPI_WTIME()
      if (commRank == 0) then
         do i = 1, size(fd)
            write (fd(i), '(A,F16.3)') "Computation time of allocation:   ", &
               comptime_part
         end do
      end if

      comptime_total = comptime_total + comptime_part

      ! init and sort for flash
      call DefFlashWells_allocate

      comptime_total = comptime_total + (MPI_WTIME() - comptime_start)
      comptime_start = MPI_WTIME()

      ! *** Time steps *** !

      TimeCurrent = 0.d0
      TimeOutput = 0.d0

      Delta_t = TimeStepInit

   end subroutine NN_init_phase2

   subroutine NN_flash_control_volumes(n, state, rocktypes, volume)

   integer, intent(in) :: n
   type(Type_IncCVReservoir), intent(inout) :: state(n)
   integer, intent(in) :: rocktypes(IndThermique + 1, n)
   double precision, intent(in) :: volume(n)

   integer :: k

   do k = 1, n
       call DefFlash_Flash_cv(state(k), rocktypes(:, k), volume(k))
   end do

   end subroutine NN_flash_control_volumes

   !> \brief Main surboutine, after each Newton iteration
   !! execute the flash to determine the phases
   !! which are actualy present, and
   !! the mode of the well (flowrate or pressure).
   subroutine NN_flash_all_control_volumes() &
       bind(C, name="NN_flash_all_control_volumes")

   call NN_flash_control_volumes(NbNodeLocal_Ncpus(commRank + 1), IncNode, NodeRocktypeLocal, PoroVolDarcyNode)
   call NN_flash_control_volumes(NbFracLocal_Ncpus(commRank + 1), IncFrac, FracRocktypeLocal, PoroVolDarcyFrac)
   call NN_flash_control_volumes(NbCellLocal_Ncpus(commRank + 1), IncCell, CellRocktypeLocal, PoroVolDarcyCell)

   ! choose between linear or non-linear update of the Newton unknown Pw
   ! The next subroutines also compute the mode of the wells ('pressure' or 'flowrate')
   call DefFlashWells_NewtonFlashLinWells

   end subroutine NN_flash_all_control_volumes

   subroutine NN_main_make_timestep(initial_time_step)

      real(c_double), optional, intent(in) :: initial_time_step
      KSPConvergedReason :: ksp_converged_reason

      if(present(initial_time_step)) Delta_t = max(0.d0, initial_time_step)

      ! init start time
      comptime_start = MPI_WTIME()

      NewtonConv = .false.
      KspConv = .false.

      ! save current status, copy Inc* to Inc*PrevisousTimeStep
      call IncCV_SaveIncPreviousTimeStep

      call IncCVWells_PressureDrop
      
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
            call IncCVWells_PressureDrop

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
            if (NewtonIter==1) call Residu_reset_history
            call Residu_compute(Delta_t)

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
                        write (j, '(A,ES12.5)', advance='no') "  ", NewtonResConvInit%values(k)
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

            if (NewtonIter > 1 .and. NewtonResNormRel(NewtonIter) < NewtonTol) then

               NewtonConv = .true.
               KspConv = .true.
               !write(*,*) 'Newton converged with', NewtonIter, 'iterations (dt=', Delta_t, ')'
               exit
            !else
               !write(*,*) 'Newton did NOT converged after', NewtonIter, 'iterations (dt=', Delta_t, ')'
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
            ksp_converged_reason = SolvePetsc_KspSolve()

            if (ksp_converged_reason < 0) then ! not converge
               write(*,*)
               write(*,*)
               write(*,*) 'Iterative solver did not converge with reason', ksp_converged_reason
               !call SolvePetsc_dump_system('ksp')
               ! call MPI_Abort(ComPASS_COMM_WORLD, errcode, Ierr)
               KspConv = .false. ! solver not converge
               KspNbFailure = KspNbFailure + 1
               ! restart a new Newton with a smaller time step
               Delta_t = Delta_t*0.5d0
               ! load status
               call IncCV_LoadIncPreviousTimeStep
               call IncCVWells_PressureDrop
               ! print residu if Ksp failure
               KspNiter = SolvePetsc_KspSolveIterationNumber()
               call SolvePetsc_Ksp_history(KspHistory)
               if (commRank == 0) then
                  do i = 1, size(fd)
                     j = fd(i)
                     write (j, *) ""
                     do k = 1, KspNiter
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
                KspNiter = SolvePetsc_KspSolveIterationNumber()
                if (commRank == 0) then
                   do i = 1, size(fd)
                      write (fd(i), '(A,I5)', advance="no") "    Nb of Ksp iter:  ", KspNiter
                   end do
                end if

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

               call NN_flash_all_control_volumes

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
            call IncCVWells_PressureDrop

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

     ! ! *** Report *** !
     ! if (commRank == 0) then
     !    do i = 1, size(fd)
     !       j = fd(i)
     !       write (j, *) ""
     !       write (j, *) ""
     !       write (j, *) "Final Report"

     !       write (j, *) ""
     !       write (j, *) "    *Final time:  ", TimeFinal/OneSecond

     !       write (j, *) ""
     !       write (j, *) "    *Mesh:"
     !       write (j, '(A,I0)') "        -Nb of cells:  ", NbCell
     !       write (j, '(A,I0)') "        -Nb of faces:  ", NbFace
     !       write (j, '(A,I0)') "        -Nb of nodes:  ", NbNode
     !       write (j, '(A,I0)') "        -Nb of fracs:  ", NbFrac

     !       write (j, *) ""
     !       write (j, *) "    *Newton/Ksp Iterations:"
     !       write (j, '(A,I0)') "        -Total nb of Newton iters:  ", NewtonNiterTotal
     !       write (j, '(A,I0)') "        -Total nb of Ksp iters:  ", KspNiterTotal

     !       write (j, *) ""
     !       write (j, *) "    *Newton/Ksp Failures:"
     !       write (j, '(A,I0)') "        -Total nb of Newton failures:  ", NewtonNbFailure
     !       write (j, '(A,I0)') "        -Total nb of Ksp failures:  ", KspNbFailure

     !       write (j, *) ""
     !       write (j, *) "    *Total simulation time:  ", comptime_total
     !    end do
     ! end if

      ! *** Free *** !

      deallocate (NewtonIncreNode)
      deallocate (NewtonIncreFrac)
      deallocate (NewtonIncreCell)
      deallocate (NewtonIncreWellInj)
      deallocate (NewtonIncreWellProd)
      if (allocated(KspHistory)) then
        deallocate (KspHistory)
      end if
      if (allocated(fd)) then
         deallocate (fd)
      end if

      call SolvePetsc_free
      call Jacobian_free
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
     !    write (*, *) ""
     !    write (*, *) "NN Finalize"
      end if
      call PetscFinalize(Ierr)

   end subroutine NN_finalize

end module NN
