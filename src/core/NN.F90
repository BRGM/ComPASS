!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module NN

   use iso_c_binding, only: c_bool, c_int, c_double
   use StringWrapper, only: cpp_string_wrapper, fortran_string

   use mpi, only: MPI_Abort, MPI_Barrier, MPI_WTIME
   use CommonMPI, only: &
     commSize, commRank, Ncpus, ComPASS_COMM_WORLD, CommonMPI_init, CommonMPI_abort
   use InteroperabilityStructures, only: cpp_array_wrapper
   use SchemeParameters, only: &
      NewtonNiterMax, KspNiterMax, KspTol, OneSecond, TimeFinal, TimeStepInit

   use GlobalMesh, only: &
      GlobalMesh_free, CellbyCell, &
      NbDirNodeP, NbDirNodeT, NbNode, NbCell, NbFace, NbFrac, NbWellInj, NbWellProd

   use CommonType, only: ModelConfiguration
   use DefModel, only: &
      NbComp, NbPhase, IndThermique, NbIncTotalMax, get_model_configuration

   use IncCVReservoir, only: &
      Type_IncCVReservoir, IncNode, IncCell, IncFrac
   use MeshSchema, only: &
      PermCellLocal, PermFracLocal, CondThermalCellLocal, CondThermalFracLocal, &
      NodeRocktypeLocal, CellRocktypeLocal, FracRocktypeLocal, &
      NbNodeLocal_Ncpus, NbCellLocal_Ncpus, NbFracLocal_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus, &
      MeshSchema_make, MeshSchema_free
   use VAGFrac, only: &
      PoroVolDarcyNode, PoroVolDarcyCell, PoroVolDarcyFrac, &
      VAGFrac_TransDarcy, VAGFrac_TransFourier, &
      VAGFrac_VolsDarcy, VAGFrac_VolsFourier, VAGFrac_free

    use LocalMesh, only: LocalMesh_Make, LocalMesh_Free
    use NumbyContext, only: NumbyContext_make, NumbyContext_free
    use IncCV, only: IncCV_allocate, IncCV_free
    use DirichletContribution, only: DirichletContribution_allocate, DirichletContribution_free
    use NeumannContribution, only: NeumannContribution_allocate, NeumannContribution_free
    use IncPrimSecd, only: IncPrimSecd_allocate, IncPrimSecd_free
    use Loisthermohydro, only: LoisThermoHydro_allocate, LoisThermoHydro_free
    use Flux, only: Flux_allocate, Flux_free
    use Jacobian, only: Jacobian_StrucJacA, Jacobian_StrucJacBigA, Jacobian_free
    use SolvePetsc, only: SolvePetsc_Init, SolvePetsc_free, SolvePetsc_SyncMat
    use DefFlash, only: DefFlash_Flash_cv
    use DefFlashWells, only: DefFlashWells_allocate, DefFlashWells_NewtonFlashLinWells

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
      NN_main_summarize_timestep, &
      NN_flash_all_control_volumes, &
      NN_init_warmup, &
      NN_init_phase2, &
      NN_finalize
   
   private :: &
       NN_flash_control_volumes, &
       NN_init_output_streams, &
       NN_partition_mesh

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

   subroutine NN_partition_mesh(ProcbyCell)
      integer(c_int), dimension(:), intent(in) :: ProcbyCell
    
      if (.NOT. commRank == 0) &
         call CommonMPI_abort("Mesh is supposed to be partitioned by master process.")

      call LocalMesh_Make(ProcbyCell)
      call GlobalMesh_free

   end subroutine NN_partition_mesh

   subroutine NN_init_warmup(LogFile) &
      bind(C, name="NN_init_warmup")
   
      type(cpp_string_wrapper), intent(in) :: Logfile
   
      ! initialisation petsc/MPI
      call PetscInitialize(PETSC_NULL_CHARACTER, Ierr); CHKERRQ(Ierr)
   
      ! cf. https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInitializeFortran.html
      call PetscInitializeFortran(Ierr); CHKERRQ(Ierr)

      ! init mpi, communicator/commRank/commSize
      call CommonMPI_init(PETSC_COMM_WORLD)
   
      call NN_init_output_streams(fortran_string(Logfile))
   
      comptime_readmesh = MPI_WTIME()
   
   end subroutine NN_init_warmup

   subroutine NN_init_phase2_summary() &
      bind(C, name="NN_init_phase2_summary")

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

         ! comptime_readmesh = MPI_WTIME() - comptime_readmesh
         ! do i = 1, size(fd)
            ! write (fd(i), '(A,F16.3)') "Computation time warm up and reading mesh: ", &
               ! comptime_readmesh
         ! end do

      end if

   end subroutine NN_init_phase2_summary

   subroutine NN_init_phase2_partition(colors) &
      bind(C, name="NN_init_phase2_partition")
      type(cpp_array_wrapper), intent(in) :: colors
      integer(c_int), pointer :: cell_colors(:)

      if (.NOT. commRank == 0) &
         call CommonMPI_abort("Mesh is supposed to be partitioned by master process.")
     if(colors%n/=NbCell) &
        call CommonMPI_abort("wrong number of cells")
        
     call c_f_pointer(colors%p, cell_colors, [colors%n])
     call NN_partition_mesh(cell_colors)
          
    end subroutine NN_init_phase2_partition

   subroutine NN_init_phase2(activate_cpramg, activate_direct_solver) &
      bind(C, name="NN_init_phase2")
       
      logical(c_bool), intent(in), value :: activate_cpramg, activate_direct_solver

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

      ! *** Numeratation derived from model *** !

      call NumbyContext_make(get_model_configuration())

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

      ! allocate IncPrimSecd
      call IncPrimSecd_allocate

      ! allocate Loisthermohydro
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
                (NbIncTotalMax, NbNodeLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreFrac &
                (NbIncTotalMax, NbFracLocal_Ncpus(commRank + 1)))
      allocate (NewtonIncreCell &
                (NbIncTotalMax, NbCellLocal_Ncpus(commRank + 1)))
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

   subroutine NN_main_summarize_timestep() &
      bind(C, name="NN_main_summarize_timestep")

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

   subroutine NN_finalize() &
      bind(C, name="NN_finalize")

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
      call IncPrimSecd_free
      call LoisThermoHydro_free
      call IncCV_free
      call DirichletContribution_free
      call NeumannContribution_free
      ! call DefFlash_free
      call NumbyContext_free
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
