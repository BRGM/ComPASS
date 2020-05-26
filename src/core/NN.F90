!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!

module NN

   use iso_c_binding, only: c_bool, c_int, c_double, c_f_pointer
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
      NodeRocktypeLocal, CellRocktypeLocal, FracRocktypeLocal, &
      NbNodeLocal_Ncpus, NbCellLocal_Ncpus, NbFracLocal_Ncpus, &
      NbWellInjLocal_Ncpus, NbWellProdLocal_Ncpus, &
      MeshSchema_make, MeshSchema_free
   use VAGFrac, only: &
      PoroVolDarcy, VAGFrac_free

   use LocalMesh, only: LocalMesh_Make, LocalMesh_Free
   use NumbyContext, only: NumbyContext_make, NumbyContext_free
   use IncCV, only: IncCV_allocate, IncCV_free
   use DirichletContribution, only: DirichletContribution_allocate, DirichletContribution_free
   use NeumannContribution, only: NeumannContribution_allocate, NeumannContribution_free
   use IncPrimSecd, only: IncPrimSecd_allocate, IncPrimSecd_free
   use Loisthermohydro, only: LoisThermoHydro_allocate, LoisThermoHydro_free
   use Flux, only: Flux_allocate, Flux_free
   use Jacobian, only: Jacobian_StrucJacA, Jacobian_StrucJacBigA, Jacobian_free
   use SolvePetsc, only: SolvePetsc_Init, SolvePetsc_free
   use DefFlash, only: DefFlash_Flash_cv
   use DefFlashWells, only: DefFlashWells_allocate, DefFlashWells_NewtonFlashLinWells

#ifdef COMPASS_PETSC_VERSION_LESS_3_6
#include <finclude/petscdef.h>
#else
#include <petsc/finclude/petsc.h>
#endif

   use petsc

   implicit none

   integer, allocatable, dimension(:) :: fd

   ! ksp variables
   double precision, allocatable, dimension(:) :: KspHistory

   ! Newton increment (NbInc,NbNodelocal)
   real(c_double), dimension(:, :), allocatable :: &
      NewtonIncreNode, &
      NewtonIncreFrac, &
      NewtonIncreCell

   ! as well have only one unknown (pressure), only a vector (NbWellLocal)
   real(c_double), dimension(:), allocatable :: &
      NewtonIncreWellInj, &
      NewtonIncreWellProd

   !  double precision :: f, CC(NbComp), SS(NbPhase), dPf, dTf, dCf(NbComp), dSf(NbPhase)

   public :: &
      NN_wait_for_debug, &
      NN_flash_all_control_volumes, &
      NN_init_warmup, &
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

      integer :: Ierr

      ! initialisation petsc/MPI
      call PetscInitialize(PETSC_NULL_CHARACTER, Ierr); CHKERRQ(Ierr)

      ! cf. https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Sys/PetscInitializeFortran.html
      call PetscInitializeFortran(Ierr); CHKERRQ(Ierr)

      ! init mpi, communicator/commRank/commSize
      call CommonMPI_init(PETSC_COMM_WORLD)

      call NN_init_output_streams(fortran_string(Logfile))

   end subroutine NN_init_warmup

   subroutine NN_wait_for_debug(debugwait) &
      bind(C, name="NN_wait_for_debug")
      logical(c_bool), intent(in), value :: debugwait

      do while (debugwait .EQV. .TRUE.)

      end do

   end subroutine NN_wait_for_debug

   subroutine NN_init_phase2_partition(colors) &
      bind(C, name="NN_init_phase2_partition")
      type(cpp_array_wrapper), intent(in) :: colors
      integer(c_int), pointer :: cell_colors(:)

      if (.NOT. commRank == 0) &
         call CommonMPI_abort("Mesh is supposed to be partitioned by master process.")
      if (colors%n /= NbCell) &
         call CommonMPI_abort("wrong number of cells")

      call c_f_pointer(colors%p, cell_colors, [colors%n])
      call NN_partition_mesh(cell_colors)

   end subroutine NN_init_phase2_partition

   subroutine NN_init_phase2_build_local_mesh() &
      bind(C, name="init_phase2_build_local_mesh")

      call MeshSchema_make

      ! free some tmps in LocalMesh
      if (commRank == 0) then
         call LocalMesh_Free
      end if

   end subroutine NN_init_phase2_build_local_mesh

   subroutine NN_init_phase2_setup_contexts() &
      bind(C, name="init_phase2_setup_contexts")

      call NumbyContext_make(get_model_configuration())

   end subroutine NN_init_phase2_setup_contexts

   subroutine NN_init_phase2_setup_solvers() &
      bind(C, name="init_phase2_setup_solvers")

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

      ! init and sort for flash
      call DefFlashWells_allocate

   end subroutine NN_init_phase2_setup_solvers

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

      call NN_flash_control_volumes(NbNodeLocal_Ncpus(commRank + 1), IncNode, NodeRocktypeLocal, PoroVolDarcy%nodes)
      call NN_flash_control_volumes(NbFracLocal_Ncpus(commRank + 1), IncFrac, FracRocktypeLocal, PoroVolDarcy%fractures)
      call NN_flash_control_volumes(NbCellLocal_Ncpus(commRank + 1), IncCell, CellRocktypeLocal, PoroVolDarcy%cells)

      ! choose between linear or non-linear update of the Newton unknown Pw
      ! The next subroutines also compute the mode of the wells ('pressure' or 'flowrate')
      call DefFlashWells_NewtonFlashLinWells

   end subroutine NN_flash_all_control_volumes

   subroutine NN_finalize() &
      bind(C, name="NN_finalize")

      deallocate (NewtonIncreNode)
      deallocate (NewtonIncreFrac)
      deallocate (NewtonIncreCell)
      deallocate (NewtonIncreWellInj)
      deallocate (NewtonIncreWellProd)
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

   end subroutine NN_finalize

end module NN
