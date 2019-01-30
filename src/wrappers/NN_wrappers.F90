!
! This file is part of ComPASS.
!
! ComPASS is free software: you can redistribute it and/or modify it under both the terms
! of the GNU General Public License version 3 (https://www.gnu.org/licenses/gpl.html),
! and the CeCILL License Agreement version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
!


    module NNWrapper

       use, intrinsic :: iso_c_binding

       use NN
       use StringWrapper

       implicit none

       public :: &
          !NN_init_from_C, &
          !NN_init_warmup_and_read_mesh_from_C, &
          NN_init_warmup_from_C, &
          !NN_init_read_mesh_from_C, &
          !NN_init_build_grid_from_C, &
          NN_init_phase2_from_C, &
          !NN_main_from_C, &
          ! NN_main_make_timestep_from_C, &
          !NN_main_output_visu_from_C, &
          NN_main_summarize_time_step_from_C, &
          NN_finalize_from_C

    contains

       !subroutine NN_init_from_C(MeshFile, LogFile, OutputDir) &
       !   bind(C, name="NN_init")
       !
       !   type(cpp_string_wrapper), intent(in) :: MeshFile, LogFile, OutputDir
       !
       !   call NN_init(fortran_string(MeshFile), fortran_string(LogFile), fortran_string(OutputDir))
       !
       !end subroutine NN_init_from_C
       !
       !subroutine NN_init_warmup_and_read_mesh_from_C(MeshFile, LogFile) &
       !   bind(C, name="NN_init_warmup_and_read_mesh")
       !
       !   type(cpp_string_wrapper), intent(in) :: MeshFile, LogFile
       !
       !   call NN_init_warmup_and_read_mesh(fortran_string(MeshFile), fortran_string(LogFile))
       !
       !end subroutine NN_init_warmup_and_read_mesh_from_C
       
       subroutine NN_init_phase2_from_C(OutputDir, activate_cpramg, activate_direct_solver) &
          bind(C, name="NN_init_phase2")
       
          type(cpp_string_wrapper), intent(in) :: OutputDir
          logical(c_bool), intent(in), value :: activate_cpramg, activate_direct_solver
       
          call NN_init_phase2(fortran_string(OutputDir), activate_cpramg, activate_direct_solver)
       
       end subroutine NN_init_phase2_from_C
       
       subroutine NN_init_warmup_from_C(Logfile) &
          bind(C, name="NN_init_warmup")
       
          type(cpp_string_wrapper), intent(in) :: Logfile
       
          call NN_init_warmup(fortran_string(Logfile))
       
       end subroutine NN_init_warmup_from_C

       !subroutine NN_init_build_grid_from_C(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz) &
       !   bind(C, name="NN_init_build_grid")
       !
       !   real(kind=c_double), value, intent(in)  :: Ox, Oy, Oz
       !   real(kind=c_double), value, intent(in)  :: lx, ly, lz
       !   integer(kind=c_int), value, intent(in)  :: nx, ny, nz
       !
       !   call NN_init_build_grid(Ox, Oy, Oz, lx, ly, lz, nx, ny, nz)
       !
       !end subroutine NN_init_build_grid_from_C

       !subroutine NN_init_read_mesh_from_C(MeshFile) &
       !   bind(C, name="NN_init_read_mesh")
       !
       !   type(cpp_string_wrapper), intent(in) :: MeshFile
       !
       !   call NN_init_read_mesh(fortran_string(MeshFile))
       !
       !end subroutine NN_init_read_mesh_from_C
       !
       !subroutine NN_main_from_C(TimeIter, OutputDir) &
       !   bind(C, name="NN_main")
       !
       !   integer(c_int), value, intent(in) :: TimeIter
       !   type(cpp_string_wrapper), intent(in) :: OutputDir
       !   integer                              :: TimeIter_tmp
       !
       !   !FIXME: using tmp variable because TimeIter is in/out in NN_main
       !   TimeIter_tmp = TimeIter
       !   call NN_main(TimeIter_tmp, fortran_string(OutputDir))
       !
       !end subroutine NN_main_from_C
       
       !subroutine NN_main_make_timestep_from_C(initial_timestep) &
       !   bind(C, name="NN_main_make_timestep")
       !
       !   real(c_double), value, intent(in) :: initial_timestep
       !
       !   call NN_main_make_timestep(initial_timestep)
       !
       !end subroutine NN_main_make_timestep_from_C

       !subroutine NN_main_output_visu_from_C(TimeIter, OutputDir) &
       !   bind(C, name="NN_main_output_visu")
       !
       !   integer(c_int), value, intent(in) :: TimeIter
       !   type(cpp_string_wrapper), intent(in) :: OutputDir
       !
       !   call NN_main_output_visu(TimeIter, fortran_string(OutputDir))
       !
       !end subroutine NN_main_output_visu_from_C

       subroutine NN_main_summarize_time_step_from_C() &
          bind(C, name="NN_main_summarize_timestep")

          call NN_main_summarize_timestep

       end subroutine NN_main_summarize_time_step_from_C

       subroutine NN_finalize_from_C() &
          bind(C, name="NN_finalize")

          call NN_finalize()

       end subroutine NN_finalize_from_C

    end module NNWrapper
