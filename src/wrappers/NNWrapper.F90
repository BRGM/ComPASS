
    module NNWrapper

      use, intrinsic :: iso_c_binding

      use NN
      use StringWrapper

      implicit none

      public :: &
        NN_init_from_C, &
      NN_main_from_C, &
      NN_finalize_from_C

    contains

      subroutine NN_init_from_C(MeshFile, LogFile, OutputDir) bind(C, name="NN_init")

        type(cpp_string_wrapper), intent(in) :: MeshFile, LogFile, OutputDir
        
        call NN_init(fortran_string(MeshFile), fortran_string(LogFile), fortran_string(OutputDir))

      end subroutine NN_init_from_C

      subroutine NN_main_from_C(TimeIter, OutputDir) bind(C, name="NN_main")
      
        integer(c_int), value,    intent(in) :: TimeIter
        type(cpp_string_wrapper), intent(in) :: OutputDir
        integer                              :: TimeIter_tmp
        
        !FIXME: using tmp variable because TimeIter is in/out in NN_main
        TimeIter_tmp = TimeIter
        call NN_main(TimeIter_tmp, fortran_string(OutputDir))
      
      end subroutine NN_main_from_C

      subroutine NN_finalize_from_C() bind(C, name="NN_finalize")

        call NN_finalize()

      end subroutine NN_finalize_from_C

    end module NNWrapper
