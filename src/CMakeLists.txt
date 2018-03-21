
# FIXME: Horrible!
# This following are case dependant configuration files which are included
# in other modules with hard coded include directives
#    conf/DefGeometry.F90 -> src/GlobalMesh.F90
#    conf/DefInitBCvalues.F90 -> src/IncCV.F90
# Following files are also configuration files but they must be compiled !
#    conf/DefFlash.F90
#    conf/DefModel.F90
#    conf/DefWell.F90

# When there is no longer hard-coded include directives this should be enough
# add_subdirectory(core)


# SCOREP is used or not used for all sources... to be refined later
# scorep scripts (scorep-mpicc, scorep-mpicxx and scorep-mpif90) shall be on your path
# scorep must also have been built with the --enable-shared option
if(DEF_SCOREP_INSTRUMENTATION)
    find_program(SCOREP_MPICC_SCRIPT scorep-mpicc)
    find_program(SCOREP_MPICXX_SCRIPT scorep-mpicxx)
    find_program(SCOREP_MPIF90_SCRIPT scorep-mpif90)
    set(CMAKE_C_COMPILER ${SCOREP_MPICC_SCRIPT})
    set(CMAKE_CXX_COMPILER ${SCOREP_MPICXX_SCRIPT})
    set(CMAKE_Fortran_COMPILER ${SCOREP_MPIF90_SCRIPT})
    # This will build a instrumented wrapper around the python interpreter
    add_subdirectory(perfanalysis)
endif(DEF_SCOREP_INSTRUMENTATION)


set( ComPASSCore_SRCS
    CommonMPI.F90
    CommonType.F90
    GlobalMesh.F90
    IncCV.F90
    Flux.F90
    Jacobian.F90 
    LocalMesh.F90
    LoisThermoHydro.F90
    MeshSchema.F90
    NumbyContext.F90
    PartitionMesh.F90
    Partition.c
    Physics.F90
    Residu.F90
    SchemeParameters.F90
    SolvePetsc.F90
    VAGFrac.F90
    NN.F90
)

# FIXME: DISPMODULE from Thibault branch ?
if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "MPI is required.")
endif(NOT MPI_Fortran_FOUND)

# FIXME: DISPMODULE from Thibault branch ?
if(DEF_DISPMODULE)
    list(APPEND ComPASSCore_SRCS
        dispmodule.f90
    )
endif(DEF_DISPMODULE)

unset(_tmp)
foreach(_src ${ComPASSCore_SRCS})
    list(APPEND _tmp core/${_src})
endforeach()
set(ComPASSCore_SRCS ${_tmp})

set( MAIN_SRCS
  ComPASS.F90
  PathUtilities.F90
)

unset(_tmp)
foreach(_src ${MAIN_SRCS})
    list(APPEND _tmp main/${_src})
endforeach()
set(MAIN_SRCS ${_tmp})

# Configuration to be built
set( CONFIGURATIONS
    water2ph
    diphasic
)

set( CONF_SRCS
    DefFlash.F90
    DefFlashWells.F90
    DefModel.F90
    DefWell.F90
    Thermodynamics.F90
)

if(DEF_BUILD_COMPASS_v3_FORTRAN_EXECUTABLE)
foreach(_conf ${CONFIGURATIONS})
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/bin/fortran/${_conf})
    set(CONF_DIR conf/${_conf})
    set(CONF_TARGETNAME ComPASS-${_conf})
    set(CONF_ALL_SRCS ${ComPASSCore_SRCS})
    foreach(_src ${CONF_SRCS})
        list(APPEND CONF_ALL_SRCS ${CONF_DIR}/${_src})
    endforeach()
    add_executable(${CONF_TARGETNAME} ${CONF_ALL_SRCS} ${MAIN_SRCS})
    set_target_properties(${CONF_TARGETNAME}
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    
    set_target_properties(${CONF_TARGETNAME} PROPERTIES LINKER_LANGUAGE Fortran)

    target_include_directories(${CONF_TARGETNAME} PUBLIC ${CONF_DIR})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    set_target_properties(${CONF_TARGETNAME}
        PROPERTIES
        COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
    )
    target_link_libraries(${CONF_TARGETNAME} ${MPI_Fortran_LIBRARIES})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${METIS_INCLUDE_DIR})
    target_link_libraries(${CONF_TARGETNAME} ${METIS_LIBRARY})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${PETSC_INCLUDES})
    # FIXME: Is ${LAPACK_INCLUDE_DIR} defined?
    target_link_libraries(${CONF_TARGETNAME} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
    if(DEF_HDF5)
        target_include_directories(${CONF_TARGETNAME} PUBLIC ${HDF5_INCLUDE_DIRS})
        target_compile_definitions(${CONF_TARGETNAME} PUBLIC ${HDF5_Fortran_DEFINITIONS})
        target_link_libraries(${CONF_TARGETNAME} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
    endif(DEF_HDF5)
    include(${CONF_DIR}/cmake.conf)
endforeach()
endif(DEF_BUILD_COMPASS_v3_FORTRAN_EXECUTABLE)

set( WRAPPERS_SRCS
  ArrayWrapper.h
  COC.h
  COC_wrappers.h
  COC_wrappers.cpp
  CommonTypesWrapper.F90
  GlobalMeshWrapper.F90
  GlobalMesh_wrappers.h
  GlobalMesh_wrappers.cpp
  GlobalVariables_wrappers.F90
  GlobalVariables_wrappers.h
  GlobalVariables_wrappers.cpp
  IncCV_wrappers.h
  IncCV_wrappers.cpp
  IncCV_wrappers.F90
  LocalMeshWrapper.F90
  MeshUtilities.h
  MeshUtilities_wrappers.h
  MeshUtilities_wrappers.cpp
  Model_wrappers.h
  Model_wrappers.cpp
  MeshSchema_wrappers.F90
  NN_wrappers.h
  NN_wrappers.cpp
  NN_wrappers.F90
  PathUtilities.F90
  PathUtilities.cpp
  PyBuffer_wrappers.h
  PyBuffer_wrappers.cpp
  StringWrapper.h
  StringWrapper.F90
  Well.h
  Well_wrappers.h
  Well_wrappers.cpp
  Well_wrappers.F90
)

unset(_tmp)
foreach(_src ${WRAPPERS_SRCS})
    list(APPEND _tmp wrappers/${_src})
endforeach()
set(WRAPPERS_SRCS ${_tmp})

set(pyComPASS_DIRECTORY "${PROJECT_SOURCE_DIR}/python/ComPASS")
set(pyComPASS_EOS_DIRECTORY "${pyComPASS_DIRECTORY}/eos")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

foreach(_conf ${CONFIGURATIONS})
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/bin/fortran/${_conf})
    set(CONF_DIR conf/${_conf})
    set(CONF_TARGETNAME ComPASSlib-${_conf})
    set(CONF_ALL_SRCS ${ComPASSCore_SRCS})
    foreach(_src ${CONF_SRCS})
        list(APPEND CONF_ALL_SRCS ${CONF_DIR}/${_src})
    endforeach()
    add_library(${CONF_TARGETNAME} STATIC ${CONF_ALL_SRCS} ${WRAPPERS_SRCS})
    set_target_properties(${CONF_TARGETNAME}
        PROPERTIES
    CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${CONF_DIR})
    target_include_directories(${CONF_TARGETNAME} PUBLIC wrappers)
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${PYTHON_INCLUDE_DIRS})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH})
    set_target_properties(${CONF_TARGETNAME}
        PROPERTIES
        COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
    )
    target_link_libraries(${CONF_TARGETNAME} ${PYTHON_LIBRARIES})
    target_link_libraries(${CONF_TARGETNAME} ${MPI_Fortran_LIBRARIES})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${METIS_INCLUDE_DIR})
    target_link_libraries(${CONF_TARGETNAME} ${METIS_LIBRARY})
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${PETSC_INCLUDES})
    # FIXME: Is ${LAPACK_INCLUDE_DIR} defined?
    target_link_libraries(${CONF_TARGETNAME} ${PETSC_LIBRARIES} ${LAPACK_LIBRARIES})
    include(${CONF_DIR}/cmake.conf)

    if(DEF_Python_BINDINGS)
        set(eos_target ${_conf})
        pybind11_add_module(${eos_target} wrappers/ComPASS_module.cpp)
        target_compile_definitions(${eos_target} PRIVATE ComPASS_CONFIGURATION_NAME=${eos_target})
        set_target_properties(${eos_target}
            PROPERTIES
            CXX_STANDARD 14
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
            LIBRARY_OUTPUT_DIRECTORY ${pyComPASS_EOS_DIRECTORY}
            LIBRARY_OUTPUT_DIRECTORY_DEBUG ${pyComPASS_EOS_DIRECTORY}
            LIBRARY_OUTPUT_DIRECTORY_RELEASE ${pyComPASS_EOS_DIRECTORY}
        )
        include_directories(${pybind11_INCLUDE_DIR})
        target_include_directories(${eos_target} PUBLIC ${pybind11_INCLUDE_DIR})
        target_link_libraries(${eos_target} PRIVATE ${CONF_TARGETNAME})
    endif(DEF_Python_BINDINGS)
    
endforeach()

#CHECKME: Should we install libraries/binaries?
# install(TARGETS NN DESTINATION bin)
    
