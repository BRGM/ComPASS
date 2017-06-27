
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
    Residu.F90
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

if(DEF_VISU)
    list(APPEND ComPASSCore_SRCS
        VisuVTK.F90
        VisuVTK.cpp
    )
endif(DEF_VISU)

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
    # default
    # Injector_Qimposed_frac
    # Injector_Pimposed_frac
    # Productor_Pimposed_frac
    # Productor_Qimposed_frac
    monothermal_1inj_1prod 
    andra
)

set( CONF_SRCS
    DefFlash.F90
    DefModel.F90
    DefWell.F90
)

foreach(_conf ${CONFIGURATIONS})
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
    if(DEF_VISU)
        target_include_directories(${CONF_TARGETNAME} PUBLIC ${VTK_INCLUDE_DIRS})
        target_link_libraries(${CONF_TARGETNAME} ${VTK_LIBRARIES})
    endif(DEF_VISU)
endforeach()

#CHECKME: To be removed eventually (dependance on boost::variant -> C++17 variant)
#Boost headers only
find_package(Boost REQUIRED)
if(NOT ${Boost_FOUND})
  message(FATAL_ERROR "Boost libraries are mandatory when building python module, CMake will stop.")
endif(NOT ${Boost_FOUND})

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
  MeshSchema_wrappers.F90
  NN_wrappers.h
  NN_wrappers.cpp
  NN_wrappers.F90
  PathUtilities.F90
  PathUtilities.cpp
  PyBuffer_wrappers.h
  PyBuffer_wrappers.cpp
#  Score-P.cpp
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

## 
## FIXME
## WARNING: Only the first configuration is used to generate the python module.
##          All of this is to be removed when all configuration options are factorized.
##
unset(python_module_done)

foreach(_conf ${CONFIGURATIONS})
    set(CONF_DIR conf/${_conf})
    set(CONF_TARGETNAME ComPASSlib-${_conf})
    set(CONF_ALL_SRCS ${ComPASSCore_SRCS})
    foreach(_src ${CONF_SRCS})
        list(APPEND CONF_ALL_SRCS ${CONF_DIR}/${_src})
    endforeach()
    add_library(${CONF_TARGETNAME} SHARED ${CONF_ALL_SRCS} ${WRAPPERS_SRCS})
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
    target_include_directories(${CONF_TARGETNAME} PUBLIC ${Boost_INCLUDE_DIRS})
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
    if(DEF_HDF5)
        target_include_directories(${CONF_TARGETNAME} PUBLIC ${HDF5_INCLUDE_DIRS})
        target_compile_definitions(${CONF_TARGETNAME} PUBLIC ${HDF5_Fortran_DEFINITIONS})
        target_link_libraries(${CONF_TARGETNAME} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
    endif(DEF_HDF5)
    if(DEF_VISU)
        target_include_directories(${CONF_TARGETNAME} PUBLIC ${VTK_INCLUDE_DIRS})
        target_link_libraries(${CONF_TARGETNAME} ${VTK_LIBRARIES})
    endif(DEF_VISU)

    ## FIXME: Only the first configuration is used to generate the python module.
    if(DEF_Python_BINDINGS)
        if(NOT python_module_done)
          pybind11_add_module(ComPASS wrappers/ComPASS_module.cpp)
        set_target_properties(ComPASS
            PROPERTIES
            CXX_STANDARD 14
            CXX_STANDARD_REQUIRED YES
            CXX_EXTENSIONS NO
        )
	include_directories(${pybind11_INCLUDE_DIR})
	target_include_directories(ComPASS PUBLIC ${pybind11_INCLUDE_DIR})
	  target_link_libraries(ComPASS PRIVATE ${CONF_TARGETNAME})
          set(python_module_done TRUE)
          set(pyComPASS_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/python/ComPASS")
          set_target_properties(ComPASS
              PROPERTIES
              LIBRARY_OUTPUT_DIRECTORY ${pyComPASS_OUTPUT_DIRECTORY}
              LIBRARY_OUTPUT_DIRECTORY_DEBUG ${pyComPASS_OUTPUT_DIRECTORY}
              LIBRARY_OUTPUT_DIRECTORY_RELEASE ${pyComPASS_OUTPUT_DIRECTORY}
          )
        endif(NOT python_module_done)
    endif(DEF_Python_BINDINGS)
    
endforeach()

#CHECKME: Should we install libraries/binaries?
# install(TARGETS NN DESTINATION bin)
    
