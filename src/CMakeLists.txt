
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
    target_link_libraries(${CONF_TARGETNAME} ${PETSC_LIBRARIES} ${LAPACK_lapack_LIBRARY})
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

#CHECKME: Should we install libraries/binaries?
# install(TARGETS NN DESTINATION bin)
    
