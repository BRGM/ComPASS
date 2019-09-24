# As configuration files are hard-coded as includes in some .F90 files
# this CMakeLists.txt configuration does not work but shall work in a future version

add_library( ComPASSCore
    CommonMPI.F90
    CommonType.F90
    ComPASS.F90
    GlobalMesh.F90
    IncCV.F90
    IncCVWells.F90
    Flux.F90
    Jacobian.F90 
    LocalMesh.F90
    IncPrimSecd.F90
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
if(DEF_DISPMODULE)
    target_sources( ComPASSCore PUBLIC
        dispmodule.f90
    )
endif(DEF_DISPMODULE)

if(DEF_VISU)
    target_sources( ComPASSCore PUBLIC
        VisuVTK.F90
        VisuVTK.cpp
    )
endif(DEF_VISU)

target_include_directories(ComPASSCore PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(ComPASSCore PUBLIC ${METIS_INCLUDE_DIR})
target_link_libraries(ComPASSCore ${METIS_LIBRARY})
target_include_directories(ComPASSCore PUBLIC ${PETSc_INCLUDE_DIRS})
# FIXME: Is ${LAPACK_INCLUDE_DIR} defined?
target_link_libraries(ComPASSCore ${PETSc_LIBRARIES} ${LAPACK_lapack_LIBRARY})

if(DEF_HDF5)
    target_include_directories(ComPASSCore PUBLIC ${HDF5_INCLUDE_DIRS})
    target_compile_definitions(ComPASSCore PUBLIC ${HDF5_DEFINITIONS})
    message(STATUS "hdf5 definitions: " ${HDF5_DEFINITIONS})
    message(STATUS "hdf5 includes: "  ${HDF5_INCLUDE_DIRS})
    message(STATUS "hdf5 libraries: " ${HDF5_Fortran_HL_LIBRARIES})
    target_link_libraries(ComPASSCore ${HDF5_Fortran_HL_LIBRARIES})
endif(DEF_HDF5)

if(DEF_VISU)
    target_include_directories(ComPASSCore PUBLIC ${VTK_INCLUDE_DIRS})
    target_link_libraries(ComPASSCore ${VTK_LIBRARIES})
endif(DEF_VISU)
    
#CHECKME: Should we install libraries/binaries?
# install(TARGETS NN DESTINATION bin)
