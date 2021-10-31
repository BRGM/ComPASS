# FIXME: Horrible! This following are case dependant configuration files which
# are included in other modules with hard coded include directives
# conf/DefGeometry.F90 -> src/GlobalMesh.F90 conf/DefInitBCvalues.F90 ->
# src/IncCV.F90 Following files are also configuration files but they must be
# compiled ! conf/DefFlash.F90 conf/DefModel.F90 conf/DefWell.F90

# When there is no longer hard-coded include directives this should be enough
# add_subdirectory(core)

# Configuration to be built FIXME: include file liquid_water,
# linear_liquid_water, water_with_tracer have been moved to attic
set(CONFIGURATIONS brine diphasic immiscible2ph linear_water water2ph)

foreach(_conf ${CONFIGURATIONS})
  option(ComPASS_WITH_${_conf}_PHYSICS "Build ${_conf} module" ON)
endforeach()

if(DEF_FREEFLOW_STRUCTURES)
  list(APPEND CONFIGURATIONS diphasic_FreeFlowBC)
endif(DEF_FREEFLOW_STRUCTURES)
cmake_dependent_option(
  ComPASS_WITH_diphasic_FreeFlowBC_PHYSICS "Build diphasic_FreeFlowBC module"
  ON "DEF_FREEFLOW_STRUCTURES" OFF
)

set(ComPASSCommonCore_SRCS
    CommonMPI.F90
    CommonType.F90
    InteroperabilityStructures.F90
    MeshInfo.F90
    NumbyContext.F90
    Physics.F90
    SchemeParameters.F90
)

unset(_tmp)
foreach(_src ${ComPASSCommonCore_SRCS})
  list(APPEND _tmp core/${_src})
endforeach()
set(ComPASSCommonCore_SRCS ${_tmp})

set(ComPASS_NEEDS_CORE_COMPILATION OFF)
foreach(_conf ${CONFIGURATIONS})
  if(${ComPASS_WITH_${_conf}_PHYSICS})
    set(ComPASS_NEEDS_CORE_COMPILATION ON)
    break()
  endif()
endforeach()

if(${ComPASS_NEEDS_CORE_COMPILATION})
  add_library(ComPASSCommonCore STATIC ${ComPASSCommonCore_SRCS})
  set_target_properties(
    ComPASSCommonCore
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
               Fortran_MODULE_DIRECTORY
               ${PROJECT_BINARY_DIR}/bin/fortran/common
  )
  # target_include_directories(ComPASSCommonCore PUBLIC ${CONF_DIR})
  target_include_directories(ComPASSCommonCore PUBLIC wrappers)
  # target_include_directories(ComPASSCommonCore PUBLIC
  # ${CMAKE_Fortran_MODULE_DIRECTORY})
  # target_include_directories(ComPASSCommonCore PUBLIC ${PYTHON_INCLUDE_DIRS})
  target_include_directories(
    ComPASSCommonCore PUBLIC ${MPI_Fortran_INCLUDE_PATH}
  )
  # target_link_libraries(ComPASSCommonCore ${MPI_Fortran_LIBRARIES})
  target_include_directories(ComPASSCommonCore PUBLIC ${METIS_INCLUDE_DIR})
  # target_link_libraries(ComPASSCommonCore ${METIS_LIBRARY})
  target_include_directories(ComPASSCommonCore PUBLIC ${PETSc_INCLUDE_DIRS})
  # FIXME: Is ${LAPACK_INCLUDE_DIR} defined?
  # target_link_libraries(ComPASSCommonCore ${PETSc_LIBRARIES}
  # ${LAPACK_LIBRARIES})
endif(${ComPASS_NEEDS_CORE_COMPILATION})

set(ComPASSCore_SRCS
    DebugUtils.F90
    GlobalMesh.F90
    IncCV.F90
    IncCVReservoir.F90
    IncCVReservoirTypes.F90
    IncCVWells.F90
    DirichletContribution.F90
    NeumannContribution.F90
    Flux.F90
    Jacobian.F90
    LocalMesh.F90
    IncPrimSecd.F90
    IncPrimSecdTypes.F90
    LoisThermoHydro.cpp
    LoisThermoHydro.F90
    MeshSchema.F90
    Newton.F90
    Residu.F90
    Simulation.cpp
    SolvePetsc.F90
    SyncPetsc.F90
    VAGFrac.F90
    NN.F90
    physics/PhysicalConstants.F90
)

set(ComPASSWells_SRCS DefWell.F90 DefFlashWells.F90 WellState.F90)

set(ComPASSMSWells_SRCS DefMSWell.F90)

set(ComPASSFreeFlow_SRCS IncPrimSecdFreeFlow.F90 FreeFlow.F90)

# FIXME: DISPMODULE from Thibault branch ?
if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "MPI is required.")
endif(NOT MPI_Fortran_FOUND)

# FIXME: DISPMODULE from Thibault branch ?
if(DEF_DISPMODULE)
  list(APPEND ComPASSCore_SRCS dispmodule.f90)
endif(DEF_DISPMODULE)

unset(_tmp)
foreach(_src ${ComPASSCore_SRCS})
  list(APPEND _tmp core/${_src})
endforeach()
set(ComPASSCore_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSWells_SRCS})
  list(APPEND _tmp wells/${_src})
endforeach()
set(ComPASSWells_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSMSWells_SRCS})
  list(APPEND _tmp mswells/${_src})
endforeach()
set(ComPASSMSWells_SRCS ${_tmp})

unset(_tmp)
foreach(_src ${ComPASSFreeFlow_SRCS})
  list(APPEND _tmp freeflow/${_src})
endforeach()
set(ComPASSFreeFlow_SRCS ${_tmp})

set(CONF_SRCS DefFlash.F90 DefModel.F90 Equilibriums.cpp Thermodynamics.F90)

# set( WRAPPERS_CORE_SRCS MeshSchema_wrappers.F90 MeshUtilities_wrappers.cpp
# NN_wrappers.cpp NN_wrappers.F90 )

set(WRAPPERS_SRCS
    COC.h
    COC_wrappers.cpp
    CommonTypesWrapper.F90
    DebugUtils_wrappers.cpp
    Flux_wrappers.cpp
    Flux_wrappers.F90
    GlobalMeshWrapper.F90
    GlobalMesh_wrappers.cpp
    GlobalVariables_wrappers.F90
    GlobalVariables_wrappers.cpp
    IncCV_wrappers.cpp
    FreeFlow_wrappers.cpp
    IncCV_wrappers.F90
    LinearSystem_wrapper.cpp
    LinearSystemBuilder.h
    LocalMeshWrapper.F90
    MeshUtilities.h
    MeshSchema_wrappers.cpp
    MeshUtilities_wrappers.cpp
    Metis_wrapper.cpp
    Model_wrappers.h
    Model_common_wrappers.cpp
    MeshSchema_wrappers.F90
    NN_wrappers.cpp
    PyBuffer_wrappers.h
    PyBuffer_wrappers.cpp
    Residu_wrappers.cpp
    SolvePetsc_wrappers.F90
    SolvePetsc_wrappers.cpp
    SyncPetsc_wrappers.cpp
    StringWrapper.h
    StringWrapper.F90
    TimeLoop_wrappers.cpp
    VAGFrac_wrappers.cpp
    Well.h
    Well_wrappers.cpp
    Well_wrappers.F90
)

unset(_tmp)
foreach(_src ${WRAPPERS_SRCS})
  list(APPEND _tmp wrappers/${_src})
endforeach()
set(WRAPPERS_SRCS ${_tmp})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

foreach(_conf ${CONFIGURATIONS})
  if(NOT ${ComPASS_WITH_${_conf}_PHYSICS})
    message(NOTICE ">> Skipping ${_conf} physics")
    continue()
  endif()
  set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/bin/fortran/${_conf})
  set(CONF_DIR conf/${_conf})
  set(CONF_TARGETNAME ComPASSlib-${_conf})
  set(CONF_ALL_SRCS ${ComPASSCore_SRCS} ${ComPASSWells_SRCS}
                    ${ComPASSMSWells_SRCS}
  )
  if(DEF_FREEFLOW_STRUCTURES)
    foreach(_src ${ComPASSFreeFlow_SRCS})
      list(APPEND CONF_ALL_SRCS ${_src})
    endforeach()
  endif(DEF_FREEFLOW_STRUCTURES)
  foreach(_src ${CONF_SRCS})
    list(APPEND CONF_ALL_SRCS ${CONF_DIR}/${_src})
  endforeach()
  add_library(
    ${CONF_TARGETNAME} STATIC
    ${CONF_ALL_SRCS} ${WRAPPERS_SRCS}
    wrappers/conf/${_conf}/Model_specific_wrappers.cpp
  )
  set_target_properties(
    ${CONF_TARGETNAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                  "${CMAKE_BINARY_DIR}/bin"
  )
  target_include_directories(${CONF_TARGETNAME} PUBLIC ${CONF_DIR})
  target_include_directories(${CONF_TARGETNAME} PUBLIC core/include)
  target_include_directories(${CONF_TARGETNAME} PUBLIC wrappers)
  target_include_directories(
    ${CONF_TARGETNAME} PUBLIC ${PROJECT_BINARY_DIR}/bin/fortran/common
  )
  target_include_directories(
    ${CONF_TARGETNAME} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}
  )
  target_include_directories(
    ${CONF_TARGETNAME} PUBLIC ${MPI_Fortran_INCLUDE_PATH}
  )
  set_target_properties(
    ${CONF_TARGETNAME} PROPERTIES COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  )
  target_link_libraries(${CONF_TARGETNAME} pybind11::pybind11)
  target_link_libraries(${CONF_TARGETNAME} ComPASSCommonCore)
  target_link_libraries(${CONF_TARGETNAME} ${MPI_Fortran_LIBRARIES})
  target_include_directories(${CONF_TARGETNAME} PUBLIC ${METIS_INCLUDE_DIR})
  target_link_libraries(${CONF_TARGETNAME} ${METIS_LIBRARY})
  target_include_directories(${CONF_TARGETNAME} PUBLIC ${PETSc_INCLUDE_DIRS})
  # FIXME: Is ${LAPACK_INCLUDE_DIR} defined?
  target_link_libraries(
    ${CONF_TARGETNAME} ${PETSc_LIBRARIES} ${LAPACK_LIBRARIES}
  )
  include(${CONF_DIR}/cmake.conf)

  set(eos_target ${_conf})
  pybind11_add_module(${eos_target} wrappers/ComPASS_module.cpp)
  target_compile_definitions(
    ${eos_target} PRIVATE ComPASS_CONFIGURATION_NAME=${eos_target}
  )
  target_link_libraries(${eos_target} PRIVATE ${CONF_TARGETNAME})

  install(
    TARGETS ${eos_target}
    DESTINATION ${COMPASS_PYTHON_PACKAGE_NAME}/eos
    EXPORT ComPASS
  )

endforeach()

# add experimental tests to compilation
if(ComPASS_COMPILES_EXPERIMENTAL_PARTS)
  add_subdirectory(experimental)
endif(ComPASS_COMPILES_EXPERIMENTAL_PARTS)
