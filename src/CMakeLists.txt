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

option(ComPASS_WITH_ALL_PHYSICS
       "Build all physics modules - this will overide other physics options"
       OFF
)
foreach(_conf ${CONFIGURATIONS})
  option(ComPASS_WITH_${_conf}_PHYSICS "Build ${_conf} module" OFF)
endforeach()

include(all_sources.cmake)

set(ComPASS_NEEDS_CORE_COMPILATION OFF)
foreach(_conf ${CONFIGURATIONS})
  if(${ComPASS_WITH_${_conf}_PHYSICS})
    set(ComPASS_NEEDS_CORE_COMPILATION ON)
    break()
  endif()
endforeach()
if(ComPASS_WITH_ALL_PHYSICS)
  set(ComPASS_NEEDS_CORE_COMPILATION ON)
endif()

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

# FIXME: DISPMODULE from Thibault branch ?
if(NOT MPI_Fortran_FOUND)
  message(FATAL_ERROR "MPI is required.")
endif(NOT MPI_Fortran_FOUND)

# FIXME: DISPMODULE from Thibault branch ?
if(DEF_DISPMODULE)
  list(APPEND ComPASSCore_SRCS dispmodule.f90)
endif(DEF_DISPMODULE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ComPASS_AT_LEAST_ONE_PHYSIC FALSE)

foreach(_conf ${CONFIGURATIONS})
  if(NOT ComPASS_WITH_ALL_PHYSICS AND NOT ${ComPASS_WITH_${_conf}_PHYSICS})
    message(NOTICE ">> Skipping ${_conf} physics")
    continue()
  endif()
  set(ComPASS_AT_LEAST_ONE_PHYSIC TRUE)
  set(CONF_DIR conf/${_conf})
  set(CONF_TARGETNAME ComPASSlib-${_conf})
  set(CONF_ALL_SRCS ${ComPASSCore_SRCS} ${ComPASSWells_SRCS}
                    ${ComPASSMSWells_SRCS}
  )
  if(${_conf} STREQUAL diphasic)
    foreach(_src ${ComPASSFreeFlow_SRCS})
      list(APPEND CONF_ALL_SRCS ${_src})
    endforeach()
  endif()
  foreach(_src ${CONF_SRCS})
    list(APPEND CONF_ALL_SRCS ${CONF_DIR}/${_src})
  endforeach()
  add_library(
    ${CONF_TARGETNAME} STATIC
    ${CONF_ALL_SRCS} ${WRAPPERS_SRCS}
    wrappers/conf/${_conf}/Model_specific_wrappers.cpp
  )
  target_compile_definitions(
    ${CONF_TARGETNAME} PUBLIC -DComPASS_WITH_${_conf}_PHYSICS
  )
  if(${_conf} STREQUAL diphasic)
    target_compile_definitions(
      ${CONF_TARGETNAME} PUBLIC -D_WITH_FREEFLOW_STRUCTURES_
    )
  endif()
  set_target_properties(
    ${CONF_TARGETNAME}
    PROPERTIES Fortran_MODULE_DIRECTORY
               ${PROJECT_BINARY_DIR}/bin/fortran/${_conf}
               RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
               COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  )
  target_include_directories(
    ${CONF_TARGETNAME}
    PUBLIC ${CONF_DIR}
    PUBLIC core/include
    PUBLIC wrappers
    PUBLIC ${PROJECT_BINARY_DIR}/bin/fortran/common
    PUBLIC ${PROJECT_BINARY_DIR}/bin/fortran/${_conf}
    PUBLIC ${MPI_Fortran_INCLUDE_PATH}
    PUBLIC ${METIS_INCLUDE_DIR}
    PUBLIC ${PETSc_INCLUDE_DIRS}
  )
  target_link_libraries(${CONF_TARGETNAME} pybind11::pybind11)
  target_link_libraries(${CONF_TARGETNAME} ComPASSCommonCore)
  target_link_libraries(${CONF_TARGETNAME} ${MPI_Fortran_LIBRARIES})
  target_link_libraries(${CONF_TARGETNAME} ${METIS_LIBRARY})
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

if(NOT ComPASS_AT_LEAST_ONE_PHYSIC)
  message(NOTICE "")
  message(NOTICE "You must compile at least one physics among:")
  foreach(_conf ${CONFIGURATIONS})
    message(NOTICE "  - " ${_conf})
  endforeach()
  message(NOTICE "")
  message(
    FATAL_ERROR
      "Activate at least one physics with -DComPASS_WITH_name_PHYSICS=ON"
  )
endif()

# add experimental tests to compilation
if(ComPASS_COMPILES_EXPERIMENTAL_PARTS)
  add_subdirectory(experimental)
endif(ComPASS_COMPILES_EXPERIMENTAL_PARTS)
