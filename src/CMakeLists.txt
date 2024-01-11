# Configuration to be built
set(CONFIGURATIONS brine diphasic immiscible2ph linear_water water2ph)

option(ComPASS_WITH_ALL_PHYSICS "build all physics modules" OFF)
foreach(_conf ${CONFIGURATIONS})
  option(ComPASS_WITH_${_conf}_PHYSICS "build ${_conf} module" OFF)
endforeach()

include(all_sources.cmake)

set(ComPASS_AT_LEAST_ONE_PHYSIC FALSE)
foreach(_conf ${CONFIGURATIONS})
  if(${ComPASS_WITH_${_conf}_PHYSICS})
    set(ComPASS_AT_LEAST_ONE_PHYSIC TRUE)
    break()
  endif()
endforeach()
if(ComPASS_WITH_ALL_PHYSICS)
  set(ComPASS_AT_LEAST_ONE_PHYSIC TRUE)
endif()

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

if(${ComPASS_AT_LEAST_ONE_PHYSIC})
  add_library(ComPASSCommonCore STATIC ${ComPASSCommonCore_SRCS})
  set_target_properties(
    ComPASSCommonCore
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
               Fortran_MODULE_DIRECTORY
               ${PROJECT_BINARY_DIR}/bin/fortran/common
  )
  target_include_directories(
    ComPASSCommonCore
    PUBLIC wrappers
    PUBLIC ${MPI_Fortran_INCLUDE_PATH}
    PUBLIC ${METIS_INCLUDE_DIR}
    PUBLIC ${PETSc_INCLUDE_DIRS}
  )
endif(${ComPASS_AT_LEAST_ONE_PHYSIC})

set(PHYSICS_TARGETS)

foreach(_conf ${CONFIGURATIONS})
  if(NOT ComPASS_WITH_ALL_PHYSICS AND NOT ${ComPASS_WITH_${_conf}_PHYSICS})
    message(NOTICE ">> Skipping ${_conf} physics")
    continue()
  endif()
  set(CONF_DIR conf/${_conf})
  set(CONF_TARGETNAME ComPASSlib-${_conf})
  list(APPEND PHYSICS_TARGETS ${CONF_TARGETNAME})
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
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               Fortran_MODULE_DIRECTORY
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

  set(physics_target ${_conf})
  pybind11_add_module(${physics_target} wrappers/ComPASS_module.cpp)
  target_compile_definitions(
    ${physics_target} PRIVATE ComPASS_CONFIGURATION_NAME=${physics_target}
  )
  target_link_libraries(${physics_target} PRIVATE ${CONF_TARGETNAME})

  install(
    TARGETS ${physics_target}
    DESTINATION ${COMPASS_PYTHON_PACKAGE_NAME}/physics
    EXPORT ComPASS
  )

endforeach()

# add experimental tests to compilation
if(ComPASS_COMPILES_EXPERIMENTAL_PARTS)
  add_subdirectory(experimental)
endif(ComPASS_COMPILES_EXPERIMENTAL_PARTS)

# generate version info FIXME: version_info is regenerated each time (even all
# physics targets are up to date)
add_custom_command(
  OUTPUT ${PROJECT_SOURCE_DIR}/ComPASS/version_info
  COMMAND ${Python_EXECUTABLE} cmake/generate_version_info.py >
          ComPASS/version_info
  DEPENDS ${PHYSICS_TARGETS}
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)
add_custom_target(
  version_info_file ALL SOURCES ${PROJECT_SOURCE_DIR}/ComPASS/version_info
)
