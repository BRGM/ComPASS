cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)

project(ComPASS Fortran C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE
  )
endif(NOT CMAKE_BUILD_TYPE)

# make sure our own modules are found
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/modules"
)

include(ComPASSCMakeUtils)
include(CMakeDependentOption)

message(STATUS "%% PETSc:" ${PETSC_DIR})

if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
  message(
    FATAL_ERROR
      "-------------------------------
In-source compilation forbidden.
Build a directory build: mkdir build.
Then: cd build; cmake ..
--------------------------------"
  )
endif(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)

# Cf. cmake documentation on RPATH and the related FAQ
# https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling#common-questions
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(COMPASS_CONFIG_DIRECTORY
    ${PROJECT_SOURCE_DIR}/config
    CACHE PATH "path to configuration scripts"
)

# compiled python modules will be put in
# ${CMAKE_INSTALL_PREFIX}/${COMPASS_PYTHON_PACKAGE_NAME} FIXME: could be
# retrieve from setup.py
set(COMPASS_PYTHON_PACKAGE_NAME ComPASS)

# will not override a user-provided installation prefix (e.g. using cmake
# -DCMAKE_INSTALL_PREFIX=...)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(COMPASS_PYTHON_PACKAGE_DIRECTORY)
    set(CMAKE_INSTALL_PREFIX
        ${COMPASS_PYTHON_PACKAGE_DIRECTORY}
        CACHE
          PATH
          "Default installation prefix, will be through COMPASS_PYTHON_PACKAGE_DIRECTORY variable only if no installation prefix is given."
          FORCE
    )
  else()
    # if no specific instalation prefix was provided and if the
    # COMPASS_PYTHON_PACKAGE_DIRECTORY variable is not set default to the source
    # directory cmake export will be put there and python modules will be put in
    # ${CMAKE_INSTALL_PREFIX}/${COMPASS_PYTHON_PACKAGE_NAME}
    set(CMAKE_INSTALL_PREFIX
        ${CMAKE_CURRENT_SOURCE_DIR}
        CACHE
          PATH
          "Default installation prefix, will be through COMPASS_PYTHON_PACKAGE_DIRECTORY variable only if no installation prefix is given."
          FORCE
    )
  endif()
endif()

set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/thirdparty)

# The following variables customize behavior of FetchContent They are all
# intended for the developer and are not set by the project Just uncomment one
# of the following lines according to your need.
# https://cmake.org/cmake/help/latest/module/FetchContent.html When this option
# is enabled, no attempt is made to download or update any content. When no
# changes have been made to any content details, turning this option ON can
# significantly speed up the configure stage.

# set(FETCHCONTENT_FULLY_DISCONNECTED ON CACHE BOOL "no attempt is made to
# download or update any content" FORCE )

# Alternative to the previous option: if content has not been downloaded
# previously, it will still be downloaded when this option is enabled

# set(FETCHCONTENT_UPDATES_DISCONNECTED ON CACHE BOOL "will download any content
# only once" FORCE )

find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
find_package(pybind11 2.6 QUIET) # try locally
if(NOT pybind11_FOUND) # else fetch it from github
  message(STATUS "Could NOT find pybind11 (fetching it from GitHub)")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.6.2
  )
  if(${CMAKE_VERSION} VERSION_LESS 3.14)
    FetchContent_GetProperties(pybind11)
    if(NOT pybind11_POPULATED)
      FetchContent_Populate(pybind11)
      add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
    endif()
  else()
    FetchContent_MakeAvailable(pybind11)
  endif()
endif()

set(PYBIND11_PYTHON_VERSION
    3
    CACHE STRING "Python version to use for compiling modules"
)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-g -fcheck=all -Wall -Wextra -fbacktrace"
      CACHE STRING "Useful gfortran debug compilation options" FORCE
  )
endif()

# ########## Configuration options

option(DEF_THERMIQUE "Thermique" ON)
if(DEF_THERMIQUE)
  add_definitions(-D_THERMIQUE_)
endif(DEF_THERMIQUE)

option(DEF_DISPMODULE "Third Party debug using displaymodule" OFF)
if(DEF_DISPMODULE)
  add_definitions(-D_DISPMODULE_)
endif(DEF_DISPMODULE)

option(DEF_FREEFLOW_STRUCTURES "Use FreeFlow Structures" OFF)
if(DEF_FREEFLOW_STRUCTURES)
  add_definitions(-D_WIP_FREEFLOW_STRUCTURES_)
endif(DEF_FREEFLOW_STRUCTURES)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-D_DEBUG_)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

option(COMPASS_USE_GDB_FORTRAN "Use gdb debugger" OFF)
if(COMPASS_USE_GDB_FORTRAN)
  set(CMAKE_Fortran_FLAGS_DEBUG
      "-g -fcheck=all -Wall -Wextra -fbacktrace -ggdb3 -fno-eliminate-unused-debug-symbols -O0"
      CACHE STRING "Useful gfortran debug compilation options with gdb." FORCE
  )
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-D_COMPASS_FORTRAN_DO_NOT_USE_ONLY_)
  endif(CMAKE_BUILD_TYPE MATCHES Debug)
endif(COMPASS_USE_GDB_FORTRAN)

get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS ${DirDefs}")

option(ComPASS_COMPILES_EXPERIMENTAL_PARTS "Compiles experimental parts" OFF)

# -- scorep intstrumentation ---------------------------------------------------

# If scorep is used to profile code we need to build a wrapper around the python
# interpreter
set(ComPASS_HAS_SCOREP_INSTRUMENTATION FALSE)
foreach(_lang C CXX Fortran)
  get_filename_component(_lang_compiler ${CMAKE_${_lang}_COMPILER} NAME)
  if(${_lang_compiler} MATCHES "^scorep-.*")
    set(ComPASS_HAS_SCOREP_INSTRUMENTATION TRUE)
  endif()
endforeach()
if(${ComPASS_HAS_SCOREP_INSTRUMENTATION})
  add_subdirectory(miscellaneous/scorep)
endif()

# ##############################################################################

# #CHECKME ???? -> portability + test constants are compiled as double
# set(CMAKE_Fortran_FLAGS "-I${LIB_PREFIX}/include/shared -fdefault-real-8
# -fdefault-double-8" CACHE STRING "Fortran FLAGS") # make sure that the
# constants (ex. 1.0) are compiled as double (ex. 1.d0) -fdefault-real-8
# -fdefault-double-8 # add_definitions(-fdefault-real-8 -fdefault-double-8)

# #CHECKME: is this the place ? portability -> Windows ? add_definitions(-O0 -g)

# CHECKME: CMake -> cmake ?
set(CMAKE_FIND_DIR ${CMAKE_SOURCE_DIR}/CMake)

find_package(MPI REQUIRED)
find_package(METIS REQUIRED)
find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  # Try to use pkg-config
  pkg_check_modules(PETSc PETSc)
endif(PKG_CONFIG_FOUND)
if(NOT PETSc_FOUND OR NOT PKG_CONFIG_FOUND)
  # else defaults to local script
  find_package(PETSc REQUIRED)
endif(NOT PETSc_FOUND OR NOT PKG_CONFIG_FOUND)

if(NOT PETSc_FOUND)
  message(
    FATAL_ERROR
      "PETSc could not be found"
      "(use a PETSc version lower or equal to 3.5.4 or greater or equal to 3.8.0)"
  )
endif(NOT PETSc_FOUND)

if(${PETSc_VERSION} VERSION_LESS 3.6)
  add_definitions(-DCOMPASS_PETSC_VERSION_LESS_3_6)
elseif(${PETSc_VERSION} VERSION_LESS 3.8)
  message(
    FATAL_ERROR
      "Your PETSc version is "
      ${PETSc_VERSION}
      ".
    Please use a PETSc version lower or equal to 3.5.4 or greater or equal to 3.8.0."
  )
endif(${PETSc_VERSION} VERSION_LESS 3.6)

find_package(LAPACK REQUIRED)

enable_testing()

# CHECKME: Should languages be enabled here or at the beginning as project
# properties enable_language(Fortran) enable_language(C) enable_language(CXX)

include_directories(${THIRD_PARTY_DIR}/mapbox/include)

add_subdirectory(src)

# -- targets and installation

# The Following is deactivated use: distutils or wheels to install ComPASS
# package install(CODE "message(\"Installing ComPASS python environment.\")")
# make module available through user site packages install(CODE
# "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" configure_user_site.py --add
# ComPASS \"${PROJECT_SOURCE_DIR}\" WORKING_DIRECTORY
# \"${COMPASS_CONFIG_DIRECTORY}\")")

install(EXPORT ComPASS DESTINATION share/cmake/)
# the usual way is to pass the python package directory through
# COMPASS_PYTHON_PACKAGE_DIRECTORY targets (python modules) are also exported
# directly from the build tree to use them you will need to reference the build
# tree cf.
# https://gitlab.kitware.com/cmake/community/wikis/doc/tutorials/Exporting-and-Importing-Targets
export(EXPORT ComPASS FILE share/cmake/compass-exports.cmake)

# The following is an option to revert to old build/install strategy Please be
# aware it will add a .pth file to your local python user site (to know where it
# is type: python -m site --user-site) and you are responsible for this file
# (and possible side effects)
option(ComPASS_USES_OLD_INSTALL
       "Use old manual build/install strategy through cmake." OFF
)
if(ComPASS_USES_OLD_INSTALL)
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -m pip list OUTPUT_VARIABLE PIP_OUTPUT
  )
  string(FIND ${PIP_OUTPUT} ComPASS ComPASS_IS_INSTALLED_WITH_PIP)
  if(${ComPASS_IS_INSTALLED_WITH_PIP} GREATER_EQUAL 0)
    message(
      WARNING
        "\r
      It is possible that ComPASS has already been installed using pip.\r
      Before reverting to old style installation please uninstall ComPASS using pip:\r\r
      ${PYTHON_EXECUTABLE} -m pip uninstall ComPASS\r"
    )
  endif(${ComPASS_IS_INSTALLED_WITH_PIP} GREATER_EQUAL 0)
  set(CMAKE_INSTALL_PREFIX
      ${CMAKE_SOURCE_DIR}
      CACHE PATH "CMake install prefix" FORCE
  )
  install(
    CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" configure_user_site.py --add ComPASS \"${CMAKE_SOURCE_DIR}\" WORKING_DIRECTORY \"${COMPASS_CONFIG_DIRECTORY}\")"
  )
endif(ComPASS_USES_OLD_INSTALL)
