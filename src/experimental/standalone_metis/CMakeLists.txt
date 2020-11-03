cmake_minimum_required(VERSION 3.14)
cmake_policy(VERSION 3.14)

project(standalone-metis C CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE
  )
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "../../../cmake/modules")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 QUIET) # try locally
if(NOT pybind11_FOUND) # else fetch it from github
  message(STATUS "Could NOT find pybind11 (fetching it from GitHub)")
  include(FetchContent)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.4.3
  )
  FetchContent_MakeAvailable(pybind11)
endif()

set(PYBIND11_PYTHON_VERSION
    3
    CACHE STRING "Python version to use for compiling modules"
)

find_package(METIS REQUIRED)

pybind11_add_module(metis metis.cpp)
target_link_libraries(metis PUBLIC ${METIS_LIBRARY})
target_include_directories(metis PUBLIC ${METIS_INCLUDE_DIR})
