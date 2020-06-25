cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)

project(ComPASS-experimental Fortran CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE
        STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
        FORCE
  )
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(pack)
add_subdirectory(templates)
add_subdirectory(sites)
# add_subdirectory(petsc_from_python)
