cmake_minimum_required(VERSION 3.8)

project(ComPASS-experimental CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

foreach(file Pack Template_utils)
    add_executable(test-${file} test-${file}.cpp)
endforeach()
