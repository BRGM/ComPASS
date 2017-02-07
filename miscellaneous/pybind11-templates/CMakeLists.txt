cmake_minimum_required(VERSION 2.8.12)

project(bindings Fortran CXX)

find_package(pybind11 REQUIRED)

add_library(simplecall_fortranlib simplecall.F90)
pybind11_add_module(simplecall simplecall.cpp)
target_link_libraries(simplecall simplecall_fortranlib)

add_library(passingstring_fortranlib passingstring.F90)
pybind11_add_module(passingstring passingstring.cpp)
target_link_libraries(passingstring passingstring_fortranlib)

set_target_properties(simplecall passingstring
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
    LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}
    LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}
)
