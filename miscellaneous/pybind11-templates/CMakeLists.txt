cmake_minimum_required(VERSION 2.8.12)

project(simplecall Fortran CXX)

find_package(pybind11 REQUIRED)
message(STATUS "Message: " ${pybind11_FOUND})

add_library(simplecall_fortranlib simplecall.F90)
pybind11_add_module(simplecall simplecall.cpp simplecall.F90)
target_link_libraries(simplecall simplecall_fortranlib)

