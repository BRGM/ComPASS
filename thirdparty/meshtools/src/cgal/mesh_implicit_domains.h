#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

py::list mesh_implicit_domains_boundaries();
