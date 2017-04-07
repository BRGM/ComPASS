#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void add_NN_wrappers(py::module& module);
