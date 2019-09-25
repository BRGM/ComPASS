#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename PetscObject>
auto cast_to_PETSc(py::object obj) {
    // cython object has a handle attribute that stores the PETSc object adress
    auto handle = py::cast<long>(obj.attr("handle"));
    return reinterpret_cast<PetscObject>((void *) handle);
}
