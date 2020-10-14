//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "Simulation.h"
#include "StateObjects.h"
namespace py = pybind11;

py::function py_fill_kr_arrays;

extern "C" {

void fill_kr_arrays(const std::size_t n, const int np, X* p_states,
                    int* p_rocktypes, double* p_kr, double* p_dkrdS) {
   StateArray states{reinterpret_cast<X*>(p_states), n};
   auto simulation = ComPASS::get_simulation();
   py::array_t<int, py::array::c_style> rocktypes{n, p_rocktypes, simulation};
   py::array_t<double, py::array::c_style> kr{
       {n, X::Model::np}, p_kr, simulation};
   py::array_t<double, py::array::c_style> dkrdS{
       {n, X::Model::np, X::Model::np}, p_dkrdS, simulation};
   py_fill_kr_arrays(states, rocktypes, kr, dkrdS);
}
}

void add_kr_pyinternals(py::module& module) {
   module.def("set_fill_kr_arrays",
              [](py::function& f) { py_fill_kr_arrays = f; });
}
