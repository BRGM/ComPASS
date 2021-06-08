//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include "LoisThermoHydro.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cassert>
#include <limits>

#include "Simulation.h"
namespace py = pybind11;

py::object py_fill_kr_arrays = py::none{};
py::object py_fill_phase_pressure_arrays = py::none{};

void fill_kr_arrays(const std::size_t n, const int np, X* p_states,
                    int* p_rocktypes, double* p_kr, double* p_dkrdS) {
   assert(n < std::numeric_limits<py::ssize_t>::max());
   assert(X::Model::np < std::numeric_limits<py::ssize_t>::max());
   StateArray states{p_states, n};
   auto simulation = ComPASS::get_simulation();
   py::array_t<int, py::array::c_style> rocktypes{static_cast<py::ssize_t>(n),
                                                  p_rocktypes, simulation};
   py::array_t<double, py::array::c_style> kr{
       {static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(X::Model::np)},
       p_kr,
       simulation};
   py::array_t<double, py::array::c_style> dkrdS{
       {static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(X::Model::np),
        static_cast<py::ssize_t>(X::Model::np)},
       p_dkrdS,
       simulation};
   py_fill_kr_arrays(states, rocktypes, kr, dkrdS);
}

void fill_phase_pressure_arrays(const std::size_t n, const int np, X* p_states,
                                int* p_rocktypes, double* p_dpadS) {
   assert(n < std::numeric_limits<py::ssize_t>::max());
   assert(X::Model::np < std::numeric_limits<py::ssize_t>::max());
   StateArray states{p_states, n};
   auto simulation = ComPASS::get_simulation();
   py::array_t<int, py::array::c_style> rocktypes{static_cast<py::ssize_t>(n),
                                                  p_rocktypes, simulation};
   py::array_t<double, py::array::c_style> dpadS{
       {static_cast<py::ssize_t>(n), static_cast<py::ssize_t>(X::Model::np)},
       p_dpadS,
       simulation};
   py_fill_phase_pressure_arrays(states, rocktypes, dpadS);
}

void update_phase_pressures(X& x) {
   using Phase_vector = X::Model::Phase_vector;
   int rt[1] = {1};
   Phase_vector _;  // dummy vector to store derivatives
   fill_phase_pressure_arrays(1, X::Model::np, &x, rt, _.data());
}

void add_petrophysics_pyinternals(py::module& module) {
   module.def("set_fill_kr_arrays",
              [](py::function f) { py_fill_kr_arrays = f; });
   module.def("get_fill_kr_arrays", []() { return py_fill_kr_arrays; });
   module.def("set_fill_phase_pressure_arrays",
              [](py::function f) { py_fill_phase_pressure_arrays = f; });
   module.def("get_fill_phase_pressure_arrays",
              []() { return py_fill_phase_pressure_arrays; });
}
