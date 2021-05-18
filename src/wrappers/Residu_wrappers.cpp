//
// This file is part of ComPASS.
//
// ComPASS is free software: you can redistribute it and/or modify it under both
// the terms of the GNU General Public License version 3
// (https://www.gnu.org/licenses/gpl.html), and the CeCILL License Agreement
// version 2.1 (http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html).
//

#include <cassert>
#include <limits>
#include <vector>

#include "XArrayWrapper.h"

// Fortran functions
extern "C" {
void Residu_associate_pointers(XArrayWrapper<double>, XArrayWrapper<double>,
                               XArrayWrapper<double>, XArrayWrapper<double>,
                               XArrayWrapper<double>, XArrayWrapper<double>,
                               XArrayWrapper<double>, XArrayWrapper<double>);
void Residu_update_accumulation();
// FIXME: all fortran functions could be grouped in a single header file
//        along with some python wrappers for direct binding cf. issue #75
std::size_t number_of_nodes();
std::size_t number_of_cells();
std::size_t number_of_fractures();
std::size_t number_of_own_nodes();
std::size_t number_of_own_cells();
std::size_t number_of_own_fractures();
std::size_t nb_injectors();
std::size_t number_of_own_injectors();
std::size_t nb_producers();
std::size_t number_of_own_producers();
int number_of_components();
}

struct Residuals {
   std::vector<double> values;
   std::vector<double> accumulations;
   std::size_t fractures_offset;
   std::size_t cells_offset;
   std::size_t injectors_offset;
   std::size_t producers_offset;
   auto nodes_size() const { return fractures_offset; }
   auto nodes_begin() { return values.data(); }
   auto nodes_end() { return values.data() + nodes_size(); }
   auto fractures_size() const { return cells_offset - fractures_offset; }
   auto fractures_begin() { return values.data() + fractures_offset; }
   auto fractures_end() { return values.data() + fractures_size(); }
   auto cells_size() const { return injectors_offset - cells_offset; }
   auto cells_begin() { return values.data() + cells_offset; }
   auto cells_end() { return values.data() + cells_size(); }
   auto injectors_size() const { return producers_offset - injectors_offset; }
   auto injectors_begin() { return values.data() + injectors_offset; }
   auto injectors_end() { return values.data() + injectors_size(); }
   auto producers_size() const { return values.size() - producers_offset; }
   auto producers_begin() { return values.data() + producers_offset; }
   auto producers_end() { return values.data() + producers_size(); }
   auto nodes_accumulation_size() const { return nodes_size(); }
   auto nodes_accumulation_begin() { return accumulations.data(); }
   auto nodes_accumulation_end() {
      return accumulations.data() + fractures_offset;
   }
   auto fractures_accumulation_size() const { return fractures_size(); }
   auto fractures_accumulation_begin() {
      return accumulations.data() + nodes_accumulation_size();
   }
   auto fractures_accumulation_end() {
      return accumulations.data() + fractures_accumulation_size();
   }
   auto cells_accumulation_size() const { return cells_size(); }
   auto cells_accumulation_begin() {
      return accumulations.data() + cells_offset;
   }
   auto cells_accumulation_end() {
      return accumulations.data() + cells_accumulation_size();
   }
   auto own_nodes_size() const {
      assert(number_of_own_nodes() >= 0);
      assert(number_of_own_nodes() <= nodes_size());
      return number_of_own_nodes();
   }
   auto own_fractures_size() const {
      assert(number_of_own_fractures() >= 0);
      assert(number_of_own_fractures() <= fractures_size());
      return number_of_own_fractures();
   }
   auto own_cells_size() const {
      assert(number_of_own_cells() >= 0);
      assert(number_of_own_cells() <= cells_size());
      return number_of_own_cells();
   }
   auto own_injectors_size() const {
      assert(number_of_own_injectors() >= 0);
      assert(number_of_own_injectors() <= injectors_size());
      return number_of_own_injectors();
   }
   auto own_producers_size() const {
      assert(number_of_own_producers() >= 0);
      assert(number_of_own_producers() <= producers_size());
      return number_of_own_producers();
   }
   Residuals() { reset(); }
   void reset() {
      const auto sou = npv();
      fractures_offset = sou * number_of_nodes();
      cells_offset = fractures_offset + sou * number_of_fractures();
      injectors_offset = cells_offset + sou * number_of_cells();
      producers_offset = injectors_offset + nb_injectors();
      const auto n = producers_offset + nb_producers();
      if (values.size() != n) {
         values.resize(n);
         values.shrink_to_fit();
         accumulations.resize(injectors_offset);
         accumulations.shrink_to_fit();
      }
      const auto p = values.data();
      const auto ap = accumulations.data();
      Residu_associate_pointers(
          {nodes_begin(), nodes_size()}, {fractures_begin(), fractures_size()},
          {cells_begin(), cells_size()}, {injectors_begin(), injectors_size()},
          {producers_begin(), producers_size()},
          {nodes_accumulation_begin(), nodes_accumulation_size()},
          {fractures_accumulation_begin(), fractures_accumulation_size()},
          {cells_accumulation_begin(), cells_accumulation_size()});
   }
   // number of primary variables
   static std::size_t npv() {
      assert(number_of_components() == ComPASS_NUMBER_OF_COMPONENTS);
      assert(number_of_components() > 0);
#ifdef _THERMIQUE_
      return static_cast<std::size_t>(number_of_components() + 1);
#else
      return static_cast<std::size_t>(number_of_components());
#endif
   }
};

#include <pybind11/numpy.h>

#include "Residu_wrappers.h"

namespace {
template <typename N>
bool is_valid_pysize(const N n) {
   return (n >= 0) && (n < std::numeric_limits<py::ssize_t>::max());
}

template <typename N>
auto cast_to_pyarray_shape(const N n) {
   assert(is_valid_pysize(n));
   return static_cast<py::ssize_t>(n);
}

template <typename N1, typename N2>
auto cast_to_pyarray_shape(const N1 n1, const N2 n2) {
   assert(is_valid_pysize(n1));
   assert(is_valid_pysize(n2));
   return std::array<py::ssize_t, 2>{
       {static_cast<py::ssize_t>(n1), static_cast<py::ssize_t>(n2)}};
}

}  // namespace

void add_Residu_wrappers(py::module& module) {
   py::class_<Residuals>(module, "Residuals")
       .def(py::init())
       .def("reset", &Residuals::reset)
       .def_static("npv", &Residuals::npv)
       .def_property_readonly(
           "nodes",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.nodes_size(), Residuals::npv()),
                  self.nodes_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_nodes",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_nodes_size(),
                                        Residuals::npv()),
                  self.nodes_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "fractures",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.fractures_size(),
                                        Residuals::npv()),
                  self.fractures_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_fractures",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_fractures_size(),
                                        Residuals::npv()),
                  self.fractures_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "cells",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.cells_size(), Residuals::npv()),
                  self.cells_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_cells",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_cells_size(),
                                        Residuals::npv()),
                  self.cells_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "nodes_accumulation",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.nodes_size(), Residuals::npv()),
                  self.nodes_accumulation_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "fractures_accumulation",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.fractures_size(),
                                        Residuals::npv()),
                  self.fractures_accumulation_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "cells_accumulation",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.cells_size(), Residuals::npv()),
                  self.cells_accumulation_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_injectors",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_injectors_size()),
                  self.injectors_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_producers",
           [](Residuals& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_producers_size()),
                  self.producers_begin()};
           },
           py::keep_alive<0, 1>());

   module.def("Residu_update_accumulation", &Residu_update_accumulation);
}
