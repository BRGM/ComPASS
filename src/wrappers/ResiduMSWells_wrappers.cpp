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
void ResiduMSWells_associate_pointers(XArrayWrapper<double>,
                                      XArrayWrapper<double>);
void ResiduMSWells_update_accumulation();
// FIXME: all fortran functions could be grouped in a single header file
//        along with some python wrappers for direct binding cf. issue #75
std::size_t number_of_mswell_nodes();
std::size_t number_of_own_mswell_nodes();
int number_of_components();
}

struct ResidualsMSWells {
   std::vector<double> values;
   std::vector<double> accumulations;
   std::size_t total_size;
   auto mswell_nodes_size() const { return total_size; }
   auto mswell_nodes_begin() { return values.data(); }
   auto mswell_nodes_end() { return values.data() + mswell_nodes_size(); }
   auto mswell_nodes_accumulation_size() const { return mswell_nodes_size(); }
   auto mswell_nodes_accumulation_begin() { return accumulations.data(); }
   auto mswell_nodes_accumulation_end() {
      return accumulations.data() + mswell_nodes_size();
   }
   auto own_mswell_nodes_size() const {
      assert(number_of_own_mswell_nodes() >= 0);
      assert(number_of_own_mswell_nodes() <= mswell_nodes_size());
      return number_of_own_mswell_nodes();
   }

   ResidualsMSWells() { reset(); }
   void reset() {
      const auto sou = npv();
      total_size = sou * number_of_mswell_nodes();
      const auto n = total_size;
      if (values.size() != n) {
         values.resize(n);
         values.shrink_to_fit();
         accumulations.resize(n);
         accumulations.shrink_to_fit();
      }
      const auto p = values.data();
      const auto ap = accumulations.data();
      ResiduMSWells_associate_pointers(
          {mswell_nodes_begin(), mswell_nodes_size()},
          {mswell_nodes_accumulation_begin(),
           mswell_nodes_accumulation_size()});
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

void add_ResiduMSWells_wrappers(py::module& module) {
   py::class_<ResidualsMSWells>(module, "ResidualsMSWells")
       .def(py::init())
       .def("reset", &ResidualsMSWells::reset)
       .def_static("npv", &ResidualsMSWells::npv)
       .def_property_readonly(
           "mswell_nodes",
           [](ResidualsMSWells& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.mswell_nodes_size(),
                                        ResidualsMSWells::npv()),
                  self.mswell_nodes_begin()};
           },
           py::keep_alive<0, 1>())
       .def_property_readonly(
           "own_mswell_nodes",
           [](ResidualsMSWells& self) {
              return py::array_t<double>{
                  cast_to_pyarray_shape(self.own_mswell_nodes_size(),
                                        ResidualsMSWells::npv()),
                  self.mswell_nodes_begin()};
           },
           py::keep_alive<0, 1>());
   module.def("ResiduMSWells_update_accumulation",
              &ResiduMSWells_update_accumulation);
}
