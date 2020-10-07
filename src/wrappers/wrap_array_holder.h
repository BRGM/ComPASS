#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

template <typename Array_type>
auto wrap_array_holder(py::module& module, const char* name)
    -> py::class_<Array_type> {
   auto pyclass = py::class_<Array_type>(module, "States");
   pyclass.def("size", [](const Array_type& self) { return self.length; });
   pyclass.def_property_readonly("shape", [](const Array_type& self) {
      return py::make_tuple(self.length);
   });
   pyclass.def("__getitem__", &Array_type::operator[],
               py::return_value_policy::reference);
   pyclass.def("__iter__", [](Array_type& self) {
      return py::make_iterator(self.pointer, self.pointer + self.length);
   });
   pyclass.def("fill", &Array_type::fill);
   pyclass.def("set", &Array_type::fill);
   pyclass.def(
       "set", [](Array_type& self, py::array_t<bool, py::array::c_style>& where,
                 const typename Array_type::wrapped_type& x) {
          assert(where.ndim() == 1);
          assert(self.length == where.size());
          auto use = where.data(0);
          for (auto& value : self) {
             if (*use) value = x;
             ++use;
          }
       });
   pyclass.def("copy", &Array_type::copy);
   return pyclass;
}
