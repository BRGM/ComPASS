#pragma once

#include <pybind11/numpy.h>

namespace py = pybind11;

template <typename X, typename ArrayType, typename AttributeType,
          typename PyClass>
auto add_attribute_array(
    PyClass& states, const char* name, std::size_t offset,
    std::vector<std::size_t> shape = std::vector<std::size_t>{},
    std::vector<std::size_t> strides = std::vector<std::size_t>{}) {
   states.def_property_readonly(
       name,
       [=](py::object& self) {
          auto wrapper = self.cast<ArrayType*>();  // C++ pointer to underlying
                                                   // holder array instance
          auto attribute_position = reinterpret_cast<const AttributeType*>(
              reinterpret_cast<const unsigned char*>(wrapper->pointer) +
              offset);
          auto final_shape = std::vector<std::size_t>{{wrapper->length}};
          std::copy(shape.begin(), shape.end(),
                    std::back_inserter(final_shape));
          auto final_strides = std::vector<std::size_t>(1, sizeof(X));
          std::copy(strides.begin(), strides.end(),
                    std::back_inserter(final_strides));
          return py::array_t<AttributeType, py::array::c_style>{
              final_shape, final_strides, attribute_position, self};
       },
       py::return_value_policy::reference_internal  // py::keep_alive<0, 1>():
                                                    // because the StateArray
                                                    // instance must be kept
                                                    // alive as long as the
                                                    // attribute is used -
                                                    // should be ok this is the
                                                    // default for def_property
   );
};
